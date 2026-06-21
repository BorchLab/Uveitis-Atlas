# R/11_annotate_celltypes.R
# ------------------------------------------------------------------------------
# Multi-modal cell type annotation using weighted consensus of:
#   (A) Reference labels  (Azimuth, SingleR-HPCA, SingleR-Monaco)  – weight 1 ea.
#   (B) FindAllMarkers gene overlap with canonical marker sets     – weight 2
#   (C) Cluster-level average expression of key discriminating genes – weight 3
#   (D) VDJ lineage detection rate (TCR/BCR)                       – weight 5
#
# Key design decisions:
#   - CD8A/CD8B (not CD4) used to distinguish CD8 vs CD4 T cells, because CD4
#     transcript has poor detection / high dropout in scRNA-seq.
#   - CD3D/CD3E used as pan-T cell markers.
#   - CD19/MS4A1/CD79A used for B cells.
#   - pDC (plasmacytoid dendritic cell) is a dendritic cell, NOT a plasma cell.
#     "Plasma" is reserved for antibody-secreting cells (plasmablasts).
#   - All reference labels are normalised to Azimuth-style canonical names
#     before consensus voting, so that e.g. "T_cell:CD4+_effector_memory"
#     (HPCA) and "Effector memory CD4 T cells" (Monaco) both map to "CD4 TEM".
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(leidenAlg)
})

# ==============================================================================
# 1. CANONICAL MARKERS – used for FindAllMarkers overlap AND expression scoring
# ==============================================================================
CANONICAL_MARKERS <- list(
  "T cell"              = c("CD3D", "CD3E", "CD3G", "TRAC"),
  "CD8 T cell"          = c("CD8A", "CD8B"),
  "gdT"                 = c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2"),
  "NK"                  = c("NKG7", "GNLY", "KLRD1", "NCAM1"),
  "B cell"              = c("CD19", "MS4A1", "CD79A", "CD79B"),
  "Plasma"              = c("MZB1", "JCHAIN", "XBP1", "SDC1"),
  "Monocyte/Macrophage" = c("CD14", "LYZ", "S100A8", "S100A9", "CD68"),
  "cDC"                 = c("FCER1A", "CLEC10A", "CD1C"),
  "pDC"                 = c("LILRA4", "IRF7", "IL3RA", "CLEC4C"),
  "Platelet"            = c("PPBP", "PF4", "GP9"),
  "RBC"                 = c("HBA1", "HBA2", "HBB"),
  "Neutrophil"          = c("CSF3R", "FCGR3B", "CXCR2", "G0S2"),
  "Mast cell"           = c("TPSAB1", "TPSB2", "KIT", "CPA3"),
  "Basophil"            = c("CLC", "HDC", "GATA2"),
  "Proliferating"       = c("MKI67", "TOP2A", "STMN1")
)

# ==============================================================================
# 2. LABEL NORMALISATION DICTIONARY
#    Maps reference-specific labels → canonical Azimuth-style names.
#    Applied BEFORE consensus voting so that equivalent labels from different
#    references are counted as the same cell type.
# ==============================================================================
LABEL_NORMALIZE <- c(
  # ---- HPCA / SingleR labels ----
  "T_cell:CD4+_central_memory"    = "CD4 TCM",
  "T_cell:CD4+_effector_memory"   = "CD4 TEM",
  "T_cell:CD8+_central_memory"    = "CD8 TCM",
  "T_cell:CD8+_effector_memory"   = "CD8 TEM",
  "T_cell:CD4+_naive"             = "CD4 Naive",
  "T_cell:CD8+_naive"             = "CD8 Naive",
  "T_cell:CD4+_Th1"               = "CD4 TEM",
  "T_cell:CD4+_Th2"               = "CD4 TEM",
  "T_cell:CD4+_Th17"              = "CD4 TEM",
  "T_cell:gamma-delta"            = "gdT",
  "T_cell:CD4+_regulatory"        = "Treg",
  "T_cells:Treg"                  = "Treg",
  "NK_cell"                       = "NK",
  "Pre-B_cell_CD34-"              = "B naive",
  "B_cell"                        = "B memory",
  "Macrophage"                    = "Macrophage",
  "Monocyte"                      = "CD14 Mono",

  # ---- Monaco labels ----
  "Central memory CD4 T cells"    = "CD4 TCM",
  "Effector memory CD4 T cells"   = "CD4 TEM",
  "Central memory CD8 T cells"    = "CD8 TCM",
  "Effector memory CD8 T cells"   = "CD8 TEM",
  "Naive CD4 T cells"             = "CD4 Naive",
  "Naive CD8 T cells"             = "CD8 Naive",
  "Terminal effector CD4 T cells"  = "CD4 CTL",
  "Terminal effector CD8 T cells"  = "CD8 TEM",
  "T regulatory cells"            = "Treg",
  "Follicular helper T cells"     = "Tfh",
  "Vd2 gd T cells"                = "gdT",
  "Non-Vd2 gd T cells"            = "gdT",
  "MAIT cells"                    = "MAIT",
  "Natural killer cells"          = "NK",

  "Classical monocytes"           = "CD14 Mono",
  "Non classical monocytes"       = "CD16 Mono",
  "Non-classical monocytes"       = "CD16 Mono",
  "Intermediate monocytes"        = "CD14 Mono",

  "Naive B cells"                 = "B naive",
  "Non-switched memory B cells"   = "B memory",
  "Switched memory B cells"       = "B memory",
  "Exhausted B cells"             = "B intermediate",

  "Plasmacytoid dendritic cells"  = "pDC",
  "Plasmacytoid dendritic cell"   = "pDC",
  "Myeloid dendritic cell"        = "cDC2",
  "Myeloid dendritic cells"       = "cDC2",

  "Plasma"                        = "Plasma",
  "Plasmablasts"                  = "Plasma",
  "Plasma cell"                   = "Plasma",
  "Plasma cells"                  = "Plasma",
  "Low-density basophils"         = "Basophil"
)

# ==============================================================================
# 3. BROAD MAPPING
#    NOTE: pDC lives under "DC", NOT "Plasma" or "Monocyte/Macrophage".
#    "Plasma" is strictly for antibody-secreting cells.
# ==============================================================================
BROAD_MAPPING <- list(
  "T cell" = c(
    "CD4 T", "CD8 T", "T cell", "Treg", "gdT", "dnT", "MAIT", "Tfh",
    "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating",
    "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating",
    "T Proliferating"
  ),
  "NK" = c(
    "NK", "NK_CD56bright", "NK Proliferating"
  ),
  "B cell" = c(
    "B cell", "B naive", "B intermediate", "B memory"
  ),
  "Plasma" = c(
    "Plasma", "Plasma cell"
  ),
  "Monocyte/Macrophage" = c(
    "CD14 Mono", "CD16 Mono", "Macrophage"
  ),
  "DC" = c(
    "cDC", "cDC1", "cDC2", "pDC", "ASDC"
  ),
  "Granulocyte" = c(
    "Neutrophil", "Mast cell", "Basophil", "Eosinophil"
  ),
  "Platelet" = c(
    "Platelet"
  ),
  "Eryth" = c(
    "Eryth"
  ),
  "Other" = c(
    "HSPC", "ILC", "Doublet"
  )
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Normalise a reference label to canonical form
normalize_label <- function(label) {
  if (is.na(label) || label == "") return(NA_character_)
  if (label %in% names(LABEL_NORMALIZE)) return(unname(LABEL_NORMALIZE[label]))
  label
}

#' Map a (normalised) fine label to a broad category
#' Order of fallback checks matters: "plasmacytoid" BEFORE "plasma"
map_to_broad <- function(label) {
  if (is.na(label) || label == "") return("Other")
  normed <- normalize_label(label)
  lbl <- if (!is.na(normed)) normed else label

  # Exact match in BROAD_MAPPING values
  for (broad in names(BROAD_MAPPING)) {
    if (lbl %in% BROAD_MAPPING[[broad]]) return(broad)
  }
  # Case-insensitive exact match
  lbl_lower <- tolower(lbl)
  for (broad in names(BROAD_MAPPING)) {
    if (lbl_lower %in% tolower(BROAD_MAPPING[[broad]])) return(broad)
  }

  # Keyword fallback – SPECIFIC patterns before GENERAL ones
  # CRITICAL: "plasmacytoid" must be checked BEFORE "plasma"
  if (grepl("plasmacytoid|\\bpdc\\b", lbl_lower))                            return("DC")
  if (grepl("dendritic|\\bcdc|\\basdc", lbl_lower))                          return("DC")
  if (grepl("\\bcd4|\\bcd8|\\bt cell|\\btreg|\\bmait|\\bgdt|\\bdnt", lbl_lower)) return("T cell")
  if (grepl("\\bnk\\b|natural killer", lbl_lower))                           return("NK")
  if (grepl("\\bb cell|\\bb naive|\\bb mem", lbl_lower))                     return("B cell")
  if (grepl("\\bplasma\\b|plasmablast", lbl_lower))                          return("Plasma")
  if (grepl("mono|macro|\\bcd14\\b|\\bcd16\\b", lbl_lower))                  return("Monocyte/Macrophage")
  if (grepl("neutro|granulo|\\bmast\\b|basoph|eosinoph", lbl_lower))         return("Granulocyte")
  if (grepl("platelet|megakaryocyte", lbl_lower))                            return("Platelet")
  if (grepl("eryth|\\brbc\\b|erythro", lbl_lower))                          return("Eryth")
  "Other"
}

#' Get mode (most frequent value) of a vector
get_mode <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

#' VDJ detection rate per cluster, with TCR vs BCR chain breakdown.
#' Uses CTgene column (concatenated V genes from scRepertoire) to
#' distinguish TCR (TRAV/TRBV/TRDV/TRGV) from BCR (IGHV/IGLV/IGKV).
get_vdj_metrics <- function(meta, cluster_col) {
  has_ctaa  <- "CTaa"  %in% colnames(meta)
  has_ctgene <- "CTgene" %in% colnames(meta)

  if (!has_ctaa && !has_ctgene) {
    out <- data.frame(cl = character(), n = integer(), v = numeric(),
                      tcr = numeric(), bcr = numeric(),
                      stringsAsFactors = FALSE)
    colnames(out) <- c(cluster_col, "n_cells", "vdj_prop", "tcr_prop", "bcr_prop")
    return(out)
  }

  meta %>%
    group_by(!!sym(cluster_col)) %>%
    summarise(
      n_cells  = n(),
      vdj_prop = if (has_ctaa) sum(!is.na(CTaa) & CTaa != "") / n() else 0,
      tcr_prop = if (has_ctgene) sum(grepl("TR[ABDG]V", CTgene, perl = TRUE), na.rm = TRUE) / n() else NA_real_,
      bcr_prop = if (has_ctgene) sum(grepl("IG[HKL]V", CTgene, perl = TRUE), na.rm = TRUE) / n() else NA_real_,
      .groups  = "drop"
    )
}

# ==============================================================================
# EXPRESSION-BASED EVIDENCE
# Compute cluster-level average expression of key discriminating markers.
# A gene is "enriched" in a cluster when its mean expression exceeds the
# cross-cluster median for that gene (relative enrichment).
# ==============================================================================
compute_expression_evidence <- function(obj, cluster_col) {
  key_genes <- unique(unlist(CANONICAL_MARKERS))
  key_genes <- intersect(key_genes, rownames(obj))
  if (length(key_genes) == 0) return(NULL)

  avg <- AverageExpression(obj, features = key_genes,
                           group.by = cluster_col, assays = "RNA")$RNA

  # Seurat v5 may prefix column names (e.g. "RNA_0" instead of "0").
  # Strip common prefixes so column names match bare cluster IDs.
  known_clusters <- as.character(sort(unique(obj[[cluster_col, drop = TRUE]])))
  cnames <- colnames(avg)
  if (!all(known_clusters %in% cnames)) {
    # Try "RNA_0" → "0" (underscore-separated prefix)
    stripped <- sub("^[A-Za-z]+_", "", cnames)
    if (all(known_clusters %in% stripped) && !anyDuplicated(stripped)) {
      colnames(avg) <- stripped
    } else {
      # Try "g0" → "0" (Seurat v5 prepends "g" when cluster IDs start with a digit)
      stripped2 <- sub("^g", "", cnames)
      if (all(known_clusters %in% stripped2) && !anyDuplicated(stripped2)) {
        colnames(avg) <- stripped2
      }
    }
  }

  gene_medians <- apply(avg, 1, median)
  enriched <- sweep(avg, 1, gene_medians, ">")  # genes × clusters (logical)

  list(avg = avg, enriched = enriched, genes = key_genes)
}

#' Score broad cell type from expression evidence for one cluster.
#' Uses biologically motivated rules:
#'   CD3D/CD3E  → T cell (pan-T)
#'   NKG7/GNLY without CD3D → NK
#'   CD19/MS4A1/CD79A (≥2) → B cell
#'   MZB1 without CD19 → Plasma (plasma cells lose surface CD19)
#'   LILRA4/IRF7 → pDC → DC category
#'   FCER1A/CLEC10A → cDC → DC category
#'   CD14 or LYZ+S100A8 → Monocyte/Macrophage
score_cluster_expression <- function(expr_data, cl) {
  broad_types <- c("T cell", "NK", "B cell", "Plasma",
                   "Monocyte/Macrophage", "DC", "Granulocyte",
                   "Platelet", "Eryth", "Other")
  scores <- setNames(rep(0, length(broad_types)), broad_types)
  if (is.null(expr_data)) return(scores)

  cl_char <- as.character(cl)
  if (!(cl_char %in% colnames(expr_data$enriched))) return(scores)

  enr <- function(g) {
    if (g %in% rownames(expr_data$enriched)) expr_data$enriched[g, cl_char] else FALSE
  }

  # T cell: CD3D or CD3E enriched
  if (enr("CD3D") || enr("CD3E")) scores["T cell"] <- 3

  # gdT: TRDC enriched → T cell signal, but only if CD3 is also present
  # (some NK subsets express delta chain genes without CD3)
  if ((enr("TRDC") || enr("TRDV1") || enr("TRDV2")) && (enr("CD3D") || enr("CD3E"))) {
    scores["T cell"] <- max(scores["T cell"], 3)
  }


  # NK: NKG7/GNLY with CD3 truly absent (raw avg < 1.0)
  # Enrichment-based CD3 check is unreliable — CD8 TEM clusters can have
  # CD3D below median if outnumbered by CD4 clusters, falsely triggering NK.
  nk_cytotoxic <- enr("NKG7") || enr("GNLY") || enr("KLRD1")
  cd3_abs <- 0
  if (cl_char %in% colnames(expr_data$avg)) {
    cd3_abs <- max(
      if ("CD3D" %in% rownames(expr_data$avg)) expr_data$avg["CD3D", cl_char] else 0,
      if ("CD3E" %in% rownames(expr_data$avg)) expr_data$avg["CD3E", cl_char] else 0
    )
  }
  if (nk_cytotoxic && cd3_abs < 1.0) scores["NK"] <- 3
  # FCGR3A without myeloid markers (CD14, LYZ) → likely NK, not CD16 Mono
  if (enr("FCGR3A") && !enr("CD14") && !enr("LYZ") && !enr("CD3D")) {
    scores["NK"] <- max(scores["NK"], 2)
  }

  # B cell: ≥2 of CD19, MS4A1, CD79A (enrichment OR raw avg > 0.3)
  b_enr_count <- sum(c(enr("CD19"), enr("MS4A1"), enr("CD79A")))
  if (b_enr_count >= 2) {
    scores["B cell"] <- 3
  } else if (cl_char %in% colnames(expr_data$avg)) {
    .avg_b <- function(g) if (g %in% rownames(expr_data$avg)) expr_data$avg[g, cl_char] else 0
    b_avg_count <- sum(c(.avg_b("CD19") > 0.3, .avg_b("MS4A1") > 0.3, .avg_b("CD79A") > 0.3))
    if (b_avg_count >= 2) scores["B cell"] <- 3
  }

  # Plasma: MZB1+ but CD19- (plasma cells downregulate CD19/CD20)
  if (enr("MZB1") && !enr("CD19")) scores["Plasma"] <- 3

  # pDC → DC (LILRA4/IRF7, NOT CD14, NOT CD3D)
  if ((enr("LILRA4") || enr("IRF7")) && !enr("CD14") && !enr("CD3D")) {
    scores["DC"] <- 3
  }

  # cDC → DC (FCER1A/CLEC10A, NOT CD14, NOT CD3D)
  if ((enr("FCER1A") || enr("CLEC10A")) && !enr("CD14") && !enr("CD3D")) {
    scores["DC"] <- max(scores["DC"], 3)
  }

  # Monocyte: CD14+ or (LYZ+ & S100A8+) — requires true myeloid markers
  if (enr("CD14") || (enr("LYZ") && enr("S100A8"))) {
    scores["Monocyte/Macrophage"] <- 3
  }

  # CD16 Mono: FCGR3A+ with at least one myeloid marker (CD14 or LYZ),
  # without T/B/NK cytotoxic markers — prevents NK→Mono misclassification
  if (enr("FCGR3A") && (enr("LYZ") || enr("S100A8")) &&
      !enr("CD3D") && !enr("CD19") && !enr("NKG7")) {
    scores["Monocyte/Macrophage"] <- max(scores["Monocyte/Macrophage"], 2)
  }

  # Granulocyte: neutrophils, mast cells, basophils
  # Neutrophils: CSF3R or FCGR3B without lymphoid markers
  if ((enr("CSF3R") || enr("FCGR3B")) && !enr("CD3D") && !enr("CD19") && !enr("NKG7")) {
    scores["Granulocyte"] <- 3
  }
  # Mast cells: TPSAB1 or TPSB2 or KIT+CPA3
  if (enr("TPSAB1") || enr("TPSB2") || (enr("KIT") && enr("CPA3"))) {
    scores["Granulocyte"] <- max(scores["Granulocyte"], 3)
  }
  # Basophils: CLC or HDC without T/B markers
  if ((enr("CLC") || enr("HDC")) && !enr("CD3D") && !enr("CD19")) {
    scores["Granulocyte"] <- max(scores["Granulocyte"], 2)
  }

  # Platelet: PPBP or PF4, without lymphoid/myeloid markers
  if ((enr("PPBP") || enr("PF4") || enr("GP9")) &&
      !enr("CD3D") && !enr("CD14") && !enr("CD19")) {
    scores["Platelet"] <- 3
  }

  # Erythrocyte: HBA1/HBA2/HBB without other lineage markers
  if ((enr("HBA1") || enr("HBA2") || enr("HBB")) &&
      !enr("CD3D") && !enr("CD14") && !enr("CD19")) {
    scores["Eryth"] <- 3
  }

  # pDC high-confidence: IRF7 is extremely specific for pDC.
  # Raw average > 5 is strong evidence regardless of enrichment threshold,
  # since pDC clusters are small and IRF7 median can be ~0.
  if (cl_char %in% colnames(expr_data$avg)) {
    irf7_avg <- if ("IRF7" %in% rownames(expr_data$avg)) expr_data$avg["IRF7", cl_char] else 0
    if (irf7_avg > 5) {
      scores["DC"] <- max(scores["DC"], 4)  # strong override
    }
  }

  scores
}

#' Within T-cell clusters, determine CD4 vs CD8 from expression.
#' CD8A/CD8B enriched → CD8; otherwise CD4 (since CD4 transcript is
#' poorly detected in scRNA-seq, CD3+/CD8- is the best proxy for CD4).
determine_t_subtype <- function(expr_data, cl) {
  if (is.null(expr_data)) return("CD4")
  cl_char <- as.character(cl)
  if (!(cl_char %in% colnames(expr_data$enriched))) return("CD4")
  enr <- function(g) {
    if (g %in% rownames(expr_data$enriched)) expr_data$enriched[g, cl_char] else FALSE
  }
  if (enr("CD8A") || enr("CD8B")) return("CD8")
  "CD4"
}

# ==============================================================================
# DATA-DRIVEN SNN CLUSTER REFINEMENT
# Evaluates each cluster for heterogeneity using:
#   1. Multi-reference annotation distribution (Azimuth, SingleR-HPCA, Monaco)
#   2. Classical marker expression conflicts
#   3. SNN sub-clustering to split where evidence supports it
#
# Sub-cluster identity is determined by weighted consensus of ALL available
# reference annotations (not just Azimuth), mirroring the main annotation logic.
# ==============================================================================
refine_clusters_by_evidence <- function(obj, annotation_map, snn_graph, expr_data,
                                        cluster_col, ref_cols, cfg) {
  min_maj_frac        <- cfg$min_majority_fraction   %||% 0.70
  min_sec_frac        <- cfg$min_secondary_fraction  %||% 0.20
  ambig_ratio         <- cfg$min_ambiguity_ratio     %||% 2.0
  split_res           <- cfg$split_resolution        %||% 0.5
  max_total           <- cfg$max_total_clusters      %||% 25L
  min_subcluster_size <- cfg$min_subcluster_size     %||% 20L

  # Convert cluster column to character before splitting so new IDs don't
  # produce NA from invalid factor levels
  obj[[cluster_col]] <- as.character(obj[[cluster_col, drop = TRUE]])
  meta <- obj[[]]  # snapshot AFTER character conversion
  clusters <- as.character(annotation_map$cluster)
  avail_refs <- ref_cols[ref_cols %in% colnames(meta)]

  if (length(avail_refs) == 0) {
    log_message("Refinement: no reference annotation columns available, skipping")
    return(annotation_map)
  }
  log_message(sprintf("Refinement using %d reference annotations: %s",
                      length(avail_refs), paste(avail_refs, collapse = ", ")))

  # Pre-fetch normalised expression data matrix for sub-cluster scoring
  counts <- Seurat::GetAssayData(obj, layer = "data")
  key_genes <- intersect(unique(unlist(CANONICAL_MARKERS)), rownames(counts))
  gene_medians <- if (!is.null(expr_data)) {
    apply(expr_data$avg[intersect(key_genes, rownames(expr_data$avg)), , drop = FALSE], 1, median)
  } else {
    setNames(rep(0, length(key_genes)), key_genes)
  }

  new_rows <- list()
  clusters_to_remove <- c()

  # Ensure cluster column is character (not factor) to avoid NA on rbind
  annotation_map$cluster <- as.character(annotation_map$cluster)

  # Track max cluster ID across ALL iterations (not recomputed each time)
  max_cluster_id <- max(as.integer(annotation_map$cluster), na.rm = TRUE)

  for (i in seq_len(nrow(annotation_map))) {
    cl <- annotation_map$cluster[i]
    cl_char <- as.character(cl)
    cell_idx <- which(meta[[cluster_col]] == cl)
    if (length(cell_idx) < 20) next

    # ------------------------------------------------------------------
    # Check 1: Multi-reference annotation heterogeneity
    # Pool normalised labels from ALL available references, then check
    # whether multiple broad lineages are represented.
    # ------------------------------------------------------------------
    needs_split <- FALSE

    # Collect per-cell broad lineage calls from every reference
    ref_broad_calls <- lapply(avail_refs, function(rc) {
      raw_labels <- meta[[rc]][cell_idx]
      raw_labels <- raw_labels[!is.na(raw_labels) & raw_labels != ""]
      sapply(raw_labels, function(x) map_to_broad(normalize_label(x)))
    })
    all_broad_calls <- unlist(ref_broad_calls)

    if (length(all_broad_calls) > 10) {
      broad_tab <- sort(table(all_broad_calls), decreasing = TRUE)
      broad_fracs <- broad_tab / sum(broad_tab)
      if (broad_fracs[1] < min_maj_frac && length(broad_fracs) >= 2 && broad_fracs[2] >= min_sec_frac) {
        needs_split <- TRUE
        log_message(sprintf(
          "  Cluster %s: multi-ref heterogeneity (top=%s %.0f%%, 2nd=%s %.0f%%, n_refs=%d)",
          cl, names(broad_fracs)[1], broad_fracs[1]*100,
          names(broad_fracs)[2], broad_fracs[2]*100, length(avail_refs)))
      }
    }

    # Also check within-reference disagreement (different references calling
    # different broad lineages for the same cluster)
    if (!needs_split && length(avail_refs) >= 2) {
      ref_modes <- sapply(avail_refs, function(rc) {
        raw_labels <- meta[[rc]][cell_idx]
        raw_labels <- raw_labels[!is.na(raw_labels) & raw_labels != ""]
        if (length(raw_labels) == 0) return(NA_character_)
        mode_lbl <- names(sort(table(raw_labels), decreasing = TRUE))[1]
        map_to_broad(normalize_label(mode_lbl))
      })
      ref_modes <- ref_modes[!is.na(ref_modes)]
      if (length(unique(ref_modes)) > 1) {
        needs_split <- TRUE
        log_message(sprintf("  Cluster %s: cross-reference disagreement (%s)",
                            cl, paste(paste0(avail_refs[!is.na(ref_modes)], "=", ref_modes), collapse = ", ")))
      }
    }

    # ------------------------------------------------------------------
    # Check 2: Classical marker conflicts
    # ------------------------------------------------------------------
    if (!needs_split && !is.null(expr_data) && cl_char %in% colnames(expr_data$avg)) {
      avg_expr <- expr_data$avg[, cl_char]
      t_sig  <- max(avg_expr[intersect(c("CD3D", "CD3E"), names(avg_expr))], 0)
      nk_sig <- max(avg_expr[intersect(c("NKG7", "GNLY", "KLRD1"), names(avg_expr))], 0)
      b_sig  <- max(avg_expr[intersect(c("CD19", "MS4A1", "CD79A"), names(avg_expr))], 0)
      mono_sig <- max(avg_expr[intersect(c("CD14", "LYZ", "S100A8"), names(avg_expr))], 0)

      sigs <- c(T_cell = t_sig, NK = nk_sig, B_cell = b_sig, Mono = mono_sig)
      top_sigs <- sort(sigs[sigs > 0.5], decreasing = TRUE)
      if (length(top_sigs) >= 2 && top_sigs[1] / top_sigs[2] < 3) {
        needs_split <- TRUE
        log_message(sprintf("  Cluster %s: marker conflict (%s=%.2f, %s=%.2f)",
                            cl, names(top_sigs)[1], top_sigs[1],
                            names(top_sigs)[2], top_sigs[2]))
      }
    }

    # ------------------------------------------------------------------
    # Check 3: NK/T cell mixing — CD8 TEM clusters often contain NK cells
    # that share NKG7/GNLY expression. Detect by checking if a substantial
    # fraction of cells lack CD3 expression despite cluster-level CD3 avg > 0.
    # ------------------------------------------------------------------
    if (!needs_split && !is.null(counts) && annotation_map$celltype_broad[i] == "T cell") {
      nkg7_genes <- intersect(c("NKG7", "GNLY"), rownames(counts))
      cd3_genes  <- intersect(c("CD3D", "CD3E"), rownames(counts))
      if (length(nkg7_genes) > 0 && length(cd3_genes) > 0 && length(cell_idx) > 50) {
        nkg7_expr <- Matrix::colSums(counts[nkg7_genes, cell_idx, drop = FALSE])
        cd3_expr  <- Matrix::colSums(counts[cd3_genes, cell_idx, drop = FALSE])
        # Fraction of cells that are NKG7+ but CD3-
        nk_like <- sum(nkg7_expr > 0 & cd3_expr == 0) / length(cell_idx)
        if (nk_like > 0.10) {
          needs_split <- TRUE
          log_message(sprintf("  Cluster %s: NK/T mixing (%.0f%% NKG7+/CD3- cells)",
                              cl, nk_like * 100))
        }
      }
    }

    # ------------------------------------------------------------------
    # Check 4: Cluster flagged as Mixed, Unresolved, or Other
    # ------------------------------------------------------------------
    unresolved_types <- c("Mixed", "Other")
    unresolved_labels <- c("Mixed", "Unresolved", "Other", "Unassigned")
    if (!needs_split && (annotation_map$celltype_broad[i] %in% unresolved_types ||
                         annotation_map$celltype[i] %in% unresolved_labels)) {
      needs_split <- TRUE
      log_message(sprintf("  Cluster %s: flagged as %s/%s by consensus — attempting split",
                          cl, annotation_map$celltype[i], annotation_map$celltype_broad[i]))
    }

    if (!needs_split) next

    # Guard: don't exceed max total clusters (count only clusters with cells,
    # not ghost clusters from prior refinement runs that have 0 cells)
    live_clusters <- sum(sapply(annotation_map$cluster, function(c) {
      sum(meta[[cluster_col]] == c) > 0
    }))
    current_total <- live_clusters + length(new_rows) - length(clusters_to_remove)
    if (current_total >= max_total) {
      log_message(sprintf("  Skipping split of cluster %s: would exceed %d total clusters (live=%d)", cl, max_total, live_clusters))
      next
    }

    # ------------------------------------------------------------------
    # SNN sub-clustering
    # ------------------------------------------------------------------
    n_vertices <- igraph::vcount(snn_graph)
    if (max(cell_idx) > n_vertices) {
      log_message(sprintf("  Cluster %s: cell indices exceed graph size, skipping", cl))
      next
    }

    subg <- tryCatch(
      igraph::induced_subgraph(snn_graph, cell_idx),
      error = function(e) NULL
    )
    if (is.null(subg) || igraph::vcount(subg) < 20) next

    sub_memb <- tryCatch(
      leidenAlg::leiden.community(subg, resolution = split_res)$membership,
      error = function(e) NULL
    )
    if (is.null(sub_memb)) next

    sub_clusters <- unique(sub_memb)
    if (length(sub_clusters) < 2) next

    # ------------------------------------------------------------------
    # Absorb sub-clusters below min_subcluster_size into the largest
    # sub-cluster (mirrors Seurat's group.singletons behavior). This
    # prevents Leiden micro-clusters / singletons from either (a) creating
    # tiny display clusters downstream, or (b) forcing us to reject an
    # otherwise-good split because of a handful of outlier cells.
    # ------------------------------------------------------------------
    sub_sizes <- table(sub_memb)
    tiny_sub  <- names(sub_sizes)[sub_sizes < min_subcluster_size]
    if (length(tiny_sub) > 0 && length(tiny_sub) < length(sub_sizes)) {
      largest_sub <- as.integer(names(sub_sizes)[which.max(sub_sizes)])
      n_absorbed  <- sum(sub_memb %in% as.integer(tiny_sub))
      sub_memb[sub_memb %in% as.integer(tiny_sub)] <- largest_sub
      sub_clusters <- unique(sub_memb)
      log_message(sprintf(
        "  Cluster %s: absorbed %d cells from %d sub-cluster(s) below size %d into largest sub-cluster",
        cl, n_absorbed, length(tiny_sub), min_subcluster_size))
    }
    if (length(sub_clusters) < 2) next  # absorption collapsed to one sub-cluster

    # ------------------------------------------------------------------
    # Validate each sub-cluster using multi-reference consensus + expression
    # ------------------------------------------------------------------
    sub_labels <- list()
    valid_split <- TRUE
    for (sc in sub_clusters) {
      sc_cells <- cell_idx[sub_memb == sc]
      # All sub-clusters should be >= min_subcluster_size after absorption.
      # Keep a defensive floor in case min_subcluster_size was set very low.
      if (length(sc_cells) < min_subcluster_size) { valid_split <- FALSE; break }

      sc_scores <- setNames(rep(0, length(names(BROAD_MAPPING))), names(BROAD_MAPPING))

      # (A) Reference votes: weight 1 per reference (same as main annotation)
      ref_fine_pool <- c()
      for (rc in avail_refs) {
        rc_labels <- meta[[rc]][sc_cells]
        rc_labels <- rc_labels[!is.na(rc_labels) & rc_labels != ""]
        if (length(rc_labels) == 0) next
        rc_mode <- names(sort(table(rc_labels), decreasing = TRUE))[1]
        rc_normed <- normalize_label(rc_mode)
        rc_broad <- map_to_broad(rc_normed)
        if (rc_broad %in% names(sc_scores)) sc_scores[rc_broad] <- sc_scores[rc_broad] + 1
        # Collect all normalised fine labels for consensus fine label later
        rc_all_normed <- sapply(rc_labels, normalize_label, USE.NAMES = FALSE)
        ref_fine_pool <- c(ref_fine_pool, rc_all_normed[!is.na(rc_all_normed)])
      }

      # (B) Expression scoring: weight 3 (same as main annotation)
      if (length(key_genes) > 0 && length(sc_cells) > 0) {
        sc_avg <- Matrix::rowMeans(counts[key_genes, sc_cells, drop = FALSE])
        sc_enr <- sc_avg[intersect(names(sc_avg), names(gene_medians))] >
                  gene_medians[intersect(names(sc_avg), names(gene_medians))]

        enr <- function(g) if (g %in% names(sc_enr)) sc_enr[g] else FALSE
        if (enr("CD3D") || enr("CD3E")) sc_scores["T cell"] <- sc_scores["T cell"] + 3
        if ((enr("NKG7") || enr("GNLY")) && !enr("CD3D")) sc_scores["NK"] <- sc_scores["NK"] + 3
        if (sum(c(enr("CD19"), enr("MS4A1"), enr("CD79A"))) >= 2) sc_scores["B cell"] <- sc_scores["B cell"] + 3
        if (enr("MZB1") && !enr("CD19")) sc_scores["Plasma"] <- sc_scores["Plasma"] + 3
        if (enr("CD14") || (enr("LYZ") && enr("S100A8"))) sc_scores["Monocyte/Macrophage"] <- sc_scores["Monocyte/Macrophage"] + 3
        if ((enr("LILRA4") || enr("FCER1A")) && !enr("CD14") && !enr("CD3D")) sc_scores["DC"] <- sc_scores["DC"] + 3
      }

      top_two_sc <- sort(sc_scores, decreasing = TRUE)[1:2]
      if (!is.na(top_two_sc[2]) && top_two_sc[2] > 0 &&
          top_two_sc[1] / max(top_two_sc[2], 1) < ambig_ratio) {
        valid_split <- FALSE
        break
      }

      sc_broad <- names(which.max(sc_scores))

      # (C) Fine label: consensus across ALL references (not just Azimuth)
      # Filter pooled fine labels to those matching the winning broad lineage
      sc_fine <- sc_broad
      if (length(ref_fine_pool) > 0) {
        valid_fine <- ref_fine_pool[sapply(ref_fine_pool, function(x) map_to_broad(x) == sc_broad)]
        if (length(valid_fine) > 0) {
          sc_fine_raw <- names(sort(table(valid_fine), decreasing = TRUE))[1]
          sc_fine_normed <- normalize_label(sc_fine_raw)
          if (!is.na(sc_fine_normed) && sc_fine_normed != "") sc_fine <- sc_fine_normed
        }
      }

      # T cell CD4/CD8 refinement for sub-cluster
      if (sc_broad == "T cell" && !grepl("CD4|CD8|gdT|dnT|MAIT|Treg|Tfh", sc_fine)) {
        if (length(key_genes) > 0 && length(sc_cells) > 0) {
          cd8_avg <- mean(c(
            if ("CD8A" %in% key_genes) sc_avg["CD8A"] else 0,
            if ("CD8B" %in% key_genes) sc_avg["CD8B"] else 0
          ))
          sc_fine <- if (cd8_avg > gene_medians["CD8A"] %||% 0) "CD8 T" else "CD4 T"
        }
      }

      sub_labels[[as.character(sc)]] <- list(
        cells = sc_cells, broad = sc_broad, fine = sc_fine,
        top_score = top_two_sc[1], n_cells = length(sc_cells)
      )
    }

    if (!valid_split) {
      log_message(sprintf("  Cluster %s: split not validated (sub-clusters ambiguous), keeping original", cl))
      next
    }

    # Check that sub-clusters actually have different labels
    sub_fines <- sapply(sub_labels, function(x) x$fine)
    if (length(unique(sub_fines)) < 2) {
      log_message(sprintf("  Cluster %s: all sub-clusters map to same type (%s), keeping original", cl, sub_fines[1]))
      next
    }

    # Accept split: create new annotation rows and update cell labels
    clusters_to_remove <- c(clusters_to_remove, i)

    for (j in seq_along(sub_labels)) {
      sl <- sub_labels[[j]]
      new_id <- max_cluster_id + j
      obj[[cluster_col]][sl$cells, 1] <- as.character(new_id)

      new_rows[[length(new_rows) + 1]] <- data.frame(
        cluster        = as.character(new_id),
        celltype       = sl$fine,
        celltype_broad = sl$broad,
        vdj_confidence = annotation_map$vdj_confidence[i],
        top_score      = sl$top_score,
        n_ref_agree    = NA_integer_,
        n_refs_used    = length(avail_refs),
        stringsAsFactors = FALSE
      )
      log_message(sprintf("  Cluster %s -> new cluster %d: %s (%s, n=%d)",
                          cl, new_id, sl$fine, sl$broad, sl$n_cells))
    }
    max_cluster_id <- max_cluster_id + length(sub_labels)
  }

  # Apply splits
  if (length(clusters_to_remove) > 0) {
    annotation_map <- annotation_map[-clusters_to_remove, , drop = FALSE]
    new_rows_df <- do.call(rbind, new_rows)
    annotation_map <- rbind(annotation_map, new_rows_df)
    rownames(annotation_map) <- NULL  # reset to avoid duplicate row names
    obj$knn.leiden.cluster <- as.factor(obj[[cluster_col, drop = TRUE]])
    log_message(sprintf("Refinement complete: %d clusters split, %d total clusters",
                        length(clusters_to_remove), nrow(annotation_map)))
  } else {
    log_message("Refinement: no clusters needed splitting")
  }

  annotation_map
}

# ==============================================================================
# NK SUBCLUSTER RECOVERY (per-cell, TCR-aware)
# ------------------------------------------------------------------------------
# Recovers NK cells hidden inside clusters that were annotated as T cell but
# are actually a mix (CD8 T majority + NK minority). Cluster-level averages
# mask the NK fraction because CD8 TEMs co-express NKG7/GNLY/KLRD1 — the
# average looks like "cytotoxic T" rather than "T + NK". We rescore each
# cell in the cluster on CD3 vs NK-cytotoxic signatures, and exclude any
# cell with a detected TCR contig (TRAV/TRBV/TRDV/TRGV in CTgene), since
# NK cells do not rearrange TCRs.
#
# Returns a character vector of cell barcodes to reassign as NK. Cluster
# IDs are preserved; only cell-level labels change.
# ==============================================================================
recover_nk_subcluster <- function(obj,
                                  cluster_id,
                                  cluster_col   = "knn.leiden.cluster",
                                  assay         = "RNA",
                                  slot          = "data",
                                  nk_genes      = c("NKG7", "GNLY", "KLRD1"),
                                  t_genes       = c("CD3D", "CD3E"),
                                  cd3_threshold = 0.5,
                                  nk_threshold  = 2.0,
                                  min_nk_frac   = 0.10,
                                  min_nk_cells  = 20,
                                  ctgene_col    = "CTgene",
                                  ctaa_col      = "CTaa",
                                  verbose       = TRUE) {

  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop(sprintf("Cluster column '%s' not found in obj@meta.data", cluster_col))
  }

  cluster_ids <- as.character(obj@meta.data[[cluster_col]])
  cells_in    <- colnames(obj)[cluster_ids == as.character(cluster_id)]
  n_cells     <- length(cells_in)
  if (n_cells == 0) {
    if (verbose) log_message(sprintf("[NK recovery] cluster %s: no cells.", cluster_id))
    return(character(0))
  }

  # Log-normalized expression (Seurat v4 uses `slot`; v5 uses `layer`)
  expr <- tryCatch(
    GetAssayData(obj, assay = assay, slot  = slot),
    error = function(e) GetAssayData(obj, assay = assay, layer = slot)
  )

  nk_present <- intersect(nk_genes, rownames(expr))
  t_present  <- intersect(t_genes,  rownames(expr))
  if (length(nk_present) == 0 || length(t_present) == 0) {
    if (verbose) log_message(sprintf(
      "[NK recovery] cluster %s: missing marker genes (NK=%d, T=%d). Skipping.",
      cluster_id, length(nk_present), length(t_present)))
    return(character(0))
  }

  nk_per_cell  <- Matrix::colSums(expr[nk_present, cells_in, drop = FALSE])
  cd3_per_cell <- Matrix::colSums(expr[t_present,  cells_in, drop = FALSE])

  # Per-cell TCR detection: prefer chain-resolved CTgene, fall back to CTaa
  md <- obj@meta.data[cells_in, , drop = FALSE]
  has_tcr <- rep(FALSE, length(cells_in))
  names(has_tcr) <- cells_in
  if (ctgene_col %in% colnames(md)) {
    ctgene  <- md[[ctgene_col]]
    has_tcr <- !is.na(ctgene) & ctgene != "" &
               grepl("TR[ABDG]V", ctgene, perl = TRUE)
  } else if (ctaa_col %in% colnames(md)) {
    ctaa    <- md[[ctaa_col]]
    has_tcr <- !is.na(ctaa) & ctaa != ""
  } else if (verbose) {
    log_message(sprintf(
      "[NK recovery] cluster %s: no VDJ columns (%s / %s) — skipping TCR exclusion.",
      cluster_id, ctgene_col, ctaa_col))
  }

  is_nk_expr <- (cd3_per_cell < cd3_threshold) & (nk_per_cell > nk_threshold)
  is_nk      <- is_nk_expr & !has_tcr

  n_expr    <- sum(is_nk_expr)
  n_dropped <- sum(is_nk_expr & has_tcr)
  n_nk      <- sum(is_nk)
  frac      <- n_nk / n_cells

  if (verbose) {
    log_message(sprintf(
      "[NK recovery] cluster %s: n=%d, expr-gated=%d, TCR+ dropped=%d, final NK=%d (%.1f%%) [CD3<%.2f & NK>%.2f]",
      cluster_id, n_cells, n_expr, n_dropped, n_nk, frac * 100,
      cd3_threshold, nk_threshold))
  }

  if (n_nk < min_nk_cells || frac < min_nk_frac) {
    if (verbose) log_message(sprintf(
      "[NK recovery] cluster %s: below threshold (need >=%d cells and >=%.0f%%). No recovery.",
      cluster_id, min_nk_cells, min_nk_frac * 100))
    return(character(0))
  }

  names(is_nk) <- cells_in
  cells_in[is_nk]
}

# ------------------------------------------------------------------------------
# Iterate over T-cell clusters whose cluster-average NK signal is high enough
# to warrant per-cell inspection. Returns a named list of cell-barcode vectors
# keyed by cluster ID.
# ------------------------------------------------------------------------------
scan_for_nk_subclusters <- function(obj,
                                    annotation_map,
                                    expr_data,
                                    cluster_col         = "knn.leiden.cluster",
                                    nk_signal_threshold = 5.0,
                                    ...) {
  if (is.null(expr_data) || is.null(expr_data$avg)) return(list())

  t_clusters <- annotation_map$cluster[annotation_map$celltype_broad == "T cell"]
  if (length(t_clusters) == 0) return(list())

  out <- list()
  for (cid in t_clusters) {
    cl_c <- as.character(cid)
    if (!cl_c %in% colnames(expr_data$avg)) next
    nk_sig <- max(vapply(c("NKG7", "GNLY", "KLRD1"), function(g) {
      if (g %in% rownames(expr_data$avg)) expr_data$avg[g, cl_c] else 0
    }, numeric(1)))
    if (nk_sig < nk_signal_threshold) next

    bc <- recover_nk_subcluster(obj, cluster_id = cid,
                                cluster_col = cluster_col, ...)
    if (length(bc) > 0) out[[cl_c]] <- bc
  }
  out
}

renumber_clusters_sequentially <- function(obj, annotation_map,
                                            cluster_col  = "knn.leiden.cluster",
                                            zero_indexed = TRUE) {
  cluster_ids <- as.character(obj[[cluster_col, drop = TRUE]])
  size_tab <- sort(table(cluster_ids), decreasing = TRUE)
  base <- if (zero_indexed) 0L else 1L
  new_levels <- as.character(seq.int(base, base + length(size_tab) - 1L))
  old_to_new <- setNames(new_levels, names(size_tab))

  # unname() strips the lookup-table names so the resulting factor carries
  # cell-barcode names (assigned below), not stale cluster-ID names.
  new_vals <- factor(unname(old_to_new[cluster_ids]), levels = new_levels)
  names(new_vals) <- colnames(obj)
  obj[[cluster_col]] <- new_vals
  Seurat::Idents(obj) <- cluster_col

  annotation_map$cluster <- old_to_new[as.character(annotation_map$cluster)]
  annotation_map <- annotation_map[!is.na(annotation_map$cluster), , drop = FALSE]
  annotation_map <- annotation_map[order(as.integer(annotation_map$cluster)), , drop = FALSE]
  rownames(annotation_map) <- annotation_map$cluster

  log_message(sprintf("Renumbered %d clusters (size-sorted, %s-indexed)",
                      length(size_tab), if (zero_indexed) "zero" else "one"))
  list(obj = obj, annotation_map = annotation_map)
}

# ==============================================================================
# MAIN ANNOTATION FUNCTION
# ==============================================================================
annotate_celltypes <- function(cfg, target = c("all", "eye")) {
  target <- match.arg(target)
  paths  <- get_target_paths(cfg, target)
  log_message("Starting multi-modal cell type annotation (target=", target, ")...")

  # ---- Load integrated object ----
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Integrated object not found: ", obj_path)
    return(invisible(FALSE))
  }
  obj <- readRDS(obj_path)
  meta <- obj[[]]
  cluster_col <- "knn.leiden.cluster"
  clusters <- sort(unique(meta[[cluster_col]]))

  # ---- Evidence source 1: VDJ detection rate ----
  vdj_stats <- get_vdj_metrics(meta, cluster_col)

  # ---- Evidence source 2: FindAllMarkers gene lists ----
  markers_path <- file.path(paths$results_tables, "FindAllMarkers_objclusters.csv")
  marker_df <- if (file.exists(markers_path)) read.csv(markers_path) else NULL

  # ---- Evidence source 3: Cluster-level average expression ----
  log_message("Computing cluster-level expression of canonical markers...")
  expr_data <- tryCatch(
    compute_expression_evidence(obj, cluster_col),
    error = function(e) { log_message("Expression scoring skipped: ", e$message); NULL }
  )

  # Safe cluster-average lookup. Returns 0 when expr_data is missing, or when
  # the gene / cluster isn't present. Captures expr_data from this scope so
  # call sites just pass (gene, cluster_id).
  avg_at <- function(gene, cl_c) {
    if (is.null(expr_data) || !cl_c %in% colnames(expr_data$avg)) return(0)
    if (!gene %in% rownames(expr_data$avg)) return(0)
    expr_data$avg[gene, cl_c]
  }

  # ---- Reference columns (Azimuth, SingleR-HPCA, SingleR-Monaco) ----
  az_col   <- cfg$annotation$azimuth_col %||% "predicted.celltype.l2"
  ref_cols <- c(az_col, "HPCA.label", "Monaco.label")
  ref_cols <- ref_cols[ref_cols %in% colnames(meta)]

  # ---- Build weighted consensus per cluster ----
  log_message("Building weighted consensus for ", length(clusters), " clusters...")

  annotation_map <- lapply(clusters, function(cl) {
    cl_meta <- meta[meta[[cluster_col]] == cl, ]
    cl_row  <- if (nrow(vdj_stats) > 0) vdj_stats[vdj_stats[[cluster_col]] == cl, ] else NULL
    cl_vdj  <- if (!is.null(cl_row) && nrow(cl_row) > 0) cl_row$vdj_prop[1] else 0
    cl_tcr  <- if (!is.null(cl_row) && nrow(cl_row) > 0 && "tcr_prop" %in% names(cl_row)) cl_row$tcr_prop[1] else NA_real_
    cl_bcr  <- if (!is.null(cl_row) && nrow(cl_row) > 0 && "bcr_prop" %in% names(cl_row)) cl_row$bcr_prop[1] else NA_real_

    broad_types <- names(BROAD_MAPPING)
    scores <- setNames(rep(0, length(broad_types)), broad_types)

    # ------------------------------------------------------------------
    # A) Reference votes (weight = 1 each, max 3)
    # ------------------------------------------------------------------
    ref_labels_raw <- list()
    for (rc in ref_cols) {
      mode_label <- get_mode(cl_meta[[rc]])
      if (!is.na(mode_label)) {
        normed    <- normalize_label(mode_label)
        broad_ref <- map_to_broad(normed)
        if (broad_ref %in% broad_types) scores[broad_ref] <- scores[broad_ref] + 1
        ref_labels_raw[[rc]] <- mode_label
      }
    }

    # Azimuth supermajority boost (weight = 3) for rare types that
    # HPCA/Monaco systematically miscall (Platelet, Eryth, pDC, Plasmablast).
    # When >70% of cells in a cluster get the same Azimuth call AND it's one
    # of these rare types, trust Azimuth over the other references.
    if (az_col %in% colnames(cl_meta)) {
      az_labels <- cl_meta[[az_col]]
      az_labels <- az_labels[!is.na(az_labels) & az_labels != ""]
      if (length(az_labels) > 0) {
        az_tab <- sort(table(az_labels), decreasing = TRUE)
        az_top <- names(az_tab)[1]
        az_frac <- az_tab[1] / sum(az_tab)
        az_normed <- normalize_label(az_top)
        az_broad  <- map_to_broad(az_normed)
        rare_types <- c("Platelet", "Eryth", "DC", "Plasma")
        if (az_frac > 0.70 && az_broad %in% rare_types) {
          scores[az_broad] <- scores[az_broad] + 3
          log_message(sprintf("  Cluster %s: Azimuth supermajority %.0f%% %s (%s) → +3 boost",
                              cl, az_frac * 100, az_top, az_broad))
        }
      }
    }

    # ------------------------------------------------------------------
    # B) FindAllMarkers gene overlap (weight = 2)
    # ------------------------------------------------------------------
    cl_markers <- character(0)  # initialise BEFORE conditional
    if (!is.null(marker_df)) {
      cl_markers <- marker_df %>%
        filter(cluster == cl, p_val_adj < 0.05, avg_log2FC > 0.5) %>%
        pull(gene)

      marker_scores <- sapply(names(CANONICAL_MARKERS), function(ct) {
        sum(CANONICAL_MARKERS[[ct]] %in% cl_markers)
      })
      best_marker_ct <- names(which.max(marker_scores))
      if (max(marker_scores) > 0) {
        broad_marker <- if (grepl("T cell|gdT", best_marker_ct)) "T cell"
                        else if (grepl("pDC|cDC", best_marker_ct)) "DC"
                        else best_marker_ct
        if (broad_marker %in% broad_types) {
          scores[broad_marker] <- scores[broad_marker] + 2
        }
      }
    }

    # ------------------------------------------------------------------
    # C) Expression evidence (weight = 3)
    #    Uses CD3D/CD3E (T), CD8A (CD8 T), CD19/MS4A1/CD79A (B),
    #    MZB1 (Plasma), LILRA4 (pDC), FCER1A (cDC), CD14 (Mono)
    # ------------------------------------------------------------------
    expr_scores <- score_cluster_expression(expr_data, cl)
    scores <- scores + expr_scores

    # ------------------------------------------------------------------
    # D) VDJ lineage boost (weight = 5)
    #    Uses CTgene chain types (TRAV/TRBV = TCR → T cell,
    #    IGHV/IGKV/IGLV = BCR → B/Plasma) for definitive lineage ID.
    #    Falls back to expression when chain info is unavailable.
    # ------------------------------------------------------------------
    if (cl_vdj > 0.20) {
      # Primary: chain-type proportions from CTgene (definitive)
      if (!is.na(cl_tcr) && !is.na(cl_bcr) && (cl_tcr + cl_bcr) > 0.05) {
        if (cl_bcr > cl_tcr && cl_bcr > 0.10) {
          # BCR-dominant cluster
          cl_c    <- as.character(cl)
          is_plsm <- avg_at("MZB1", cl_c) > 0.5 && avg_at("CD19", cl_c) < 0.1
          if (is_plsm) scores["Plasma"] <- scores["Plasma"] + 5
          else          scores["B cell"] <- scores["B cell"] + 5
        } else if (cl_tcr > cl_bcr && cl_tcr > 0.10) {
          # TCR-dominant cluster
          scores["T cell"] <- scores["T cell"] + 5
        } else {
          # Mixed chains — boost the majority
          if (cl_tcr >= cl_bcr) scores["T cell"] <- scores["T cell"] + 3
          else                  scores["B cell"] <- scores["B cell"] + 3
        }
      } else {
        # Fallback: no chain info, use markers or expression
        is_b    <- any(c("CD19", "MS4A1", "CD79A") %in% cl_markers)
        is_plsm <- any(c("MZB1", "SDC1") %in% cl_markers)

        if (length(cl_markers) == 0 && !is.null(expr_data)) {
          cl_c <- as.character(cl)
          if (cl_c %in% colnames(expr_data$enriched)) {
            .enr <- function(g) if (g %in% rownames(expr_data$enriched)) expr_data$enriched[g, cl_c] else FALSE
            is_b    <- sum(c(.enr("CD19"), .enr("MS4A1"), .enr("CD79A"))) >= 2
            is_plsm <- .enr("MZB1") && !.enr("CD19")
          }
          if (!is_b) {
            b_det <- sum(c(avg_at("CD19", cl_c) > 0.3,
                           avg_at("MS4A1", cl_c) > 0.3,
                           avg_at("CD79A", cl_c) > 0.3))
            if (b_det >= 2) is_b <- TRUE
            if (avg_at("MZB1", cl_c) > 0.5 && avg_at("CD19", cl_c) < 0.1) is_plsm <- TRUE
          }
        }

        if (is_plsm)        scores["Plasma"] <- scores["Plasma"] + 5
        else if (is_b)      scores["B cell"] <- scores["B cell"] + 5
        else                scores["T cell"] <- scores["T cell"] + 5
      }
    }

    # ------------------------------------------------------------------
    # E) VDJ absence penalty for T cell (weight = -3)
    #    Clusters with <5% VDJ detection are very unlikely to be T cells.
    #    If chain info available, also penalize when TCR specifically is low.
    # ------------------------------------------------------------------
    low_vdj <- cl_vdj < 0.05
    low_tcr <- !is.na(cl_tcr) && cl_tcr < 0.05
    if ((low_vdj || low_tcr) && scores["T cell"] > 0) {
      scores["T cell"] <- max(scores["T cell"] - 3, 0)
    }

    # ------------------------------------------------------------------
    # Expression-based NK rescue
    #   CD16-bright NK cells are routinely mis-called "CD16 Mono" /
    #   "Non-classical monocytes" by Azimuth & SingleR because both
    #   populations express FCGR3A. When expression unambiguously says
    #   NK (cytotoxic+ CD3- myeloid-), override the reference votes.
    #
    #   Uses RAW average expression (not the enrichment matrix) because
    #   the enrichment threshold (above cross-cluster median) fails when
    #   NKG7 is widely expressed across ~70 T cell overclusters.
    # ------------------------------------------------------------------
    {
      cl_c <- as.character(cl)
      # NK cytotoxic signal: any of NKG7/GNLY/KLRD1 meaningfully expressed
      nk_signal   <- max(avg_at("NKG7", cl_c), avg_at("GNLY", cl_c), avg_at("KLRD1", cl_c))
      # T cell signal: CD3D or CD3E
      t_signal    <- max(avg_at("CD3D", cl_c), avg_at("CD3E", cl_c))
      # Myeloid signal: CD14 or LYZ
      mono_signal <- max(avg_at("CD14", cl_c), avg_at("LYZ", cl_c))

      # NK rescue: boost NK above monocyte when cytotoxic > myeloid.
      # Only boost above T cell when CD3 is LOW in absolute terms (<1.0),
      # since CD8 TEM cells co-express high NKG7/GNLY with high CD3.
      if (nk_signal > 0.5 && nk_signal > mono_signal) {
        scores["NK"] <- max(scores["NK"], scores["Monocyte/Macrophage"] + 1L)
        if (t_signal < 1.0) {
          # CD3 near-absent → true NK, not cytotoxic T cell
          scores["NK"] <- max(scores["NK"], scores["T cell"] + 1L)
        }
        log_message(sprintf(
          "  Cluster %s: NK rescue (NKG7/GNLY/KLRD1=%.2f, CD3=%.2f, myeloid=%.2f)",
          cl, nk_signal, t_signal, mono_signal))
      }
    }

    # ------------------------------------------------------------------
    # Pick broad winner — with cross-lineage mixing guardrail
    # ------------------------------------------------------------------
    broad_final <- names(which.max(scores))

    # Ratio-based cross-lineage guardrail: only flag "Mixed" when scores are
    # nearly tied (ratio < 1.5).
    top_two <- sort(scores, decreasing = TRUE)[1:2]
    if (!is.na(top_two[2]) && top_two[2] >= 2) {
      cross_lineage_pairs <- list(
        c("T cell", "Monocyte/Macrophage"),
        c("T cell", "DC"),
        c("T cell", "B cell"),
        c("B cell", "Monocyte/Macrophage"),
        c("NK", "Monocyte/Macrophage"),
        c("Granulocyte", "Monocyte/Macrophage"),
        c("Granulocyte", "T cell")
      )
      top_names <- names(top_two)
      is_cross <- any(vapply(cross_lineage_pairs, function(p) all(p %in% top_names), logical(1)))
      ratio <- top_two[1] / max(top_two[2], 1)
      if (is_cross && ratio < 1.5) {
        log_message(sprintf(
          "WARNING: Cluster %s shows ambiguous cross-lineage signal (%s=%d, %s=%d, ratio=%.2f). Flagging as Mixed.",
          cl, top_names[1], top_two[1], top_names[2], top_two[2], ratio
        ))
        broad_final <- "Mixed"
      }
    }

    # ------------------------------------------------------------------
    # Pick best fine label (normalised, matching broad winner)
    # ------------------------------------------------------------------
    all_fine <- unlist(lapply(ref_cols, function(rc) {
      sapply(cl_meta[[rc]], normalize_label, USE.NAMES = FALSE)
    }))
    valid_fine <- all_fine[!is.na(all_fine) & sapply(all_fine, map_to_broad) == broad_final]
    fine_final <- if (length(valid_fine) > 0) {
      names(sort(table(valid_fine), decreasing = TRUE))[1]
    } else {
      broad_final
    }

    # ------------------------------------------------------------------
    # For T cells: refine CD4 vs CD8 using expression when fine label
    # doesn't already specify a subtype
    # ------------------------------------------------------------------
    if (broad_final == "T cell" && !is.null(expr_data)) {
      t_sub <- determine_t_subtype(expr_data, cl)
      if (!grepl("CD4|CD8|gdT|dnT|MAIT|Treg|Tfh", fine_final)) {
        fine_final <- paste0(t_sub, " T")
      }
      # Proliferation check: MKI67/TOP2A enriched → label as "Proliferating"
      # rather than a specific CD4/CD8 subtype, since cycling cells are mixed
      cl_c <- as.character(cl)
      if (avg_at("MKI67", cl_c) > 1 || avg_at("TOP2A", cl_c) > 1) {
        fine_final <- "T Proliferating"
      }
    }

    # ------------------------------------------------------------------
    # CD8-like cluster with low TCR recovery → reconsider gdT
    #   gdT cells: TRDC/TRDV expression, often poor TRAV/TRBV capture
    #   True CD8 T: CD3D+, CD8A+, robust TRAV/TRBV recovery
    #
    # Note: CD8→NK reclassification is now handled per-cell by the post-hoc
    # recover_nk_subcluster pass (see end of annotate_celltypes), which can
    # rescue NK cells from MIXED clusters that this cluster-level check
    # misses (e.g. CD8-majority clusters with a hidden NK minority).
    # ------------------------------------------------------------------
    if (broad_final == "T cell" && grepl("CD8|gdT", fine_final)) {
      tcr_low <- (!is.na(cl_tcr) && cl_tcr < 0.15) || (is.na(cl_tcr) && cl_vdj < 0.15)
      if (tcr_low) {
        cl_c    <- as.character(cl)
        nk_sig  <- max(avg_at("NKG7", cl_c), avg_at("GNLY", cl_c), avg_at("KLRD1", cl_c))
        cd3_sig <- max(avg_at("CD3D", cl_c), avg_at("CD3E", cl_c))
        gdt_sig <- max(avg_at("TRDC", cl_c), avg_at("TRDV1", cl_c), avg_at("TRDV2", cl_c))

        if (gdt_sig > 0.3 && (cd3_sig > 0.2 || nk_sig > 0.3)) {
          # Delta chain expression present — likely gdT (poor standard TCR capture)
          fine_final <- "gdT"
          log_message(sprintf(
            "  Cluster %s: low TCR (%.1f%%) but TRDC=%.2f → confirmed gdT",
            cl, (if (!is.na(cl_tcr)) cl_tcr else cl_vdj) * 100, gdt_sig))
        } else if (nk_sig > 0.3 && cd3_sig > 0.3) {
          # Both signals present — could be NKT or mixed, flag for review
          log_message(sprintf(
            "  Cluster %s: low TCR (%.1f%%), NKG7=%.2f, CD3=%.2f → ambiguous CD8/NK, keeping %s (per-cell recovery will handle NK fraction)",
            cl, (if (!is.na(cl_tcr)) cl_tcr else cl_vdj) * 100, nk_sig, cd3_sig, fine_final))
        }
      }
    }

    # ------------------------------------------------------------------
    # Resolve still-generic broad labels to specific fine names using
    # expression evidence (CD14 vs FCGR3A for Mono, pDC/cDC1/cDC2 for DC,
    # Neutrophil/Mast/Basophil for Granulocyte).
    # ------------------------------------------------------------------
    # If fine_final is still just the broad label AND expression can
    # resolve it further, do so now. Two subtypes only:
    #   CD16 Mono : FCGR3A > CD14 and FCGR3A > 0.5
    #   CD14 Mono : everything else in the Monocyte/Macrophage class
    if (fine_final == "Monocyte/Macrophage") {
      cl_c        <- as.character(cl)
      cd14_expr   <- avg_at("CD14",   cl_c)
      fcgr3a_expr <- avg_at("FCGR3A", cl_c)
      fine_final  <- if (fcgr3a_expr > cd14_expr && fcgr3a_expr > 0.5) "CD16 Mono" else "CD14 Mono"
    }
    if (fine_final == "DC" && !is.null(expr_data)) {
      cl_c <- as.character(cl)
      if (cl_c %in% colnames(expr_data$enriched)) {
        enr_check <- function(g) if (g %in% rownames(expr_data$enriched)) expr_data$enriched[g, cl_c] else FALSE
        if (enr_check("LILRA4") || enr_check("IRF7")) fine_final <- "pDC"
        else if (enr_check("CLEC10A")) fine_final <- "cDC2"
        else if (enr_check("CLEC9A") || enr_check("XCR1")) fine_final <- "cDC1"
      }
    }
    if (fine_final == "Granulocyte" && !is.null(expr_data)) {
      cl_c <- as.character(cl)
      if (cl_c %in% colnames(expr_data$enriched)) {
        enr_check <- function(g) if (g %in% rownames(expr_data$enriched)) expr_data$enriched[g, cl_c] else FALSE
        if (enr_check("TPSAB1") || enr_check("TPSB2")) fine_final <- "Mast cell"
        else if (enr_check("CLC") || enr_check("HDC")) fine_final <- "Basophil"
        else if (enr_check("CSF3R") || enr_check("FCGR3B")) fine_final <- "Neutrophil"
      }
    }

    # ------------------------------------------------------------------
    # Concordance: how many references agree with the final broad label
    # ------------------------------------------------------------------
    ref_broads  <- sapply(ref_labels_raw, function(lbl) map_to_broad(normalize_label(lbl)))
    n_agree     <- sum(ref_broads == broad_final, na.rm = TRUE)

    data.frame(
      cluster        = as.character(cl),
      celltype       = fine_final,
      celltype_broad = broad_final,
      vdj_confidence = cl_vdj,
      top_score      = max(scores),
      n_ref_agree    = n_agree,
      n_refs_used    = length(ref_broads),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  # ---- Data-driven SNN refinement ----
  # Eye target shallow-merges cfg$eye_refinement on top of cfg$refinement,
  # so e.g. setting eye_refinement.enable: false disables refinement only
  # for the eye branch.
  ref_cfg <- cfg$refinement %||% list(enable = TRUE)
  if (target == "eye") {
    ref_cfg <- utils::modifyList(ref_cfg, cfg$eye_refinement %||% list())
  }
  if (isTRUE(ref_cfg$enable) && !is.null(obj@misc$knn_leiden$snn_graph)) {
    log_message("Running data-driven SNN cluster refinement...")
    annotation_map <- refine_clusters_by_evidence(
      obj            = obj,
      annotation_map = annotation_map,
      snn_graph      = obj@misc$knn_leiden$snn_graph,
      expr_data      = expr_data,
      cluster_col    = cluster_col,
      ref_cols       = ref_cols,
      cfg            = ref_cfg
    )
  }

  # ---- Post-hoc Mono rescue + orphan cluster recovery ----
  # Two paths, same logic (check myeloid markers, assign CD14/CD16 Mono):
  #   (1) "Unresolved" rows already in annotation_map (reference voting failed)
  #   (2) Orphan Leiden cluster IDs present in obj$knn.leiden.cluster but
  #       absent from annotation_map (SNN refinement / dedup artifact). These
  #       would otherwise produce "Unresolved" per-cell labels for ALL their
  #       cells at the downstream NA-fill step, which wipes out whole
  #       monocyte populations. Caught here and given proper mapping rows.
  #
  # Classification: CD16 Mono when FCGR3A-dominant, CD14 Mono otherwise.
  if (!is.null(expr_data)) {
    # Path 2: build stub rows for orphan Leiden cluster IDs so they get
    # processed by the same Mono-rescue loop below.
    cell_cluster_ids_pre <- as.character(obj[[cluster_col]][, 1])
    present_in_map       <- as.character(annotation_map$cluster)
    orphan_ids           <- setdiff(unique(cell_cluster_ids_pre), present_in_map)
    if (length(orphan_ids) > 0) {
      log_message(sprintf("Found %d orphan Leiden clusters missing from annotation map: %s",
                          length(orphan_ids), paste(orphan_ids, collapse = ", ")))
      orphan_stub <- data.frame(
        cluster        = orphan_ids,
        celltype       = "Unresolved",
        celltype_broad = "Other",
        vdj_confidence = NA_real_,
        top_score      = NA_real_,
        n_ref_agree    = NA_integer_,
        n_refs_used    = NA_integer_,
        stringsAsFactors = FALSE
      )
      annotation_map <- rbind(annotation_map, orphan_stub)
    }

    unres_rows <- which(is.na(annotation_map$celltype) |
                        annotation_map$celltype == "" |
                        annotation_map$celltype == "Unresolved")
    n_mono_rescued <- 0
    for (i in unres_rows) {
      cl_c <- as.character(annotation_map$cluster[i])
      if (!cl_c %in% colnames(expr_data$avg)) next
      cd14   <- avg_at("CD14", cl_c)
      lyz    <- avg_at("LYZ", cl_c)
      fcgr3a <- avg_at("FCGR3A", cl_c)
      mono_signal <- max(cd14, lyz)
      # Trigger: strong classical-Mono signature OR strong non-classical FCGR3A
      if (mono_signal > 0.5 || fcgr3a > 0.5) {
        new_ct <- if (fcgr3a > cd14 && fcgr3a > 0.5) "CD16 Mono" else "CD14 Mono"
        log_message(sprintf(
          "  Cluster %s: Mono rescue -> %s (CD14=%.2f, LYZ=%.2f, FCGR3A=%.2f)",
          cl_c, new_ct, cd14, lyz, fcgr3a))
        annotation_map$celltype[i]       <- new_ct
        annotation_map$celltype_broad[i] <- "Monocyte/Macrophage"
        n_mono_rescued <- n_mono_rescued + 1
      }
    }
    if (n_mono_rescued > 0) {
      log_message(sprintf("Mono rescue: reclassified %d clusters", n_mono_rescued))
    }
  }

  # ---- Low-quality cluster detection (eye target) -------------------------------
  # Clusters whose FindAllMarkers output is dominated by NEGATIVE log2FC (every
  # gene is below the cross-cluster median) are typically debris / low-RNA /
  # ambient-background droplets. Reference annotators still produce a label
  # (e.g. "CD8 TEM") because they classify by relative similarity, but the
  # cluster has no real biology. Flag these as "Low Quality" so they're
  # excluded from downstream interpretation.
  lq_cfg <- cfg$eye_qc$low_quality_filter %||% list()
  if (target == "eye" && isTRUE(lq_cfg$enable) && !is.null(marker_df)) {
    min_pos      <- as.integer(lq_cfg$min_positive_markers %||% 5L)
    max_neg_frac <- as.numeric(lq_cfg$max_negative_fraction %||% 0.7)
    qc_summary <- marker_df %>%
      dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        n_pos    = sum(avg_log2FC > 0.5, na.rm = TRUE),
        n_total  = dplyr::n(),
        neg_frac = sum(avg_log2FC < 0, na.rm = TRUE) / max(dplyr::n(), 1L),
        .groups  = "drop"
      )
    bad <- as.character(qc_summary$cluster[qc_summary$n_pos < min_pos |
                                           qc_summary$neg_frac > max_neg_frac])
    if (length(bad) > 0) {
      log_message(sprintf("Low-quality clusters detected (eye): %s",
                          paste(bad, collapse = ", ")))
      idx <- which(as.character(annotation_map$cluster) %in% bad)
      annotation_map$celltype[idx]       <- "Low Quality"
      annotation_map$celltype_broad[idx] <- "Other"
    }
  }

  # ---- Validate: every cluster must have a defined celltype and celltype_broad ----
  na_ct <- is.na(annotation_map$celltype) | annotation_map$celltype == ""
  na_broad <- is.na(annotation_map$celltype_broad) | annotation_map$celltype_broad == ""

  if (any(na_ct)) {
    bad_cls <- annotation_map$cluster[na_ct]
    log_message(sprintf("WARNING: %d clusters have no celltype label: %s. Assigning 'Unresolved'.",
                        sum(na_ct), paste(bad_cls, collapse = ", ")))
    annotation_map$celltype[na_ct] <- "Unresolved"
  }
  if (any(na_broad)) {
    bad_cls <- annotation_map$cluster[na_broad]
    log_message(sprintf("WARNING: %d clusters have no celltype_broad label: %s. Assigning 'Other'.",
                        sum(na_broad), paste(bad_cls, collapse = ", ")))
    annotation_map$celltype_broad[na_broad] <- "Other"
  }

  # Log the final cluster-to-celltype mapping
  log_message(sprintf("Final annotation: %d clusters annotated", nrow(annotation_map)))
  for (i in seq_len(nrow(annotation_map))) {
    log_message(sprintf("  Cluster %-4s -> %-20s (%s)",
                        annotation_map$cluster[i],
                        annotation_map$celltype[i],
                        annotation_map$celltype_broad[i]))
  }

  # ---- Apply to Seurat object ----
  annotation_map$cluster <- as.character(annotation_map$cluster)
  # Deduplicate: if any cluster ID appears more than once, keep the first row
  if (anyDuplicated(annotation_map$cluster)) {
    dup_ids <- annotation_map$cluster[duplicated(annotation_map$cluster)]
    log_message(sprintf("WARNING: Removing %d duplicate cluster rows: %s",
                        length(dup_ids), paste(unique(dup_ids), collapse = ", ")))
    annotation_map <- annotation_map[!duplicated(annotation_map$cluster), , drop = FALSE]
  }
  rownames(annotation_map) <- annotation_map$cluster

  cell_cluster_ids <- as.character(obj[[cluster_col]][, 1])
  obj$celltype       <- annotation_map[cell_cluster_ids, "celltype"]
  obj$celltype_broad <- annotation_map[cell_cluster_ids, "celltype_broad"]

  # Validate no cells fell through
  n_na_ct <- sum(is.na(obj$celltype))
  if (n_na_ct > 0) {
    unmapped_cls <- setdiff(unique(cell_cluster_ids), rownames(annotation_map))
    log_message(sprintf("WARNING: %d cells have NA celltype (unmapped clusters: %s). Assigning 'Unresolved'.",
                        n_na_ct, paste(unmapped_cls, collapse = ", ")))
    obj$celltype[is.na(obj$celltype)] <- "Unresolved"
    obj$celltype_broad[is.na(obj$celltype_broad)] <- "Other"
  }

  # ---- Per-cell NK subcluster recovery ----
  # Rescue NK cells hidden inside T-cell-majority clusters. Uses per-cell
  # CD3/NK-cytotoxic scoring and excludes cells with detected TCR contigs.
  # Reads thresholds from cfg$nk_recovery when present; otherwise defaults.
  nk_cfg <- cfg$nk_recovery %||% list()
  if (!isFALSE(nk_cfg$enable)) {
    log_message("Running per-cell NK subcluster recovery on T-cell clusters...")
    nk_recovery <- scan_for_nk_subclusters(
      obj                 = obj,
      annotation_map      = annotation_map,
      expr_data           = expr_data,
      cluster_col         = cluster_col,
      nk_signal_threshold = nk_cfg$nk_signal_threshold %||% 5.0,
      cd3_threshold       = nk_cfg$cd3_threshold       %||% 0.5,
      nk_threshold        = nk_cfg$nk_threshold        %||% 2.0,
      min_nk_frac         = nk_cfg$min_nk_frac         %||% 0.10,
      min_nk_cells        = nk_cfg$min_nk_cells        %||% 20,
      verbose             = TRUE
    )

    if (length(nk_recovery) > 0) {
      total_recovered <- sum(vapply(nk_recovery, length, integer(1)))
      log_message(sprintf("NK recovery: reassigned %d cells across %d cluster(s).",
                          total_recovered, length(nk_recovery)))

      majority_frac   <- nk_cfg$majority_frac %||% 0.50
      eye_ann         <- cfg$eye_annotation %||% list()
      force_split     <- target == "eye" && isTRUE(eye_ann$force_nk_split)
      split_min_cells <- as.integer(eye_ann$nk_split_min_cells %||% 50L)
      # Convert factor → character so we can introduce new cluster IDs from
      # force_split without first declaring factor levels.
      if (force_split) {
        obj[[cluster_col]] <- as.character(obj[[cluster_col, drop = TRUE]])
      }
      cell_cluster_id <- as.character(obj[[cluster_col]][, 1])
      for (cid in names(nk_recovery)) {
        bc        <- nk_recovery[[cid]]
        in_cl     <- cell_cluster_id == cid
        n_cluster <- sum(in_cl)
        frac_nk   <- length(bc) / max(n_cluster, 1)

        if (frac_nk >= majority_frac) {
          # Majority promotion: relabel the ENTIRE cluster (mapping row + every
          # cell in it). Only updating `bc` would leave the non-gate-passing
          # cells with stale labels, which 37's per-cell override would then
          # pull back out of the promoted NK display cluster.
          obj$celltype[in_cl]       <- "NK"
          obj$celltype_broad[in_cl] <- "NK"
          idx <- which(annotation_map$cluster == cid)
          if (length(idx) == 1) {
            annotation_map$celltype[idx]       <- "NK"
            annotation_map$celltype_broad[idx] <- "NK"
          }
          log_message(sprintf(
            "  Cluster %s: %.0f%% cells reassigned to NK (>=%.0f%%) → cluster and all %d cells relabeled NK",
            cid, frac_nk * 100, majority_frac * 100, n_cluster))
        } else if (force_split && length(bc) >= split_min_cells) {
          # Eye-target force-split: minority NK fraction is large enough to
          # warrant its own cluster ID. Allocate a new ID, move the recovered
          # cells into it, append a row to annotation_map. This makes the NK
          # subpopulation visible as a distinct cluster on UMAP rather than
          # being hidden inside the parent cytotoxic-T cluster.
          existing_ids <- suppressWarnings(as.integer(
            unique(c(as.character(annotation_map$cluster), cell_cluster_id))
          ))
          existing_ids <- existing_ids[!is.na(existing_ids)]
          new_id <- as.character(max(existing_ids, 0L) + 1L)
          obj[[cluster_col]][bc, 1] <- new_id
          obj$celltype[bc]          <- "NK"
          obj$celltype_broad[bc]    <- "NK"
          annotation_map <- rbind(annotation_map, data.frame(
            cluster        = new_id,
            celltype       = "NK",
            celltype_broad = "NK",
            vdj_confidence = NA_real_,
            top_score      = NA_real_,
            n_ref_agree    = NA_integer_,
            n_refs_used    = NA_integer_,
            stringsAsFactors = FALSE
          ))
          rownames(annotation_map) <- annotation_map$cluster
          cell_cluster_id <- as.character(obj[[cluster_col]][, 1])
          log_message(sprintf(
            "  Cluster %s: force-split NK fraction (%d cells, %.1f%%) into new cluster %s",
            cid, length(bc), frac_nk * 100, new_id))
        } else {
          # Minority recovery: only update the gate-passing cells, leave the
          # cluster label alone. 37 will route these cells to the existing NK
          # display cluster via the per-cell override.
          obj$celltype[bc]       <- "NK"
          obj$celltype_broad[bc] <- "NK"
        }
      }
      # Refresh factor levels if force_split added new IDs
      if (force_split) {
        obj[[cluster_col]] <- as.factor(obj[[cluster_col, drop = TRUE]])
        Seurat::Idents(obj) <- cluster_col
      }
    } else {
      log_message("NK recovery: no clusters met recovery criteria.")
    }
  }

  # ---- Eye mature-label preference (target == "eye") --------------------------
  # Tissue context: in eye, classical/non-classical monocytes are infiltrating
  # tissue macrophages; antibody-secreting B cells are plasma cells, not
  # memory B. Apply this remap AFTER all auto-annotation but BEFORE manual
  # overrides so user-set per-cluster overrides still win.
  eye_ann <- cfg$eye_annotation %||% list()
  if (target == "eye" && isTRUE(eye_ann$prefer_mature_labels)) {
    mzb1_thr <- as.numeric(eye_ann$mzb1_threshold %||% 0.5)
    cell_cluster_id <- as.character(obj[[cluster_col]][, 1])
    for (i in seq_len(nrow(annotation_map))) {
      cl_c   <- as.character(annotation_map$cluster[i])
      ct     <- annotation_map$celltype[i]
      in_cl  <- cell_cluster_id == cl_c
      new_ct <- NULL; new_br <- NULL

      # Mono -> Macrophage
      if (!is.na(ct) && grepl("Mono$|Monocyte", ct)) {
        new_ct <- "Macrophage"
        new_br <- "Macrophage"
      }
      # B intermediate / B memory with MZB1 signal -> Plasma
      if (!is.na(ct) && grepl("^B (memory|intermediate|cell)$|^B$", ct)) {
        if (avg_at("MZB1", cl_c) > mzb1_thr) {
          new_ct <- "Plasma"
          new_br <- "Plasma"
        }
      }

      if (!is.null(new_ct)) {
        log_message(sprintf(
          "  Eye mature remap: cluster %s: %s -> %s",
          cl_c, ct, new_ct))
        annotation_map$celltype[i]       <- new_ct
        annotation_map$celltype_broad[i] <- new_br
        if (any(in_cl)) {
          obj$celltype[in_cl]       <- new_ct
          obj$celltype_broad[in_cl] <- new_br
        }
      }
    }
  }

  # ---- Manual per-cluster overrides (config-driven) ----------------------
  # Applied last so they win over all auto-annotation (consensus voting,
  # rule-based boosts, NK rescue). Updates BOTH the annotation_map row AND
  # every cell's obj$celltype / obj$celltype_broad — otherwise 37's per-cell
  # override pathway would treat the leftover cells as mismatched and spawn
  # ghost display clusters with the old label.
  manual_ov <- if (target == "eye") {
    cfg$eye_manual_overrides %||% list()
  } else {
    cfg$manual_overrides %||% list()
  }
  if (length(manual_ov) > 0) {
    cell_cluster_id <- as.character(obj[[cluster_col]][, 1])
    for (cid in names(manual_ov)) {
      ov       <- manual_ov[[cid]]
      new_ct   <- ov$celltype
      new_br   <- ov$celltype_broad %||% ov$celltype
      idx      <- which(as.character(annotation_map$cluster) == cid)
      in_cl    <- cell_cluster_id == cid
      if (length(idx) == 1) {
        old_ct <- annotation_map$celltype[idx]
        annotation_map$celltype[idx]       <- new_ct
        annotation_map$celltype_broad[idx] <- new_br
        log_message(sprintf(
          "  Manual override: cluster %s: %s -> %s (broad: %s) [%d cells]",
          cid, old_ct, new_ct, new_br, sum(in_cl)))
      } else {
        log_message(sprintf(
          "  Manual override: cluster %s not in annotation_map — skipping mapping row update",
          cid))
      }
      if (any(in_cl)) {
        obj$celltype[in_cl]       <- new_ct
        obj$celltype_broad[in_cl] <- new_br
      }
    }
  }

  # ---- Cluster ID hygiene ----------------------------------------------------
  # Renumber sequentially (0..N-1) by descending size to remove gaps left
  # by SNN refinement / force_split allocation.
  hygiene_cfg <- cfg$cluster_hygiene %||% list()
  if (!isFALSE(hygiene_cfg$renumber)) {
    res <- renumber_clusters_sequentially(obj, annotation_map, cluster_col)
    obj            <- res$obj
    annotation_map <- res$annotation_map
  }

  # ---- Save outputs ----
  ensure_dir(paths$results_tables)
  write.csv(annotation_map,
            file.path(paths$results_tables, "cluster_celltype_mapping.csv"),
            row.names = FALSE)
  log_message("Cluster -> cell type mapping saved.")

  # Composition table (cells per sample × broad type, with metadata)
  comp_df <- obj[[]] %>%
    dplyr::count(orig.ident, celltype_broad, name = "n_cells") %>%
    group_by(orig.ident) %>%
    mutate(proportion = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    left_join(distinct(obj[[]], orig.ident, Tissue_1, Phenotype_2, Etiology, Subject),
              by = "orig.ident")
  write.csv(comp_df,
            file.path(paths$results_tables, "celltype_composition.csv"),
            row.names = FALSE)
  log_message("Cell type composition table saved.")

  saveRDS(obj, obj_path)
  log_message("Annotated object saved with weighted consensus labels.")

  invisible(TRUE)
}
