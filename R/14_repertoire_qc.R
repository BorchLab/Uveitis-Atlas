# R/14_repertoire_qc.R
#
# Repertoire-QC diagnostics that gate the compartment subset choices.
#
# Two entry points:
#   * run_ctaa_filter_diagnostic(cfg)
#       T0.1 from the 2026-05-13 plan. Cross-tabulates CTaa+ cells in the
#       myeloid parent clusters against doublet calls, lineage module
#       scores, lens scores, and lymphocyte top-clone membership, then
#       reports diversity under four candidate filter options (a/b/c/d).
#       Read-only with respect to the eye object — outputs are tables +
#       a comparison PDF used to lock cfg$compartments$myeloid_ctaa_filter.
#   * run_clone_definition_sensitivity(cfg)
#       T2.3. Compares eye-blood overlap counts under CTaa / CTnt / CTgene
#       / CTstrict for each modality so the manuscript can cite a single
#       locked definition with a sensitivity table backing it.
#
# Both functions read the eye sub-atlas at
# outputs/objects/eye/IntegratedSeuratObject.rds and never modify it.
#
# Assumptions about metadata (set up in R/02_ingest_data.R +
# R/03_annotate_AIRR.R): Subject, Tissue_1, Cohort, doublet.class,
# doublet.score, CTaa, CTstrict, CTgene, clonalFrequency,
# knn.leiden.cluster, celltype, celltype_broad. Missing columns are
# tolerated with NA fallbacks.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
})

# --- Diversity helpers (manual; avoids a scRepertoire round-trip) -----------

# Shannon diversity on a vector of clone IDs (NA / empty strings excluded).
# Returns NA when fewer than 2 clones remain.
.shannon <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) < 2L) return(NA_real_)
  p <- prop.table(table(x))
  -sum(p * log(p))
}

# Inverse Simpson on a vector of clone IDs.
.inv_simpson <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) < 2L) return(NA_real_)
  p <- prop.table(table(x))
  1 / sum(p^2)
}

# Resolve the eye-level cluster column with the same fallback the compartment
# subset uses.
.eye_cluster_col <- function(meta) {
  if ("knn.leiden.cluster" %in% colnames(meta)) "knn.leiden.cluster" else "seurat_clusters"
}

# Identify the parent cells for one compartment. Supports both the current
# cluster-ID schema (`parent_clusters`) and the post-T1.2 celltype-key
# schema (`celltype_keys`) so this diagnostic survives the Tier 1 refactor.
.compartment_parent_cells <- function(eye_obj, cfg, cmp) {
  meta <- eye_obj[[]]
  cmp_cfg <- cfg$compartments[[cmp]]
  ctb_col <- resolve_celltype_broad(meta)
  if (!is.null(cmp_cfg$celltype_keys) && !is.null(ctb_col)) {
    keys <- as.character(cmp_cfg$celltype_keys)
    colnames(eye_obj)[as.character(meta[[ctb_col]]) %in% keys]
  } else {
    cluster_col <- .eye_cluster_col(meta)
    parents <- as.character(cmp_cfg$parent_clusters)
    colnames(eye_obj)[as.character(meta[[cluster_col]]) %in% parents]
  }
}

# Per-cell lineage / lens module scores via Seurat::AddModuleScore on a
# subset of the eye object. Fallback to Matrix::colSums of log-norm
# expression when AddModuleScore fails on a tiny subset.
.score_lineage_panels <- function(eye_obj, cells, cfg) {
  panels <- list(
    T     = c("CD3D", "CD3E", "CD8A", "CD4"),
    B     = c("CD19", "MS4A1", "CD79A"),
    lens  = cfg$eye_qc$lens_filter$genes
  )
  sub <- subset(eye_obj, cells = cells)
  DefaultAssay(sub) <- "RNA"
  sub <- NormalizeData(sub, verbose = FALSE)

  scores <- data.frame(cell_id = colnames(sub), stringsAsFactors = FALSE)
  for (nm in names(panels)) {
    genes <- intersect(panels[[nm]], rownames(sub))
    if (length(genes) == 0L) {
      scores[[paste0("score_", nm)]] <- NA_real_
      next
    }
    score_name <- paste0("Score_", nm)
    sub_scored <- tryCatch(
      AddModuleScore(sub, features = list(genes), name = score_name,
                     seed = cfg$seed %||% 42L,
                     ctrl = min(50L, max(5L, length(genes)))),
      error = function(e) NULL
    )
    if (is.null(sub_scored)) {
      mat <- GetAssayData(sub, assay = "RNA", layer = "data")[genes, , drop = FALSE]
      scores[[paste0("score_", nm)]] <- Matrix::colSums(mat)
    } else {
      scores[[paste0("score_", nm)]] <- sub_scored[[paste0(score_name, "1"), drop = TRUE]]
    }
  }
  scores
}

# Tag each row with the four candidate filter decisions. Returns the
# input data.frame plus columns kept_filter_a/b/c/d (TRUE = keep).
#
#   a — drop any cell with a productive CTaa (current default)
#   b — drop CTaa+ matching a top-N expanded lymphocyte clone (carryover/doublet)
#   c — keep all; flag CTaa+ in metadata
#   d — drop CTaa+ AND lineage score above median (T or B)
.apply_filter_options <- function(per_cell, top_n = 10L) {
  has_clone <- per_cell$has_clone
  ct <- per_cell$CTstrict

  # Top-N expanded lymphocyte clones per Subject (lymphocyte cells only).
  lymph_idx <- per_cell$is_lymphocyte_clone & !is.na(ct) & nzchar(as.character(ct))
  per_cell$is_top_n_lymph_clone <- FALSE
  if (any(lymph_idx)) {
    top_clones <- per_cell %>%
      dplyr::filter(lymph_idx) %>%
      dplyr::group_by(Subject) %>%
      dplyr::count(CTstrict, sort = TRUE) %>%
      dplyr::slice_max(n, n = top_n, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(Subject, CTstrict) %>%
      as.data.frame()
    top_key  <- paste(top_clones$Subject, top_clones$CTstrict, sep = "||")
    cell_key <- paste(per_cell$Subject,   ct,                 sep = "||")
    per_cell$is_top_n_lymph_clone <- cell_key %in% top_key & has_clone
  }

  median_T <- stats::median(per_cell$score_T, na.rm = TRUE)
  median_B <- stats::median(per_cell$score_B, na.rm = TRUE)
  high_lineage <- (per_cell$score_T > median_T) | (per_cell$score_B > median_B)
  high_lineage[is.na(high_lineage)] <- FALSE

  per_cell$kept_filter_a <- !has_clone
  per_cell$kept_filter_b <- !(has_clone & per_cell$is_top_n_lymph_clone)
  per_cell$kept_filter_c <- rep(TRUE, nrow(per_cell))
  per_cell$kept_filter_d <- !(has_clone & high_lineage)
  per_cell
}

# Build the per-cell long table for the myeloid parent set + lymphocyte
# reference set used for top-N clone identification.
.build_myeloid_ctaa_table <- function(eye_obj, cfg) {
  log_message("CTaa diagnostic: collecting myeloid parent cells")
  myeloid_cells <- .compartment_parent_cells(eye_obj, cfg, "myeloid")
  bcell_cells   <- .compartment_parent_cells(eye_obj, cfg, "bcell")
  tcell_cells   <- .compartment_parent_cells(eye_obj, cfg, "tcell")
  all_cells     <- unique(c(myeloid_cells, bcell_cells, tcell_cells))

  log_message(sprintf("  myeloid=%d  bcell=%d  tcell=%d  total=%d",
                      length(myeloid_cells), length(bcell_cells),
                      length(tcell_cells),   length(all_cells)))

  meta <- eye_obj[[]]
  meta$cell_id <- rownames(meta)
  cluster_col  <- .eye_cluster_col(meta)

  pick <- function(col, default = NA_character_, fn = as.character) {
    if (col %in% colnames(meta)) fn(meta[[col]]) else rep(default, nrow(meta))
  }

  ctb_col <- resolve_celltype_broad(meta)
  ct_col  <- resolve_celltype(meta)
  base <- data.frame(
    cell_id         = meta$cell_id,
    orig.ident      = meta$orig.ident,
    Subject         = pick("Subject"),
    Tissue_1        = pick("Tissue_1"),
    Cohort          = pick("Cohort"),
    eye_cluster     = as.character(meta[[cluster_col]]),
    celltype_eye    = if (!is.null(ct_col)) as.character(meta[[ct_col]]) else NA_character_,
    celltype_broad  = if (!is.null(ctb_col)) as.character(meta[[ctb_col]]) else NA_character_,
    CTaa            = pick("CTaa"),
    CTstrict        = pick("CTstrict"),
    CTgene          = pick("CTgene"),
    clonalFrequency = if ("clonalFrequency" %in% colnames(meta)) meta$clonalFrequency else NA_integer_,
    doublet.class   = pick("doublet.class"),
    doublet.score   = if ("doublet.score" %in% colnames(meta)) meta$doublet.score else NA_real_,
    stringsAsFactors = FALSE
  )
  base$has_clone <- !is.na(base$CTaa) & nzchar(base$CTaa)
  base$clone_class <- dplyr::case_when(
    grepl("^IG", base$CTgene) ~ "BCR",
    grepl("^TR", base$CTgene) ~ "TCR",
    TRUE                      ~ NA_character_
  )
  base$compartment_parent <- dplyr::case_when(
    base$cell_id %in% myeloid_cells ~ "myeloid",
    base$cell_id %in% bcell_cells   ~ "bcell",
    base$cell_id %in% tcell_cells   ~ "tcell",
    TRUE                            ~ "other"
  )
  base$is_lymphocyte_clone <- base$compartment_parent %in% c("bcell", "tcell") &
                              !is.na(base$clone_class)

  # Score only the cells we'll actually use to keep this fast.
  scores <- .score_lineage_panels(eye_obj, all_cells, cfg)
  base   <- dplyr::left_join(base, scores, by = "cell_id")

  with_filters <- .apply_filter_options(
    dplyr::filter(base, compartment_parent %in% c("myeloid", "bcell", "tcell")),
    top_n = cfg$repertoire$top_n_lymph_for_filter %||% 10L
  )
  list(
    per_cell  = dplyr::filter(with_filters, compartment_parent == "myeloid"),
    lymph_tbl = dplyr::filter(with_filters, compartment_parent %in% c("bcell", "tcell"))
  )
}

# Summarize cells kept/dropped + diversity under each filter.
.summarize_filter_comparison <- function(per_cell, lymph_tbl) {
  filter_specs <- list(
    a = "drop any CTaa+ (current default)",
    b = "drop CTaa+ matching top-N lymphocyte clones (carryover/doublet)",
    c = "keep all, flag CTaa+ in metadata",
    d = "drop CTaa+ AND lineage_score above median"
  )
  shannon_lymph   <- .shannon(lymph_tbl$CTstrict[lymph_tbl$has_clone])
  isimpson_lymph  <- .inv_simpson(lymph_tbl$CTstrict[lymph_tbl$has_clone])

  rows <- lapply(names(filter_specs), function(letter) {
    keep_col <- paste0("kept_filter_", letter)
    kept     <- per_cell[[keep_col]]
    dropped  <- !kept
    n_drop   <- sum(dropped)
    n_drop_d <- if (n_drop > 0L)
                  sum(per_cell$doublet.class[dropped] == "doublet", na.rm = TRUE)
                else 0L
    data.frame(
      filter_option       = letter,
      description         = filter_specs[[letter]],
      n_myeloid_total     = nrow(per_cell),
      n_ctaa_positive     = sum(per_cell$has_clone),
      n_kept              = sum(kept),
      n_dropped           = n_drop,
      n_dropped_doublet   = n_drop_d,
      pct_dropped_doublet = if (n_drop > 0L) round(100 * n_drop_d / n_drop, 1) else NA_real_,
      n_top_n_lymph_in_myeloid = sum(per_cell$is_top_n_lymph_clone, na.rm = TRUE),
      shannon_myeloid_clones    = .shannon(per_cell$CTstrict[kept & per_cell$has_clone]),
      isimpson_myeloid_clones   = .inv_simpson(per_cell$CTstrict[kept & per_cell$has_clone]),
      shannon_lymph_baseline    = shannon_lymph,
      isimpson_lymph_baseline   = isimpson_lymph,
      mean_doublet_score_dropped = if (n_drop > 0L)
                                     round(mean(per_cell$doublet.score[dropped], na.rm = TRUE), 3)
                                   else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

# Tier 0.1 entry point.
run_ctaa_filter_diagnostic <- function(cfg) {
  paths_eye <- get_target_paths(cfg, "eye")
  eye_path  <- file.path(paths_eye$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(eye_path)) {
    log_message("CTaa diagnostic: eye object not found at ", eye_path, "; skipping.")
    return(invisible(FALSE))
  }
  rep_tab_dir <- file.path(cfg$paths$results_tables, "repertoire")
  ensure_dir(rep_tab_dir)

  log_message("CTaa diagnostic: loading eye object: ", eye_path)
  eye_obj <- readRDS(eye_path)

  bundle <- .build_myeloid_ctaa_table(eye_obj, cfg)
  per_cell <- bundle$per_cell
  lymph_tbl <- bundle$lymph_tbl

  per_cell_path <- file.path(rep_tab_dir, "myeloid_ctaa_per_cell.csv")
  write.csv(per_cell, per_cell_path, row.names = FALSE)
  log_message("Wrote ", per_cell_path)

  summary_df <- .summarize_filter_comparison(per_cell, lymph_tbl)
  summary_path <- file.path(rep_tab_dir, "myeloid_ctaa_filter_comparison.csv")
  write.csv(summary_df, summary_path, row.names = FALSE)
  log_message("Wrote ", summary_path)

  invisible(list(per_cell = per_cell, summary = summary_df))
}

# -----------------------------------------------------------------------------
# Tier 2.3: clone-definition sensitivity for eye<->blood overlap
# -----------------------------------------------------------------------------
#
# Recompute the eye-blood overlap counts under CTaa / CTnt / CTgene / CTstrict.
# Runs on the eye object metadata directly so the answer is independent of
# the compartment subset decisions.

run_clone_definition_sensitivity <- function(cfg) {
  # Eye<->blood overlap requires both tissues, so read the full integrated
  # object, not the eye sub-atlas (which has zero Blood cells by construction).
  full_path <- file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(full_path)) {
    log_message("Clone-definition sensitivity: full object missing at ",
                full_path, "; skipping.")
    return(invisible(FALSE))
  }
  rep_tab_dir <- file.path(cfg$paths$results_tables, "repertoire")
  ensure_dir(rep_tab_dir)

  obj <- readRDS(full_path)
  meta <- obj[[]]
  if (!all(c("Subject", "Tissue_1") %in% colnames(meta))) {
    log_message("Clone-definition sensitivity: Subject/Tissue_1 missing; skipping.")
    return(invisible(FALSE))
  }

  defs <- intersect(c("CTaa", "CTnt", "CTgene", "CTstrict"), colnames(meta))
  if (length(defs) == 0L) {
    log_message("Clone-definition sensitivity: no clone columns found; skipping.")
    return(invisible(FALSE))
  }

  meta$clone_class <- dplyr::case_when(
    grepl("^IG", meta$CTgene) ~ "BCR",
    grepl("^TR", meta$CTgene) ~ "TCR",
    TRUE                      ~ NA_character_
  )

  out_rows <- list()
  for (modality in c("BCR", "TCR")) {
    sub <- dplyr::filter(meta, clone_class == modality,
                         !is.na(Subject), !is.na(Tissue_1))
    if (nrow(sub) == 0L) next
    for (def in defs) {
      eye_clones   <- unique(stats::na.omit(sub[[def]][sub$Tissue_1 == "Eye"]))
      blood_clones <- unique(stats::na.omit(sub[[def]][sub$Tissue_1 == "Blood"]))
      shared       <- intersect(eye_clones, blood_clones)
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        modality         = modality,
        clone_definition = def,
        n_eye_clones     = length(eye_clones),
        n_blood_clones   = length(blood_clones),
        n_shared         = length(shared),
        pct_eye_shared   = if (length(eye_clones))   round(100 * length(shared) / length(eye_clones), 2)   else NA_real_,
        pct_blood_shared = if (length(blood_clones)) round(100 * length(shared) / length(blood_clones), 2) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out_rows) == 0L) {
    log_message("Clone-definition sensitivity: no clones found; skipping write.")
    return(invisible(FALSE))
  }
  out_df <- dplyr::bind_rows(out_rows)
  for (m in unique(out_df$modality)) {
    sub_m <- dplyr::filter(out_df, modality == m)
    out_path <- file.path(rep_tab_dir,
                          paste0(tolower(m), "_clone_definition_sensitivity.csv"))
    write.csv(sub_m, out_path, row.names = FALSE)
    log_message("Wrote ", out_path)
  }
  invisible(out_df)
}

# -----------------------------------------------------------------------------
# Tier 1.3: apply-time CTaa filter for the myeloid compartment
# -----------------------------------------------------------------------------
#
# Replaces the hardcoded "drop any CTaa+" logic that previously lived in
# subset_one_compartment(). Now driven by cfg$compartments$myeloid_ctaa_filter:
#   "a" — drop any cell with a productive CTaa (legacy default, conservative)
#   "b" — drop CTaa+ cells whose CTstrict matches a top-N expanded lymphocyte
#         clone (suggests doublet/carryover). Top-N is computed per Subject
#         over cells whose celltype_broad/cluster places them in the bcell
#         or tcell compartment.
#   "c" — keep all cells; tag CTaa+ in metadata as `ctaa_positive_myeloid`
#         so downstream tables can stratify on it.
#   "d" — drop CTaa+ AND lineage_score above median (T or B panel).
#
# Returns a list with:
#   $cells   — barcode vector to retain in the myeloid compartment
#   $tagged  — named logical vector flagging CTaa+ status (for option "c"
#              audit; NULL otherwise)
#   $report  — data.frame summarizing kept/dropped counts written to
#              outputs/tables/eye/myeloid/myeloid_ctaa_filter_apply.csv
apply_myeloid_ctaa_filter <- function(eye_obj, cells, cfg) {
  flt <- cfg$compartments$myeloid_ctaa_filter %||% "a"
  flt <- match.arg(as.character(flt), choices = c("a", "b", "c", "d"))
  meta <- eye_obj@meta.data

  if (!("CTaa" %in% colnames(meta))) {
    log_message("CTaa column missing; myeloid CTaa filter skipped.")
    return(list(cells = cells, tagged = NULL, report = NULL))
  }

  has_clone_all <- !is.na(meta$CTaa) & nzchar(as.character(meta$CTaa))
  has_clone     <- has_clone_all[cells]

  decision <- switch(flt,
    a = {
      keep <- !has_clone
      list(keep = keep, note = "drop any CTaa+ (legacy)")
    },
    b = {
      bkeys <- as.character(cfg$compartments$bcell$celltype_keys)
      tkeys <- as.character(cfg$compartments$tcell$celltype_keys)
      cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                       "knn.leiden.cluster" else "seurat_clusters"
      ctb_col <- resolve_celltype_broad(meta)
      lymph_mask <- if (length(c(bkeys, tkeys)) > 0L && !is.null(ctb_col)) {
        as.character(meta[[ctb_col]]) %in% c(bkeys, tkeys)
      } else {
        parents <- c(as.character(cfg$compartments$bcell$parent_clusters),
                     as.character(cfg$compartments$tcell$parent_clusters))
        as.character(meta[[cluster_col]]) %in% parents
      }
      ct_full   <- meta$CTstrict
      subj_full <- if ("Subject" %in% colnames(meta)) as.character(meta$Subject) else NA_character_
      top_n <- cfg$repertoire$top_n_lymph_for_filter %||% 10L

      top_df <- data.frame(
        Subject  = subj_full[lymph_mask & has_clone_all],
        CTstrict = ct_full[lymph_mask & has_clone_all],
        stringsAsFactors = FALSE
      ) %>%
        dplyr::filter(!is.na(CTstrict), nzchar(CTstrict)) %>%
        dplyr::count(Subject, CTstrict, sort = TRUE) %>%
        dplyr::group_by(Subject) %>%
        dplyr::slice_max(n, n = top_n, with_ties = FALSE) %>%
        dplyr::ungroup()
      top_key  <- paste(top_df$Subject, top_df$CTstrict, sep = "||")
      cell_key <- paste(subj_full[match(cells, rownames(meta))],
                        ct_full[match(cells, rownames(meta))],
                        sep = "||")
      is_top   <- cell_key %in% top_key & has_clone
      keep <- !is_top
      list(keep = keep, note = sprintf("drop CTaa+ matching top-%d lymphocyte clones", top_n))
    },
    c = list(keep = rep(TRUE, length(cells)), note = "keep all, flag in metadata"),
    d = {
      # Median lineage score on the myeloid cell pool.
      sub <- subset(eye_obj, cells = cells)
      DefaultAssay(sub) <- "RNA"
      sub <- NormalizeData(sub, verbose = FALSE)
      score_T <- .score_lineage_panels(eye_obj, cells, cfg) # reuse diagnostic helper
      sT <- score_T$score_T; sB <- score_T$score_B
      sT <- sT[match(cells, score_T$cell_id)]
      sB <- sB[match(cells, score_T$cell_id)]
      high <- (sT > stats::median(sT, na.rm = TRUE)) |
              (sB > stats::median(sB, na.rm = TRUE))
      high[is.na(high)] <- FALSE
      list(keep = !(has_clone & high),
           note = "drop CTaa+ AND lineage score > median")
    }
  )

  n_total <- length(cells)
  n_keep  <- sum(decision$keep)
  n_drop  <- n_total - n_keep
  log_message(sprintf(
    "Myeloid CTaa filter (%s — %s): keeping %d / %d cells (dropped %d)",
    flt, decision$note, n_keep, n_total, n_drop))

  report <- data.frame(
    filter_option = flt,
    description   = decision$note,
    n_total       = n_total,
    n_kept        = n_keep,
    n_dropped     = n_drop,
    stringsAsFactors = FALSE
  )

  tagged <- NULL
  if (flt == "c") {
    tagged <- has_clone
    names(tagged) <- cells
  }

  paths_myeloid <- get_target_paths(cfg, "myeloid")
  ensure_dir(paths_myeloid$results_tables)
  write.csv(report,
            file.path(paths_myeloid$results_tables, "myeloid_ctaa_filter_apply.csv"),
            row.names = FALSE)

  list(cells = cells[decision$keep], tagged = tagged, report = report)
}
