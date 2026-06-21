# R/32_escape.R
suppressPackageStartupMessages({
  library(Seurat)
  library(escape)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Helper: align a chunk matrix to a canonical row set ---
# Zero-fills any missing gene sets, reorders to canonical order
.align_chunk_rows <- function(mat, canonical_rows) {
  missing <- setdiff(canonical_rows, rownames(mat))
  if (length(missing) > 0) {
    pad <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                  dimnames = list(missing, colnames(mat)))
    mat <- rbind(mat, pad)
  }
  mat[canonical_rows, , drop = FALSE]
}

# --- Helper: load specific pathway rows for ALL cells (for visualization) ---
# Loads one chunk at a time, keeps only the requested rows → peak ~1.2 GB
.load_escape_features <- function(manifest, features, type = "unnorm") {
  parts <- list()
  for (i in seq_len(manifest$n_chunks)) {
    chunk <- readRDS(file.path(manifest$chunk_dir, paste0(type, "_", i, ".rds")))
    sub   <- .align_chunk_rows(chunk, features)[features, , drop = FALSE]
    parts[[length(parts) + 1]] <- sub
    rm(chunk, sub); gc(verbose = FALSE)
  }
  do.call(cbind, parts)
}

# --- Helper: list all available gene-set names from the first chunk ---
.escape_feature_names <- function(manifest, type = "unnorm") {
  chunk1 <- readRDS(file.path(manifest$chunk_dir, paste0(type, "_1.rds")))
  rownames(chunk1)
}

# --- Helper: load escape scores for a set of cells from chunk files ---
.load_escape_cells <- function(manifest, cells_wanted, type = "unnorm") {
  all_cells  <- manifest$cells
  chunk_size <- manifest$chunk_size
  chunk_dir  <- manifest$chunk_dir

  cell_pos     <- match(cells_wanted, all_cells)
  cell_pos     <- cell_pos[!is.na(cell_pos)]
  chunk_ids    <- unique(ceiling(cell_pos / chunk_size))

  # Canonical row set from chunk 1 — all chunks aligned to this
  canonical_rows <- .escape_feature_names(manifest, type)

  parts <- list()
  for (ci in sort(chunk_ids)) {
    chunk <- readRDS(file.path(chunk_dir, paste0(type, "_", ci, ".rds")))
    chunk_start  <- (ci - 1L) * chunk_size + 1L
    chunk_cells  <- all_cells[chunk_start:(chunk_start + ncol(chunk) - 1L)]
    keep         <- intersect(cells_wanted, chunk_cells)
    if (length(keep) > 0) {
      sub <- .align_chunk_rows(chunk[, keep, drop = FALSE], canonical_rows)
      parts[[length(parts) + 1]] <- sub
    }
    rm(chunk); gc(verbose = FALSE)
  }
  do.call(cbind, parts)
}

# --- Run ssGSEA scoring ---
run_escape_ssgsea <- function(cfg) {
  if (!isTRUE(cfg$escape$enable)) {
    message("ESCAPE step disabled (config.escape.enable = FALSE). Skipping.")
    return(invisible(TRUE))
  }

  obj_path <- file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    warning("obj object not found at: ", obj_path, " — skipping ESCAPE.")
    return(invisible(TRUE))
  }

  # Check if scoring is already complete (manifest exists with all chunks)
  manifest_path <- file.path(cfg$paths$results_objects, "EscapeChunkManifest.rds")
  if (file.exists(manifest_path)) {
    manifest <- readRDS(manifest_path)
    all_present <- all(file.exists(
      file.path(manifest$chunk_dir,
                c(paste0("unnorm_", seq_len(manifest$n_chunks), ".rds"),
                  paste0("norm_",   seq_len(manifest$n_chunks), ".rds")))
    ))
    if (all_present) {
      message("ESCAPE scoring already complete (", manifest$n_chunks,
              " chunks on disk). Skipping.")
      return(invisible(TRUE))
    }
  }

  obj <- readRDS(obj_path)

  species   <- if (!is.null(cfg$escape$species)) cfg$escape$species else "Homo sapiens"
  libraries <- if (!is.null(cfg$escape$libraries)) unlist(cfg$escape$libraries) else c("H","C2")
  method    <- if (!is.null(cfg$escape$method)) cfg$escape$method else "UCell"
  groups    <- if (!is.null(cfg$escape$groups)) as.integer(cfg$escape$groups) else 10000L
  min_size  <- if (!is.null(cfg$escape$min_size)) as.integer(cfg$escape$min_size) else 3L
  assay_out <- if (!is.null(cfg$escape$new_assay)) cfg$escape$new_assay else "escape.UCell"

  message("Fetching gene sets via escape::getGeneSets() …")
  gs <- escape::getGeneSets(species = species, library = libraries)
  gs.length <- unlist(lapply(gs, length))
  gs <- gs[-which(gs.length > 1500)]

  #Removing Cancer-Related Gene Sets and Reactome
  patterns <- "REACTOME|CLL|GLIOMA|NEUROBLASTOMA|ONCOGENESIS|GLIOBLASTOMA|MEDULLOBLASTOMA|IMMORTALIZED|MELANOMA|ONCOGENIC|CARCINOMA|AML|LEUKEMIA|SARCOMA|ASTROCYTOMA|HEPATOBLASTOMA|MYELODYSPLASTIC|CANCER|NEOPLASTIC|SEMINOMA|FUSION|CARCINOGENESIS|RHABDOMYOSARCOMA|MYELOMA|ADENOMA|MESOTHELIOMA|MESOTELIOMA|LYMPHOMA|PAPILLOMA|MYCOSIS|SEZARY|TUMOR|METASTASIS"
  remove.idx <- c(which(is.na(names(gs))), grep(patterns, names(gs)))
  gs <- gs[-remove.idx]

  # Ensure we use an expression assay; ESCAPE expects a counts or data matrix
  assay_to_use <- if ("RNA" %in% names(obj@assays)) "RNA" else "SCT"
  DefaultAssay(obj) <- assay_to_use

  # --- Chunked ESCAPE: score + normalize per chunk, save to disk ---
  # Never assembles the full matrix — avoids the ~24 GB peak that aborts R.
  cells_all  <- colnames(obj)
  n_cells    <- length(cells_all)
  chunk_size <- 50000L
  n_chunks   <- ceiling(n_cells / chunk_size)
  chunk_dir  <- file.path(cfg$paths$results_objects, "escape_chunks")
  dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)

  sf_all <- obj$nFeature_RNA
  sf_all[is.na(sf_all) | sf_all == 0] <- 1

  message("Running ESCAPE in ", n_chunks, " chunks of ", chunk_size, " cells …")

  for (i in seq_len(n_chunks)) {
    unnorm_f <- file.path(chunk_dir, paste0("unnorm_", i, ".rds"))
    norm_f   <- file.path(chunk_dir, paste0("norm_", i, ".rds"))

    if (file.exists(unnorm_f) && file.exists(norm_f)) {
      message("  Chunk ", i, "/", n_chunks, " — already on disk, skipping.")
      next
    }

    idx_start <- (i - 1L) * chunk_size + 1L
    idx_end   <- min(i * chunk_size, n_cells)
    cell_ids  <- cells_all[idx_start:idx_end]

    message("  Chunk ", i, "/", n_chunks, " (", length(cell_ids), " cells) …")
    obj_chunk <- subset(obj, cells = cell_ids)

    obj_chunk <- escape::runEscape(
      obj_chunk,
      method    = method,
      gene.sets = gs,
      groups    = groups,
      min.size  = min_size,
      new.assay.name = assay_out
    )

    chunk_mat <- LayerData(obj_chunk[[assay_out]], layer = "data")

    # Normalize this chunk immediately
    sf_chunk  <- sf_all[idx_start:idx_end]
    norm_mat  <- sweep(chunk_mat, 2, sf_chunk, FUN = "/")

    saveRDS(chunk_mat, unnorm_f)
    saveRDS(norm_mat,  norm_f)

    rm(obj_chunk, chunk_mat, norm_mat, sf_chunk); gc(verbose = FALSE)
  }

  # Save manifest — downstream functions load chunks on demand via .load_escape_cells()
  manifest <- list(
    n_chunks   = n_chunks,
    chunk_size = chunk_size,
    cells      = cells_all,
    chunk_dir  = chunk_dir
  )
  saveRDS(manifest, manifest_path)
  message("ESCAPE scoring complete. Manifest saved to: ", manifest_path)

  invisible(TRUE)
}

# --- Pathway-level differential testing ---
# Pseudobulk at the SAMPLE level (median pathway score per sample per celltype),
# then fit a limma linear model. Previously this ran FindMarkers at cell level,
# which is pseudoreplication: with ~n samples per group but ~10^4 cells, the
# p-values were wildly inflated.
#
# Design: if the contrast varies within Subject AND the `~ Subject + group`
# design is rank-sufficient (paired or partially-paired data), blocks by
# Subject. Falls back to `~ group` otherwise.
# Curated regex tags for compartment-specific gene-set axis grouping. Used by
# the F3-F5 GSEA bubble plot to color-code pathway groups (Type 17, ISG,
# antigen-presentation, exhaustion, etc.).
.escape_axis_regex <- list(
  myeloid = "IL17|IL23|TH17|CD1|LANGERHANS|MREGDC|LAM|INTERFERON|COMPLEMENT|TNFA_SIGNALING|INFLAMMATORY|ANTIGEN_PRESENTATION|MHC",
  bcell   = "MEMORY_B|PLASMA|GERMINAL_CENTER|FCRL4|ATYPICAL_B|CLASS_SWITCH|SOMATIC_HYPERMUTATION|INTERFERON",
  tcell   = "TH17|TH1|TH2|TFH|TREG|EXHAUSTION|TISSUE_RESIDENT|TRM|STEMNESS|EFFECTOR|CYTOTOXIC|IL2_STAT5|INTERFERON_GAMMA|TNFA_SIGNALING"
)

.tag_axis <- function(gene_sets, target) {
  rgx <- .escape_axis_regex[[target]]
  if (is.null(rgx)) return(rep(NA_character_, length(gene_sets)))
  m <- regmatches(gene_sets, regexpr(rgx, gene_sets))
  out <- rep(NA_character_, length(gene_sets))
  hits <- grepl(rgx, gene_sets)
  out[hits] <- m
  out
}

# Resolve an etiology -> bucket map from cfg$etiology_groups. Etiologies with
# fewer than min_samples_for_etiology distinct samples collapse into a pooled
# bucket ("<group>_pooled"). Returns named character vector: names = etiology,
# values = bucket label used as column in the pathway-by-etiology table.
.resolve_etiology_buckets <- function(meta, sample_col, cfg) {
  groups <- cfg$etiology_groups
  if (is.null(groups)) return(NULL)
  min_n  <- as.integer(groups$min_samples_for_etiology %||% 3L)
  niu    <- as.character(groups$niu     %||% character(0))
  viral  <- as.character(groups$viral   %||% character(0))
  healthy <- as.character(groups$healthy %||% character(0))

  if (!"Etiology" %in% colnames(meta)) return(NULL)

  # Per-etiology sample count
  samp_per_et <- meta |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c("Etiology", sample_col)))) |>
    dplyr::count(Etiology, name = "n_samples")
  n_lookup <- setNames(samp_per_et$n_samples, samp_per_et$Etiology)

  bucket_for <- function(et, parent_group) {
    n <- as.integer(n_lookup[et] %||% 0L)
    if (is.na(n) || n < min_n) paste0(parent_group, "_pooled") else et
  }
  out <- character(0)
  for (et in niu)     out[et] <- bucket_for(et, "NIU")
  for (et in viral)   out[et] <- bucket_for(et, "Viral")
  for (et in healthy) out[et] <- bucket_for(et, "Healthy")
  out
}

# Compute per-etiology median pathway scores (samples x pathways pseudobulk,
# medianed within etiology bucket). Returns a long data.frame:
#   pathway, etiology_bucket, parent_group, n_samples, median_score.
# Writes "Pathway_by_etiology_<target>.csv" under paths$results_tables.
#
# Within the NIU arm, an additional limma fit per NIU etiology (etiology vs
# other-NIU pooled) is run and saved as "Pathway_NIU_subcontrast_<target>.csv".
# This is the "what's specific to HLA_B27 vs other NIU" view requested in the
# F3-F5 redesign.
run_escape_etiology_breakdown <- function(cfg, target = c("myeloid", "bcell", "tcell")) {
  target <- match.arg(target)
  if (!isTRUE(cfg$escape$enable)) return(invisible(TRUE))

  suppressPackageStartupMessages({
    library(limma)
  })

  paths   <- get_target_paths(cfg, target)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  manifest_path <- file.path(cfg$paths$results_objects, "EscapeChunkManifest.rds")
  if (!file.exists(obj_path) || !file.exists(manifest_path)) {
    log_message("Etiology breakdown: missing inputs for ", target,
                " (need IntegratedSeuratObject.rds + EscapeChunkManifest.rds).")
    return(invisible(FALSE))
  }
  ensure_dir(paths$results_tables)

  manifest <- readRDS(manifest_path)
  obj      <- readRDS(obj_path)
  meta_all <- obj[[]]

  sample_col <- if ("Subject_Timepoint" %in% colnames(meta_all))
                  "Subject_Timepoint" else "orig.ident"

  bucket_map <- .resolve_etiology_buckets(meta_all, sample_col, cfg)
  if (is.null(bucket_map) || length(bucket_map) == 0) {
    log_message("Etiology breakdown: no etiology_groups config; skipping ", target)
    return(invisible(FALSE))
  }

  # Reuse the curated NIU-vs-Viral GSEA pathways as the row set so the heatmap
  # is comparable to the existing bubble plot. Fall back to the full pathway
  # set when the GSEA table is absent.
  gsea_path <- file.path(paths$results_tables, "GSEA_Autoimmune_vs_Viral.csv")
  keep_paths <- NULL
  if (file.exists(gsea_path)) {
    g <- read.csv(gsea_path, stringsAsFactors = FALSE)
    keep_paths <- unique(g$pathway[!is.na(g$axis_tag) & g$stratum == "global"])
    if (length(keep_paths) == 0)
      keep_paths <- unique(g$pathway[!is.na(g$axis_tag)])
  }

  # Cells we care about: those whose etiology maps to a bucket
  meta_all$Etiology <- as.character(meta_all$Etiology)
  meta_all$bucket   <- bucket_map[meta_all$Etiology]
  cells_keep <- rownames(meta_all)[!is.na(meta_all$bucket)]
  if (length(cells_keep) < 50) {
    log_message("Etiology breakdown: <50 cells with mapped etiology for ", target)
    return(invisible(FALSE))
  }

  log_message("Etiology breakdown (", target, "): loading UCell scores for ",
              length(cells_keep), " cells.")
  esc_mat <- .load_escape_cells(manifest, cells_keep, type = "unnorm")
  cells_keep <- intersect(cells_keep, colnames(esc_mat))
  esc_mat <- esc_mat[, cells_keep, drop = FALSE]
  if (!is.null(keep_paths)) {
    keep_paths <- intersect(keep_paths, rownames(esc_mat))
    if (length(keep_paths) >= 5) {
      esc_mat <- esc_mat[keep_paths, , drop = FALSE]
    }
  }

  cell_meta <- meta_all[cells_keep, , drop = FALSE]
  samp_vec  <- as.character(cell_meta[[sample_col]])
  samples   <- unique(samp_vec)

  # Sample-level median pathway matrix (pathway x sample)
  pb <- vapply(samples, function(s) {
    cols <- which(samp_vec == s)
    if (length(cols) < 3) rep(NA_real_, nrow(esc_mat))
    else matrixStats::rowMedians(as.matrix(esc_mat[, cols, drop = FALSE]),
                                 na.rm = TRUE)
  }, numeric(nrow(esc_mat)))
  rownames(pb) <- rownames(esc_mat)
  colnames(pb) <- samples
  pb <- pb[, !apply(pb, 2, function(x) all(is.na(x))), drop = FALSE]

  samp_df <- cell_meta |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c(sample_col, "Etiology",
                                                  "bucket", "Phenotype_2")))) |>
    as.data.frame()
  rownames(samp_df) <- as.character(samp_df[[sample_col]])
  samp_df <- samp_df[colnames(pb), , drop = FALSE]

  # ---- 1. Pathway x etiology-bucket median table ----------------------------
  long <- list()
  for (b in unique(samp_df$bucket)) {
    cols <- rownames(samp_df)[samp_df$bucket == b]
    if (length(cols) == 0) next
    med  <- matrixStats::rowMedians(pb[, cols, drop = FALSE], na.rm = TRUE)
    long[[b]] <- data.frame(
      pathway         = rownames(pb),
      etiology_bucket = b,
      parent_group    = unique(samp_df$Phenotype_2[samp_df$bucket == b])[1],
      n_samples       = length(cols),
      median_score    = med,
      stringsAsFactors = FALSE
    )
  }
  out_long <- do.call(rbind, long)
  out_path <- file.path(paths$results_tables,
                        paste0("Pathway_by_etiology_", target, ".csv"))
  write.csv(out_long, out_path, row.names = FALSE)
  log_message("  Saved: ", basename(out_path), " (", nrow(out_long), " rows, ",
              length(unique(out_long$etiology_bucket)), " etiology buckets)")

  # ---- 2. Per-NIU-etiology sub-contrast (etiology vs other-NIU pooled) -----
  niu_set <- intersect(names(bucket_map)[grepl("^NIU", bucket_map) |
                                          bucket_map %in% cfg$etiology_groups$niu],
                       cfg$etiology_groups$niu)
  niu_samp <- samp_df[samp_df$Phenotype_2 == "NIU", , drop = FALSE]
  niu_results <- list()
  for (et in unique(niu_samp$Etiology)) {
    n_in  <- sum(niu_samp$Etiology == et)
    n_out <- sum(niu_samp$Etiology != et)
    if (n_in < 2 || n_out < 2) next
    samp_sub <- niu_samp
    samp_sub$group <- factor(ifelse(samp_sub$Etiology == et, et, "other_NIU"),
                             levels = c("other_NIU", et))
    pb_sub <- pb[, rownames(samp_sub), drop = FALSE]
    keep_rows <- apply(pb_sub, 1, function(x)
      sum(!is.na(x)) >= 4 && stats::var(x, na.rm = TRUE) > 0)
    pb_sub <- pb_sub[keep_rows, , drop = FALSE]
    if (nrow(pb_sub) < 10) next

    design <- model.matrix(~ group, data = samp_sub)
    fit <- limma::eBayes(limma::lmFit(pb_sub, design), robust = TRUE)
    tab <- limma::topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
    tab$pathway      <- rownames(pb_sub)
    tab$etiology     <- et
    tab$n_in         <- n_in
    tab$n_other_NIU  <- n_out
    niu_results[[et]] <- tab
  }
  if (length(niu_results) > 0) {
    out_niu <- dplyr::bind_rows(niu_results) |>
      dplyr::rename(logFC = logFC, p_val = P.Value) |>
      dplyr::mutate(FDR = stats::p.adjust(p_val, method = "BH"),
                    direction = ifelse(logFC > 0,
                                       paste0(etiology, "_up"),
                                       "other_NIU_up"))
    out_path_niu <- file.path(paths$results_tables,
                              paste0("Pathway_NIU_subcontrast_", target, ".csv"))
    write.csv(out_niu, out_path_niu, row.names = FALSE)
    log_message("  Saved: ", basename(out_path_niu), " (",
                length(unique(out_niu$etiology)),
                " NIU etiologies tested vs other-NIU pooled)")
  } else {
    log_message("  No NIU etiology had >=2 samples + other-NIU for sub-contrast.")
  }

  invisible(TRUE)
}

# run_escape_etiology_breakdown is a compartment-only follow-up to
# run_escape_differential. For each NIU/Viral etiology bucket it computes a
# sample-level median UCell score per curated pathway (heatmap input for the
# F3-F5 panels) and, within the NIU arm, fits a limma model of each NIU
# etiology vs all other NIU subjects pooled (HLA_B27 vs other-NIU etc.).
# Small etiology groups (n < cfg$etiology_groups$min_samples_for_etiology) are
# collapsed into "<group>_pooled" so the heatmap doesn't have n=1 columns.
# Outputs go to outputs/tables/eye/<cmp>/Pathway_by_etiology_<cmp>.csv and
# Pathway_NIU_subcontrast_<cmp>.csv.

# Score a hand-curated list of gene-set modules with UCell on a compartment
# object and stash the result in a new Seurat assay (e.g. escape.APC on the
# myeloid object). Reads modules from cfg$escape$custom_modules[[set_key]]
# and writes back to the compartment IntegratedSeuratObject.rds so downstream
# panel-F viz can pull module scores by GetAssayData(obj, assay = new_assay).
#
# This is independent of the genome-wide ssGSEA chunk manifest: the assay is
# small (one row per module), per-cell, and computed in-memory. UCell ranks
# are computed within the compartment object so the scores are directly
# comparable to escape.UCell_eye on the same cells.
run_escape_custom_modules <- function(cfg,
                                      target  = c("myeloid", "bcell", "tcell"),
                                      set_key = NULL) {
  target <- match.arg(target)
  if (!isTRUE(cfg$escape$enable)) return(invisible(TRUE))

  sets <- cfg$escape$custom_modules
  if (is.null(sets) || length(sets) == 0) {
    log_message("escape$custom_modules not set; skipping custom modules.")
    return(invisible(TRUE))
  }
  if (is.null(set_key)) {
    matches <- vapply(sets, function(s) {
      isTRUE(s$target_compartment == target)
    }, logical(1))
    set_keys <- names(sets)[matches]
  } else {
    set_keys <- set_key
  }
  if (length(set_keys) == 0) {
    log_message("No custom module set targets compartment '", target,
                "'. Skipping.")
    return(invisible(TRUE))
  }
  if (!requireNamespace("escape", quietly = TRUE) ||
      !requireNamespace("GSEABase", quietly = TRUE)) {
    log_message("escape / GSEABase not installed; skipping custom modules.")
    return(invisible(FALSE))
  }

  paths    <- get_target_paths(cfg, target)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Custom modules: ", obj_path, " missing. Skipping.")
    return(invisible(FALSE))
  }
  obj <- readRDS(obj_path)

  for (key in set_keys) {
    spec <- sets[[key]]
    if (is.null(spec$modules) || length(spec$modules) == 0) next
    assay_name <- spec$new_assay %||% paste0("escape.custom.", key)
    log_message("Custom modules: scoring '", key, "' -> assay '", assay_name,
                "' on ", target)

    # Hand-curated modules are small (often 4-12 genes), well below the
    # default escape::runEscape(min.size = 5) cutoff. Take the lower of the
    # config-specified `min.size` and the smallest module size so a single
    # 4-gene module like COSTIM_LICENSING doesn't get silently dropped.
    detected <- lapply(spec$modules, function(g) intersect(g, rownames(obj)))
    detected <- detected[vapply(detected, length, integer(1)) >= 2]
    if (length(detected) == 0) {
      log_message("  No modules had >=2 genes detected on ", target,
                  "; skipping set ", key)
      next
    }
    gs  <- lapply(names(detected), function(nm)
      GSEABase::GeneSet(detected[[nm]], setName = nm))
    gsc <- GSEABase::GeneSetCollection(gs)
    user_min_size <- spec$min_size %||% 2L
    min_size <- max(2L,
                    min(as.integer(user_min_size),
                        min(vapply(detected, length, integer(1)))))
    log_message("  Scoring ", length(gsc), " modules with min.size = ",
                min_size)
    obj <- tryCatch(
      escape::runEscape(obj, method = "UCell", gene.sets = gsc,
                        min.size = min_size,
                        new.assay.name = assay_name,
                        BPPARAM = BiocParallel::SerialParam()),
      error = function(e) {
        log_message("  runEscape failed for ", key, ": ",
                    conditionMessage(e))
        obj
      })
  }
  saveRDS(obj, obj_path)
  log_message("Custom modules: saved compartment object with new assays to ",
              obj_path)
  invisible(TRUE)
}

run_escape_differential <- function(cfg, target = c("all", "eye", "myeloid", "bcell", "tcell")) {
  target <- match.arg(target)
  if (!isTRUE(cfg$escape$enable)) {
    log_message("ESCAPE disabled. Skipping differential.")
    return(invisible(TRUE))
  }

  suppressPackageStartupMessages({
    library(limma)
  })

  paths <- get_target_paths(cfg, target)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Integrated object not found at ", obj_path, ". Skipping ESCAPE differential.")
    return(invisible(TRUE))
  }

  # ESCAPE chunk manifest is generated for the full atlas; reuse it for all
  # targets. UCell scores are per-cell rank-based so subsetting yields
  # identical scores at the cells of interest.
  manifest_path <- file.path(cfg$paths$results_objects, "EscapeChunkManifest.rds")
  if (!file.exists(manifest_path)) {
    log_message("ESCAPE chunk manifest not found. Run run_escape_ssgsea() first.")
    return(invisible(TRUE))
  }

  manifest <- readRDS(manifest_path)
  obj      <- readRDS(obj_path)

  ensure_dir(paths$results_tables)
  ensure_dir(paths$viz_dir)

  # Sample identifier (prefer most granular)
  meta_all <- obj[[]]
  if ("Subject_Timepoint" %in% colnames(meta_all)) {
    sample_col <- "Subject_Timepoint"
  } else if ("orig.ident" %in% colnames(meta_all)) {
    sample_col <- "orig.ident"
  } else {
    sample_col <- "Subject"
  }
  has_subject  <- "Subject"        %in% colnames(meta_all)
  has_ct_broad <- "celltype_broad" %in% colnames(meta_all)
  is_compartment <- target %in% c("myeloid", "bcell", "tcell")

  # Compartment-scoped runs: only Autoimmune_vs_Viral, stratify by Leiden
  # substate rather than celltype_broad. Eye-only (target = "eye") shares the
  # full-atlas contrast list because the eye object spans both phenotypes.
  if (is_compartment) {
    escape_contrasts <- list(
      list(name = "Autoimmune_vs_Viral", col = "Phenotype_2", g1 = "NIU", g2 = "Viral")
    )
    stratify_col <- "knn.leiden.cluster"
  } else {
    escape_contrasts <- list(
      list(name = "Eye_vs_Blood",        col = "Tissue_1",    g1 = "Eye", g2 = "Blood"),
      list(name = "Autoimmune_vs_Viral", col = "Phenotype_2", g1 = "NIU", g2 = "Viral")
    )
    stratify_col <- "celltype_broad"
  }

  # --- Per-contrast/per-celltype sample-level test -------------------------
  for (contrast in escape_contrasts) {
    log_message("ESCAPE differential (sample-level): ", contrast$name)

    tryCatch({
      cells_keep <- rownames(meta_all)[
        meta_all[[contrast$col]] %in% c(contrast$g1, contrast$g2)
      ]
      if (length(cells_keep) < 50) {
        log_message("  Too few cells; skipping.")
        next
      }

      log_message("  Loading escape scores for ", length(cells_keep), " cells …")
      esc_mat <- .load_escape_cells(manifest, cells_keep, type = "unnorm")
      cells_keep <- intersect(cells_keep, colnames(esc_mat))
      esc_mat <- esc_mat[, cells_keep, drop = FALSE]

      cell_meta <- meta_all[cells_keep, , drop = FALSE]

      # Build stratum list: "global" + each substate/celltype_broad with enough cells.
      strata <- "global"
      if (stratify_col %in% colnames(cell_meta)) {
        strat_vec <- as.character(cell_meta[[stratify_col]])
        ct_tab    <- table(strat_vec[!is.na(strat_vec)])
        strata    <- c(strata, names(ct_tab)[ct_tab >= 50])
      }

      contrast_results <- list()
      for (stratum in strata) {
        stratum_cells <- if (stratum == "global") cells_keep
                         else cells_keep[as.character(cell_meta[[stratify_col]]) == stratum &
                                         !is.na(cell_meta[[stratify_col]])]
        if (length(stratum_cells) < 30) next

        # Aggregate per-sample median pathway score
        sub_meta <- cell_meta[stratum_cells, , drop = FALSE]
        samp_vec <- as.character(sub_meta[[sample_col]])
        samples  <- unique(samp_vec)

        # Build samples x pathways matrix (rows = pathway, cols = sample)
        pb <- vapply(samples, function(s) {
          cols <- which(samp_vec == s)
          if (length(cols) < 3) {
            rep(NA_real_, nrow(esc_mat))
          } else {
            matrixStats::rowMedians(
              as.matrix(esc_mat[, stratum_cells[cols], drop = FALSE]),
              na.rm = TRUE)
          }
        }, numeric(nrow(esc_mat)))
        rownames(pb) <- rownames(esc_mat)
        colnames(pb) <- samples
        pb <- pb[, !apply(pb, 2, function(x) all(is.na(x))), drop = FALSE]
        if (ncol(pb) < 4) next

        # Sample-level metadata
        samp_df <- sub_meta %>%
          dplyr::distinct(!!sym(sample_col), !!sym(contrast$col),
                          !!!(if (has_subject) rlang::syms("Subject") else list())) %>%
          as.data.frame()
        rownames(samp_df) <- as.character(samp_df[[sample_col]])
        samp_df <- samp_df[colnames(pb), , drop = FALSE]
        samp_df$group <- factor(as.character(samp_df[[contrast$col]]),
                                levels = c(contrast$g2, contrast$g1))
        if (length(unique(samp_df$group)) < 2 || any(table(samp_df$group) < 2)) next

        # Paired design when feasible
        use_paired <- FALSE
        if (has_subject) {
          samp_df$Subject <- factor(as.character(samp_df$Subject))
          per_sbj <- tapply(as.character(samp_df$group),
                            as.character(samp_df$Subject),
                            function(x) length(unique(x)))
          if (any(per_sbj > 1, na.rm = TRUE)) {
            mm_try <- try(model.matrix(~ Subject + group, data = samp_df),
                          silent = TRUE)
            if (!inherits(mm_try, "try-error") &&
                qr(mm_try)$rank == ncol(mm_try) &&
                nrow(mm_try) > ncol(mm_try)) {
              use_paired <- TRUE
            }
          }
        }
        design_f <- if (use_paired) ~ Subject + group else ~ group
        design   <- model.matrix(design_f, data = samp_df)
        coef_name <- tail(colnames(design), 1)

        # Drop pathways with NA variance (constant within samples)
        keep_rows <- apply(pb, 1, function(x) sum(!is.na(x)) >= 4 &&
                                               stats::var(x, na.rm = TRUE) > 0)
        pb <- pb[keep_rows, , drop = FALSE]
        if (nrow(pb) < 10) next

        fit <- limma::lmFit(pb, design)
        fit <- limma::eBayes(fit, robust = TRUE)
        tab <- limma::topTable(fit, coef = coef_name, number = Inf,
                               sort.by = "none")
        tab$pathway  <- rownames(pb)
        tab$stratum  <- stratum
        tab$n_g1     <- sum(samp_df$group == contrast$g1)
        tab$n_g2     <- sum(samp_df$group == contrast$g2)
        tab$paired   <- use_paired
        contrast_results[[stratum]] <- tab

        log_message(sprintf(
          "    [%s] stratum %-15s pathways=%d samples=%d paired=%s",
          contrast$name, stratum, nrow(pb), ncol(pb), use_paired))
      }
      rm(esc_mat); gc(verbose = FALSE)

      if (length(contrast_results) == 0) {
        log_message("  No strata produced results for ", contrast$name)
        next
      }

      out <- dplyr::bind_rows(contrast_results) %>%
        dplyr::rename(logFC = logFC, p_val = P.Value, p_adj_stratum = adj.P.Val) %>%
        dplyr::mutate(contrast = contrast$name)
      # Global FDR across all tests in this contrast (all strata pooled)
      out$FDR <- stats::p.adjust(out$p_val, method = "BH")
      # Backward-compat alias column (legacy consumers expect p_val_adj)
      out$p_val_adj <- out$FDR

      # Compartment runs add an axis_tag column for the F3-F5 GSEA bubble plot.
      if (is_compartment) {
        out$axis_tag <- .tag_axis(out$pathway, target)
        out$direction <- ifelse(out$logFC > 0,
                                paste0(contrast$g1, "_up"),
                                paste0(contrast$g2, "_up"))
        # Compartment naming convention: GSEA_<contrast>.csv per design doc.
        out_path <- file.path(paths$results_tables,
                              paste0("GSEA_", contrast$name, ".csv"))
      } else {
        out_path <- file.path(paths$results_tables,
                              paste0("ESCAPE_diff_", contrast$name, ".csv"))
      }
      write.csv(out, out_path, row.names = FALSE)
      log_message("  Saved: ", basename(out_path),
                  " (", nrow(out), " pathway x stratum tests)")

    }, error = function(e) {
      log_message("  ESCAPE differential failed for ", contrast$name, ": ",
                  conditionMessage(e))
    })
  }

  invisible(TRUE)
}
