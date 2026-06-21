# R/17_lineage_validation.R
#
# T0.2 from the 2026-05-13 plan. Score every eye Leiden cluster against a
# curated marker panel and compare the top-scoring panel to the existing
# celltype call. The panels themselves are tissue-aware (microglia, TRM,
# mregDC) so cells mis-called by blood-biased references like Azimuth PBMC
# get surfaced before compartment subsetting.
#
# Read-only with respect to the eye object — outputs are:
#   outputs/tables/eye/cluster_lineage_validation.csv
#       cluster x panel mean module score, plus existing celltype_broad,
#       top-scoring panel, agreement boolean, and a `flag` column noting
#       any disagreement that warrants a `cfg$eye_manual_overrides` entry.
#
# The marker panels follow audit R1 verbatim plus a Myeloid catch-all so
# every compartment parent gets a positive panel. Tune the lineage_keys
# mapping below if your celltype_broad labels differ from the defaults
# in cfg$eye_annotation.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Marker panels from the audit (R1). The `lineage_keys` mapping records
# which celltype_broad labels are *expected* to top-score on each panel,
# so disagreement = top-scoring panel implies a different lineage than
# the call. Adjust the keys here when celltype_broad labels change.
.lineage_panels <- function(cfg) {
  list(
    Microglia = c("TMEM119", "P2RY12", "CX3CR1", "AIF1", "CSF1R"),
    TRM       = c("ITGAE", "CD69", "CXCR6", "ZNF683", "RUNX3"),
    mregDC    = c("CD1A", "CD1C", "LAMP3", "CCR7", "CD83", "FSCN1"),
    Bplasma   = c("CD19", "MS4A1", "CD79A", "SDC1", "XBP1", "MZB1"),
    Tcell     = c("CD3D", "CD3E", "CD8A", "CD4"),
    Lens      = cfg$eye_qc$lens_filter$genes %||% character(0),
    Myeloid   = c("LYZ", "CD14", "FCGR3A", "C1QA", "C1QB", "APOE")
  )
}

.lineage_expected_broad <- function() {
  # panel -> celltype_broad values that should top-score on it.
  list(
    Microglia = c("Macrophage", "Mono", "Mac", "Monocyte"),
    TRM       = c("CD8_T", "CD4_T", "T cell", "T"),
    mregDC    = c("DC", "mregDC", "cDC"),
    Bplasma   = c("B cell", "Plasma", "B", "Plasma_cell", "Memory_B"),
    Tcell     = c("CD4_T", "CD8_T", "T cell", "T", "gd_T"),
    Lens      = c("Lens"),
    Myeloid   = c("Monocyte", "Macrophage", "Mac", "Mono")
  )
}

# Score one panel via AddModuleScore. Returns per-cell numeric vector
# (names = barcodes). Falls back to log-norm colSums for tiny panels
# where AddModuleScore's control selection fails.
.score_one_panel <- function(obj, genes, name, seed = 42L) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) < 2L) {
    return(stats::setNames(rep(NA_real_, ncol(obj)), colnames(obj)))
  }
  score_name <- paste0("__panel_", name)
  scored <- tryCatch(
    AddModuleScore(obj, features = list(genes), name = score_name,
                   seed = seed, ctrl = min(50L, max(5L, length(genes)))),
    error = function(e) NULL
  )
  if (is.null(scored)) {
    mat <- GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
    return(stats::setNames(Matrix::colSums(mat), colnames(obj)))
  }
  scored[[paste0(score_name, "1"), drop = TRUE]]
}

# Top-level entry point.
run_lineage_validation <- function(cfg) {
  paths_eye <- get_target_paths(cfg, "eye")
  eye_path  <- file.path(paths_eye$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(eye_path)) {
    log_message("Lineage validation: eye object not found at ", eye_path, "; skipping.")
    return(invisible(FALSE))
  }
  log_message("Lineage validation: loading eye object")
  eye_obj <- readRDS(eye_path)

  DefaultAssay(eye_obj) <- "RNA"
  eye_obj <- NormalizeData(eye_obj, verbose = FALSE)

  panels <- .lineage_panels(cfg)
  panels <- panels[lengths(panels) > 0L]
  log_message("Scoring ", length(panels), " panels")

  scores_df <- data.frame(cell_id = colnames(eye_obj), stringsAsFactors = FALSE)
  for (pn in names(panels)) {
    scores_df[[paste0("score_", pn)]] <-
      .score_one_panel(eye_obj, panels[[pn]], pn, seed = cfg$seed %||% 42L)
  }

  meta <- eye_obj[[]]
  cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                   "knn.leiden.cluster" else "seurat_clusters"
  scores_df$eye_cluster <- as.character(meta[[cluster_col]])
  ctb_col <- resolve_celltype_broad(meta)
  scores_df$celltype_broad <- if (!is.null(ctb_col)) as.character(meta[[ctb_col]])
                              else                    NA_character_

  # Aggregate per cluster: mean module score per panel + dominant celltype_broad.
  per_cluster <- scores_df %>%
    dplyr::group_by(eye_cluster) %>%
    dplyr::summarise(
      n_cells           = dplyr::n(),
      celltype_broad    = {
        tbl <- table(celltype_broad, useNA = "ifany")
        if (length(tbl) == 0L) NA_character_ else names(which.max(tbl))[1]
      },
      dplyr::across(dplyr::starts_with("score_"), \(x) round(mean(x, na.rm = TRUE), 4)),
      .groups = "drop"
    )

  # Identify the top-scoring panel per cluster.
  score_cols <- grep("^score_", colnames(per_cluster), value = TRUE)
  per_cluster$top_panel <- apply(per_cluster[, score_cols, drop = FALSE], 1, function(r) {
    if (all(is.na(r))) return(NA_character_)
    sub("^score_", "", score_cols[which.max(r)])
  })
  per_cluster$top_panel_score <- apply(per_cluster[, score_cols, drop = FALSE], 1, function(r)
    if (all(is.na(r))) NA_real_ else round(max(r, na.rm = TRUE), 4))

  # Agreement check.
  expected <- .lineage_expected_broad()
  per_cluster$agrees <- mapply(function(panel, broad) {
    if (is.na(panel) || is.na(broad)) return(NA)
    exp_set <- expected[[panel]]
    if (is.null(exp_set)) return(NA)
    any(tolower(broad) %in% tolower(exp_set))
  }, per_cluster$top_panel, per_cluster$celltype_broad)

  per_cluster$flag <- dplyr::case_when(
    is.na(per_cluster$agrees) ~ "review_no_expected_mapping",
    per_cluster$agrees        ~ "ok",
    TRUE                      ~ "DISAGREES_consider_manual_override"
  )

  ensure_dir(paths_eye$results_tables)
  out_path <- file.path(paths_eye$results_tables, "cluster_lineage_validation.csv")
  write.csv(per_cluster, out_path, row.names = FALSE)
  log_message("Wrote ", out_path)

  # Surface clusters that need attention in the log so the user sees them
  # without opening the CSV.
  disagreers <- per_cluster$eye_cluster[per_cluster$flag == "DISAGREES_consider_manual_override"]
  if (length(disagreers) > 0L) {
    log_message("Lineage disagreement on clusters: ",
                paste(disagreers, collapse = ", "),
                " — consider cfg$eye_manual_overrides entries.")
  } else {
    log_message("Lineage validation: no cluster-level disagreements.")
  }
  invisible(per_cluster)
}
