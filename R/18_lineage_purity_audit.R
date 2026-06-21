# R/18_lineage_purity_audit.R
#
# Tier 0 audit (2026-05-15). Quantify per-cluster canonical-lineage purity
# across the eye, bcell, myeloid, and tcell objects. Read-only.
#
# Method:
#   * Score 15 curated marker panels per cell as "panel-positive" if
#     detection_rate >= panel_threshold (default 0.25 of panel genes with
#     log-norm > 0). PTPRC handled separately as a single-gene CD45 gate.
#   * Per-cell purity = is_cd45_pos AND (any expected immune lineage positive)
#                       AND (no non-immune panel positive)
#                       AND (no off-target immune lineage positive).
#   * Per-cluster purity = fraction of pure cells.
#
# Outputs:
#   outputs/tables/qc/cluster_lineage_purity.csv
#   outputs/viz/qc/cluster_lineage_purity_heatmap.pdf
#
# Entry point: run_lineage_purity_audit(cfg). Gated by
# cfg$steps$lineage_purity_audit in run_pipeline.R (Tier 0 diagnostics).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

OBJECT_PATHS <- list(
  eye     = "outputs/objects/eye/IntegratedSeuratObject.rds",
  bcell   = "outputs/objects/eye/bcell/IntegratedSeuratObject.rds",
  myeloid = "outputs/objects/eye/myeloid/IntegratedSeuratObject.rds",
  tcell   = "outputs/objects/eye/tcell/IntegratedSeuratObject.rds"
)

OUT_TABLE <- "outputs/tables/qc/cluster_lineage_purity.csv"
OUT_VIZ   <- "outputs/viz/01_qc/cluster_lineage_purity_heatmap.pdf"

PANEL_THRESHOLD <- 0.25   # min detection-rate within panel to call cell positive

# ---- Marker panels ----------------------------------------------------------
# Immune lineages
PANELS_IMMUNE <- list(
  T_cell      = c("CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2","BCL11B"),
  NK          = c("NCAM1","NKG7","GNLY","KLRD1","KLRF1","FCGR3A","KLRC1"),
  B_cell      = c("CD79A","CD79B","MS4A1","CD19","BANK1","FCRLA","PAX5"),
  Plasma      = c("JCHAIN","MZB1","XBP1","PRDM1","SDC1","DERL3","FKBP11","TENT5C"),
  Mono_Mac    = c("CD14","CD68","LYZ","C1QA","C1QB","C1QC","CSF1R","MARCO","VCAN","S100A12","FCN1"),
  cDC         = c("FCER1A","CD1C","CD1A","CLEC10A","CLEC9A","IRF8","IRF4","CD207"),
  mregDC      = c("LAMP3","CCR7","FSCN1","IDO1","CD274","CCL19","CCL22"),
  pDC         = c("CLEC4C","LILRA4","IL3RA","TCF4","RUNX2"),
  Granulocyte = c("S100A8","S100A9","FCGR3B","CSF3R","MPO","ELANE","CXCR2")
)

# Non-immune (any positive => contamination for an immune compartment)
PANELS_NONIMMUNE <- list(
  Stromal_fibro  = c("COL1A1","COL1A2","COL3A1","DCN","LUM","AEBP1","PDGFRA","FAP"),
  Epithelial     = c("EPCAM","KRT5","KRT14","KRT8","KRT18","KRT19","KRT3","KRT12","KRT15"),
  Endothelial    = c("PECAM1","CDH5","VWF","KDR","CLDN5"),
  Lens           = c("CRYAA","CRYAB","BFSP1","BFSP2","MIP","LIM2","CRYBB1","CRYBB2","CRYGD","LGSN"),
  Ocular_pigment = c("TYR","TYRP1","DCT","PMEL","MITF","BEST1","SILV"),
  Neural         = c("RHO","RCVRN","NEFL","NEFM","NRL","MAP1A","SCN7A","CDH2","NRN1")
)

PANELS_ALL <- c(PANELS_IMMUNE, PANELS_NONIMMUNE)

# Compartment -> expected immune panel names
EXPECTED <- list(
  eye     = names(PANELS_IMMUNE),
  bcell   = c("B_cell","Plasma"),
  myeloid = c("Mono_Mac","cDC","mregDC","pDC","Granulocyte"),
  tcell   = c("T_cell","NK")
)

# log_message / ensure_dir / %||% come from R/01_setup_utils.R.

# Returns named numeric vector of per-cell detection rate (fraction of panel
# genes with log-norm > 0). Genes not in the assay are skipped.
panel_detection_rate <- function(data_mat, panel) {
  present <- intersect(panel, rownames(data_mat))
  if (length(present) == 0) {
    return(rep(NA_real_, ncol(data_mat)))
  }
  if (length(present) == 1) {
    # detection rate is just 0/1 per cell for a single gene
    return(as.numeric(data_mat[present, ] > 0))
  }
  m <- data_mat[present, , drop = FALSE]
  Matrix::colSums(m > 0) / length(present)
}

# Returns vector of which panels were present in the assay (i.e. detection
# rate not NA). Used for the panel_coverage report.
audit_one_object <- function(obj_path, compartment) {
  log_message("=== ", compartment, " ===")
  log_message("loading ", obj_path)
  obj <- readRDS(obj_path)
  DefaultAssay(obj) <- "RNA"

  # Normalize only if data layer is empty (assume Seurat v5)
  needs_norm <- tryCatch({
    d <- GetAssayData(obj, assay = "RNA", layer = "data")
    nrow(d) == 0L || all(d@x == 0)
  }, error = function(e) TRUE)
  if (isTRUE(needs_norm)) {
    log_message("normalizing RNA assay")
    obj <- NormalizeData(obj, verbose = FALSE)
  }

  data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
  meta <- obj[[]]

  cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                   "knn.leiden.cluster"
                 else if ("seurat_clusters" %in% colnames(meta))
                   "seurat_clusters"
                 else stop("No cluster column found for ", compartment)
  clusters <- as.character(meta[[cluster_col]])

  # CD45+ per cell
  if (!"PTPRC" %in% rownames(data_mat)) {
    log_message("WARNING: PTPRC absent from RNA assay — pct_cd45_pos will be NA")
    is_cd45 <- rep(NA, ncol(obj))
  } else {
    is_cd45 <- as.numeric(data_mat["PTPRC", ]) > 0
  }

  # Per-cell, per-panel detection rate and positivity
  log_message("scoring ", length(PANELS_ALL), " panels for ", ncol(obj), " cells")
  panel_pos <- vapply(PANELS_ALL, function(panel) {
    rate <- panel_detection_rate(data_mat, panel)
    if (all(is.na(rate))) return(rep(NA, ncol(obj)))
    rate >= PANEL_THRESHOLD
  }, FUN.VALUE = logical(ncol(obj)))
  storage.mode(panel_pos) <- "logical"
  # panel_pos: cells x panels
  colnames(panel_pos) <- names(PANELS_ALL)

  expected_panels <- EXPECTED[[compartment]]
  immune_panels   <- names(PANELS_IMMUNE)
  nonimmune_panels <- names(PANELS_NONIMMUNE)

  any_expected_pos <- rowSums(panel_pos[, expected_panels, drop = FALSE], na.rm = TRUE) > 0
  any_nonimmune_pos <- rowSums(panel_pos[, nonimmune_panels, drop = FALSE], na.rm = TRUE) > 0
  offtarget_immune <- setdiff(immune_panels, expected_panels)
  any_offtarget_pos <- if (length(offtarget_immune) > 0)
                         rowSums(panel_pos[, offtarget_immune, drop = FALSE], na.rm = TRUE) > 0
                       else
                         rep(FALSE, nrow(panel_pos))

  cd45_for_purity <- if (all(is.na(is_cd45))) rep(TRUE, length(is_cd45)) else is_cd45
  is_pure <- cd45_for_purity & any_expected_pos & !any_nonimmune_pos & !any_offtarget_pos

  # Per-cluster aggregation
  per_cell <- data.frame(
    cluster = clusters,
    is_cd45 = is_cd45,
    any_expected_pos = any_expected_pos,
    any_nonimmune_pos = any_nonimmune_pos,
    any_offtarget_pos = any_offtarget_pos,
    is_pure = is_pure,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  per_cell <- cbind(per_cell, as.data.frame(panel_pos))

  agg <- per_cell %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      pct_cd45_pos          = round(100 * mean(is_cd45, na.rm = TRUE), 1),
      pct_pure              = round(100 * mean(is_pure), 1),
      pct_contam_nonimmune  = round(100 * mean(any_nonimmune_pos), 1),
      pct_lineage_mismatch  = round(100 * mean(any_offtarget_pos & !any_nonimmune_pos), 1),
      dplyr::across(dplyr::all_of(names(PANELS_ALL)),
                    \(x) round(100 * mean(x, na.rm = TRUE), 1),
                    .names = "pct_{.col}_pos"),
      .groups = "drop"
    )

  # Top panel per cluster
  panel_pct_cols <- paste0("pct_", names(PANELS_ALL), "_pos")
  panel_mat <- as.matrix(agg[, panel_pct_cols])
  colnames(panel_mat) <- names(PANELS_ALL)
  top_idx <- apply(panel_mat, 1, function(r) {
    if (all(is.na(r))) return(NA_integer_)
    which.max(r)
  })
  agg$top_panel <- ifelse(is.na(top_idx), NA_character_, names(PANELS_ALL)[top_idx])
  agg$top_panel_pct <- vapply(seq_len(nrow(agg)), function(i) {
    if (is.na(top_idx[i])) NA_real_ else panel_mat[i, top_idx[i]]
  }, numeric(1))

  agg$expected_lineage <- paste(expected_panels, collapse = "|")
  agg$compartment <- compartment
  agg$flag <- dplyr::case_when(
    agg$pct_pure < 25 ~ "severe_contamination",
    agg$pct_pure < 50 ~ "low_purity",
    TRUE              ~ "ok"
  )

  # Reorder columns
  lead_cols <- c("compartment","cluster","n_cells","expected_lineage",
                 "pct_pure","pct_cd45_pos","pct_contam_nonimmune",
                 "pct_lineage_mismatch","top_panel","top_panel_pct","flag")
  agg <- agg[, c(lead_cols, panel_pct_cols)]

  agg <- agg[order(as.numeric(agg$cluster), agg$cluster), ]

  rm(obj, data_mat, panel_pos, per_cell); gc(verbose = FALSE)
  agg
}

# ---- Run --------------------------------------------------------------------

run_lineage_purity_audit <- function(cfg = NULL) {
  ensure_dir("outputs/tables/qc")
  ensure_dir("outputs/viz/01_qc")

  all_rows <- list()
  for (comp in names(OBJECT_PATHS)) {
    p <- OBJECT_PATHS[[comp]]
    if (!file.exists(p)) {
      log_message("SKIP ", comp, ": object not found at ", p)
      next
    }
    all_rows[[comp]] <- audit_one_object(p, comp)
  }
  if (length(all_rows) == 0) stop("No objects loaded; nothing to audit.")
  out <- do.call(rbind, all_rows)
  write.csv(out, OUT_TABLE, row.names = FALSE)
  log_message("wrote ", OUT_TABLE)

  # Console summary
  cat("\n========= FLAGGED CLUSTERS =========\n")
  flagged <- subset(out, flag != "ok")
  if (nrow(flagged) == 0) {
    cat("None — every cluster >=50% pure. Unexpected; sanity-check thresholds.\n")
  } else {
    print(flagged[, c("compartment","cluster","n_cells","pct_pure",
                      "pct_cd45_pos","pct_contam_nonimmune",
                      "pct_lineage_mismatch","top_panel","flag")],
          row.names = FALSE)
  }
  cat("\n========= PER-COMPARTMENT TOTALS =========\n")
  tot <- out %>%
    dplyr::group_by(compartment) %>%
    dplyr::summarise(
      n_clusters       = dplyr::n(),
      n_flagged        = sum(flag != "ok"),
      n_cells_total    = sum(n_cells),
      n_cells_in_flagged = sum(n_cells[flag != "ok"]),
      pct_cells_flagged = round(100 * sum(n_cells[flag != "ok"]) / sum(n_cells), 1),
      .groups = "drop"
    )
  print(tot, row.names = FALSE)

  # Heatmap viz
  panel_pct_cols <- paste0("pct_", names(PANELS_ALL), "_pos")
  long <- out %>%
    dplyr::select(compartment, cluster, pct_pure, dplyr::all_of(panel_pct_cols)) %>%
    tidyr::pivot_longer(dplyr::all_of(panel_pct_cols),
                        names_to = "panel", values_to = "pct_pos") %>%
    dplyr::mutate(
      panel = sub("^pct_", "", sub("_pos$", "", panel)),
      panel_type = ifelse(panel %in% names(PANELS_IMMUNE), "immune", "non-immune"),
      panel = factor(panel, levels = c(names(PANELS_IMMUNE), names(PANELS_NONIMMUNE)))
    )

  # one facet per compartment
  p <- ggplot(long, aes(x = factor(cluster,
                                    levels = sort(unique(as.numeric(cluster)))),
                         y = panel)) +
    geom_tile(aes(fill = pct_pos), color = "grey90") +
    geom_text(aes(label = ifelse(pct_pos >= 25, sprintf("%.0f", pct_pos), "")),
              size = 2.4, color = "black") +
    scale_fill_gradient2(low = "white", mid = "#fee08b", high = "#d73027",
                         midpoint = 50, limits = c(0, 100),
                         name = "% cells\npanel-positive") +
    facet_grid(panel_type ~ compartment, scales = "free", space = "free") +
    labs(x = "Cluster (knn.leiden.cluster)", y = NULL,
         title = "Tier 0 lineage purity audit",
         subtitle = sprintf("Cell-level panel positivity (detection rate >= %.0f%%). Each tile = %% cells positive in that cluster.",
                            100 * PANEL_THRESHOLD)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 0),
          panel.grid.major = element_blank(),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey95", color = NA))

  ggsave(OUT_VIZ, p, width = 11, height = 8)
  log_message("wrote ", OUT_VIZ)

  invisible(out)
}
