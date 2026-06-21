# R/16_compartment_lineage_gate.R
#
# Canonical-marker lineage gate applied at the compartment-subset boundary
# (R/22_subset_compartments.R). For each candidate cell of a compartment, score
# 15 marker panels by within-panel detection rate; keep the cell only if its
# top-scoring panel is in the compartment's expected_panels list and it does
# not co-express any non-immune panel at the same threshold. Optional strict
# mode also drops cells positive for an off-target immune panel.
#
# This catches the cross-lineage contamination that the consensus broad-label
# annotator (PBMC-trained) misses — e.g., mregDCs routed into bcell, CD3+ T
# cells routed into myeloid, CD45-low debris of any origin.
#
# Panels are intentionally distinct from R/17_lineage_validation.R: tighter
# Plasma panel (no XBP1/PRDM1 — both leak into stressed T cells), tighter pDC
# panel (no TCF4/RUNX2 — both leak into B cells), and a CD16 fix
# (FCGR3A moved out of the NK panel into Mono_Mac since it dominates in
# CD16+ monocytes and falsely flags monocytes as NK contamination).
#
# Public functions:
#   apply_compartment_lineage_gate(eye_obj, cells, cmp, cfg) ->
#       list(cells = filtered_cells, audit = per-cell pass/fail tibble)
#
# Audit rows are accumulated by 22_subset_compartments.R and written to
#   outputs/tables/qc/compartment_gate_report.csv

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
})

# ---- Panel definitions ------------------------------------------------------
.gate_panels_immune <- function() list(
  T_cell      = c("CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2","BCL11B"),
  NK          = c("NKG7","GNLY","KLRD1","KLRF1","NCAM1","KLRC1"),
  B_cell      = c("CD79A","CD79B","MS4A1","CD19","BANK1","FCRLA","PAX5"),
  Plasma      = c("SDC1","DERL3","MZB1","JCHAIN","IGHG1","IGHA1","IGKC","IGLC2"),
  Mono_Mac    = c("CD14","CD68","LYZ","C1QA","C1QB","C1QC","CSF1R","MARCO",
                  "VCAN","S100A12","FCN1","FCGR3A"),
  cDC         = c("FCER1A","CD1C","CD1A","CLEC10A","CLEC9A","IRF8","IRF4","CD207"),
  mregDC      = c("LAMP3","CCR7","FSCN1","IDO1","CD274","CCL19","CCL22"),
  pDC         = c("CLEC4C","LILRA4","IL3RA","GZMB","ITM2C"),
  Granulocyte = c("S100A8","S100A9","FCGR3B","CSF3R","MPO","ELANE","CXCR2")
)

.gate_panels_nonimmune <- function() list(
  Stromal_fibro  = c("COL1A1","COL1A2","COL3A1","DCN","LUM","AEBP1","PDGFRA","FAP"),
  Epithelial     = c("EPCAM","KRT5","KRT14","KRT8","KRT18","KRT19","KRT3","KRT12","KRT15"),
  Endothelial   = c("PECAM1","CDH5","VWF","KDR","CLDN5"),
  Lens           = c("CRYAA","CRYAB","BFSP1","BFSP2","MIP","LIM2","CRYBB1","CRYBB2","CRYGD","LGSN"),
  Ocular_pigment = c("TYR","TYRP1","DCT","PMEL","MITF","BEST1","SILV"),
  Neural         = c("RHO","RCVRN","NEFL","NEFM","NRL","MAP1A","SCN7A","CDH2","NRN1")
)

# ---- Cell scoring -----------------------------------------------------------
# Per-cell within-panel detection rate. Returns a (cells x panels) numeric
# matrix of fractions; NA entries indicate the panel had no genes present in
# the assay.
.gate_score_cells <- function(data_mat, cells, panels) {
  data_sub <- data_mat[, cells, drop = FALSE]
  out <- matrix(NA_real_, nrow = length(cells), ncol = length(panels),
                dimnames = list(cells, names(panels)))
  for (pn in names(panels)) {
    g <- intersect(panels[[pn]], rownames(data_sub))
    if (length(g) == 0L) next
    if (length(g) == 1L) {
      out[, pn] <- as.numeric(data_sub[g, ] > 0)
    } else {
      m <- data_sub[g, , drop = FALSE]
      out[, pn] <- Matrix::colSums(m > 0) / length(g)
    }
  }
  out
}

# ---- Public entry point -----------------------------------------------------
# Applies the gate to the candidate cells of one compartment. Reads the gate
# config under cfg$compartment_lineage_gate; if absent or disabled the
# function is a no-op (returns cells unchanged with an empty audit).
apply_compartment_lineage_gate <- function(eye_obj, cells, cmp, cfg) {
  gate_cfg <- cfg$compartment_lineage_gate
  if (is.null(gate_cfg) || !isTRUE(gate_cfg$enable)) {
    return(list(cells = cells, audit = NULL))
  }
  if (length(cells) == 0L) {
    return(list(cells = cells, audit = NULL))
  }

  threshold   <- as.numeric(gate_cfg$panel_threshold %||% 0.25)
  # strict_mode can be a single logical (applied to all compartments) OR a
  # named list / character vector of compartments where strict mode should
  # be on. The latter is useful when one compartment (e.g. myeloid) has
  # doublet contamination that off-target-immune filtering catches, while
  # other compartments (bcell/tcell) lose too many legitimate cells under
  # the same strictness.
  sm <- gate_cfg$strict_mode
  strict_mode <- if (is.null(sm)) FALSE
                 else if (is.logical(sm) && length(sm) == 1L) isTRUE(sm)
                 else if (is.list(sm)) isTRUE(sm[[cmp]])
                 else cmp %in% as.character(sm)
  expected    <- gate_cfg$expected_panels[[cmp]]
  if (is.null(expected)) {
    log_message(sprintf(
      "Gate %s: no expected_panels entry; skipping (cells unchanged)", cmp))
    return(list(cells = cells, audit = NULL))
  }
  expected <- as.character(expected)
  # Off-target immune panels to ignore for THIS compartment when strict_mode is
  # on. Use case: pDCs co-express Plasma markers (JCHAIN/MZB1/IGKC) as part of
  # their biology, so the Plasma panel should not be treated as off-target for
  # the myeloid compartment.
  offtarget_exclude <- as.character(gate_cfg$offtarget_immune_exclude[[cmp]] %||% character())

  panels_immune    <- .gate_panels_immune()
  panels_nonimmune <- .gate_panels_nonimmune()
  panels_all       <- c(panels_immune, panels_nonimmune)

  unknown <- setdiff(expected, names(panels_immune))
  if (length(unknown) > 0L) {
    stop(sprintf("Gate %s: expected_panels contains unknown panel(s): %s",
                 cmp, paste(unknown, collapse = ",")), call. = FALSE)
  }

  # Use RNA log-normalized expression. Compartment object is not created yet;
  # we read from the eye object directly to avoid double-loading.
  DefaultAssay(eye_obj) <- "RNA"
  data_mat <- tryCatch(GetAssayData(eye_obj, assay = "RNA", layer = "data"),
                       error = function(e) GetAssayData(eye_obj, assay = "RNA", slot = "data"))
  if (is.null(data_mat) || nrow(data_mat) == 0L) {
    log_message(sprintf("Gate %s: RNA `data` layer empty; normalizing first.", cmp))
    eye_obj <- NormalizeData(eye_obj, verbose = FALSE)
    data_mat <- GetAssayData(eye_obj, assay = "RNA", layer = "data")
  }

  log_message(sprintf("Gate %s: scoring %d cells across %d panels (threshold=%.2f, strict=%s)",
                      cmp, length(cells), length(panels_all),
                      threshold, ifelse(strict_mode, "TRUE", "FALSE")))

  scores <- .gate_score_cells(data_mat, cells, panels_all)
  is_pos <- scores >= threshold
  is_pos[is.na(is_pos)] <- FALSE

  # Top panel per cell
  top_idx <- max.col(scores, ties.method = "first")
  top_panel <- colnames(scores)[top_idx]
  top_score <- scores[cbind(seq_along(top_idx), top_idx)]
  top_panel[is.na(top_score)] <- NA_character_

  nonimmune_panels <- names(panels_nonimmune)
  offtarget_immune <- setdiff(names(panels_immune), c(expected, offtarget_exclude))

  pos_nonimmune <- rowSums(is_pos[, nonimmune_panels, drop = FALSE]) > 0
  pos_offtarget <- if (length(offtarget_immune) > 0L)
                     rowSums(is_pos[, offtarget_immune, drop = FALSE]) > 0
                   else
                     rep(FALSE, length(cells))
  pos_expected  <- rowSums(is_pos[, expected, drop = FALSE]) > 0

  top_in_expected <- !is.na(top_panel) & top_panel %in% expected

  # PASS rules:
  #   top_panel must be in expected; AND
  #   no non-immune panel above threshold; AND
  #   in strict mode, no off-target immune panel above threshold.
  fail_reason <- rep(NA_character_, length(cells))
  fail_reason[!top_in_expected] <- paste0("top_panel_not_expected(", top_panel[!top_in_expected], ")")
  fail_reason[is.na(fail_reason) & pos_nonimmune] <- "nonimmune_panel_positive"
  if (strict_mode) {
    fail_reason[is.na(fail_reason) & pos_offtarget] <- "offtarget_immune_panel_positive"
  }
  if (any(!pos_expected)) {
    fail_reason[is.na(fail_reason) & !pos_expected] <- "no_expected_panel_positive"
  }
  pass <- is.na(fail_reason)

  audit <- data.frame(
    cell_id     = cells,
    compartment = cmp,
    top_panel   = top_panel,
    top_score   = round(top_score, 3),
    pos_expected = pos_expected,
    pos_nonimmune = pos_nonimmune,
    pos_offtarget_immune = pos_offtarget,
    pass        = pass,
    fail_reason = fail_reason,
    stringsAsFactors = FALSE
  )

  log_message(sprintf(
    "Gate %s: kept %d / %d cells (%.1f%%); top fail reasons: %s",
    cmp, sum(pass), length(cells), 100 * mean(pass),
    paste(head(sort(table(fail_reason), decreasing = TRUE), 3) |>
            (\(t) sprintf("%s=%d", names(t), as.integer(t)))(),
          collapse = ", ")))

  list(cells = cells[pass], audit = audit)
}
