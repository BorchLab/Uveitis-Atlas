# R/59_bcell_baseline_selection.R
# BASELINe selection sigma per (subject x phenotype x substate) using shazam.
# Heavy compute; gated by cfg$steps$bcell_baseline_selection.
#
# Output:
#   outputs/tables/eye/bcell/architecture/baseline_selection.csv
#
# Inputs come from R/56's .bcell_load_eye_airr(cfg):
#   - inputs$bm: bcell meta with Subject, knn.leiden.cluster, phenotype,
#       cell_id_unique (eye-only B cells).
#   - inputs$bcr_eye: AIRR rows for productive IGH heavy chains of eye cells,
#       with sequence_alignment, germline_alignment_d_mask, subject_id,
#       clone_id, cell_id_unique.
#
# We enrich bcr_eye with phenotype + substate via join on cell_id_unique,
# then run observedMutations -> calcBaseline -> groupBaseline ->
# summarizeBaseline grouped by subject_id x phenotype x substate. Output
# is renamed so subject_id -> subject for the CSV.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

run_bcell_baseline_selection <- function(cfg) {
  log_message("[figS6a] BASELINe selection compute")
  if (!requireNamespace("shazam", quietly = TRUE)) {
    log_message("[figS6a] shazam not installed; skipping.")
    return(invisible(FALSE))
  }

  arch <- cfg$bcr_lineage$architecture
  set.seed(arch$baseline$seed %||% 42)

  # R/56 helper provides bm + bcr_eye + paths.
  if (!exists(".bcell_load_eye_airr"))
    source(file.path("R", "56_bcell_lineage_trees.R"))
  inputs <- .bcell_load_eye_airr(cfg)
  if (is.null(inputs) || is.null(inputs$bcr_eye) ||
      nrow(inputs$bcr_eye) == 0L) {
    log_message("[figS6a] no AIRR rows from .bcell_load_eye_airr; skipping.")
    return(invisible(FALSE))
  }

  paths   <- inputs$paths
  out_dir <- file.path(paths$results_tables, "architecture")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  airr <- inputs$bcr_eye
  bm   <- inputs$bm

  # Build cell-level lookup: cell_id_unique -> phenotype, substate.
  # Substate is the human-readable label via substate_labels(); mirrors
  # the R/56 .bcell_lineage_meta pattern.
  bm$substate <- substate_labels(cfg, "bcell", bm$knn.leiden.cluster)
  meta_lookup <- bm |>
    dplyr::select(cell_id_unique, phenotype, substate) |>
    dplyr::filter(!is.na(phenotype), !is.na(substate))

  airr <- airr |>
    dplyr::inner_join(meta_lookup, by = "cell_id_unique") |>
    dplyr::filter(!is.na(phenotype), !is.na(substate),
                  phenotype %in% c("Viral", "NIU"))

  if (nrow(airr) == 0L) {
    log_message("[figS6a] no AIRR rows after phenotype/substate join; skipping.")
    return(invisible(FALSE))
  }
  log_message("[figS6a] AIRR rows after enrichment: ", nrow(airr))

  # Pick germline column: prefer d-masked (matches R/56's tree builder),
  # fall back to plain germline_alignment if d_mask is missing or all NA.
  germ_col <- if ("germline_alignment_d_mask" %in% colnames(airr) &&
                  !all(is.na(airr$germline_alignment_d_mask))) {
    "germline_alignment_d_mask"
  } else if ("germline_alignment" %in% colnames(airr) &&
             !all(is.na(airr$germline_alignment))) {
    "germline_alignment"
  } else {
    log_message("[figS6a] no usable germline column; skipping.")
    return(invisible(FALSE))
  }
  log_message("[figS6a] using germline column: ", germ_col)

  if (!"sequence_alignment" %in% colnames(airr)) {
    log_message("[figS6a] sequence_alignment column missing; skipping.")
    return(invisible(FALSE))
  }

  test_stat <- arch$baseline$test_statistic %||% "focused"

  # ---- observedMutations ---------------------------------------------------
  # combine=FALSE keeps per-region mu_freq counts on the AIRR rows (CDR vs
  # FWR stay split). calcBaseline reads the alignment/germline columns
  # directly, so the CDR/FWR split in the final summary is driven by
  # regionDefinition=IMGT_V, not by `combine`.
  obs <- tryCatch(
    shazam::observedMutations(
      airr,
      sequenceColumn   = "sequence_alignment",
      germlineColumn   = germ_col,
      regionDefinition = shazam::IMGT_V,
      frequency        = TRUE,
      combine          = FALSE),
    error = function(e) {
      log_message("[figS6a] observedMutations failed: ", conditionMessage(e))
      NULL
    })
  if (is.null(obs)) return(invisible(FALSE))

  # ---- collapseClones ------------------------------------------------------
  # Canonical Immcantation BASELINe workflow: collapse to one effective
  # sequence per clone before calcBaseline. Without this step each
  # non-collapsed sequence is treated as independent, inflating power
  # within large clones and biasing per-(subject x substate) sigma.
  # Emits new columns `clonal_sequence` and `clonal_germline`.
  #
  # shazam::collapseClones requires all sequences within a clone share the
  # same junction length. Our R/56 clonal definition allows minor junction-
  # length drift, so we sub-cluster clones by junction_length to build
  # length-homogeneous sub-clones before collapsing. The original clone_id
  # is preserved as clone_id_orig for downstream traceability.
  if (all(c("clone_id", "junction_length") %in% colnames(obs))) {
    obs$clone_id_orig <- obs$clone_id
    obs$clone_id      <- paste0(obs$clone_id_orig, "_jl", obs$junction_length)
    n_sub <- length(unique(obs$clone_id)) - length(unique(obs$clone_id_orig))
    if (n_sub > 0L)
      log_message("[figS6a] sub-clustered ", n_sub,
                  " clones by junction_length for collapseClones.")
  }

  obs <- tryCatch(
    shazam::collapseClones(
      obs,
      cloneColumn         = "clone_id",
      sequenceColumn      = "sequence_alignment",
      germlineColumn      = germ_col,
      method              = "thresholdedFreq",
      minimumFrequency    = 0.6,
      includeAmbiguous    = FALSE,
      breakTiesStochastic = FALSE),
    error = function(e) {
      log_message("[figS6a] collapseClones failed: ", conditionMessage(e))
      NULL
    })
  if (is.null(obs)) return(invisible(FALSE))
  log_message("[figS6a] rows after collapseClones (clone-level): ", nrow(obs))

  # ---- calcBaseline --------------------------------------------------------
  bl <- tryCatch(
    shazam::calcBaseline(
      obs,
      sequenceColumn   = "clonal_sequence",
      germlineColumn   = "clonal_germline",
      testStatistic    = test_stat,
      regionDefinition = shazam::IMGT_V,
      targetingModel   = shazam::HKL_S5F),
    error = function(e) {
      log_message("[figS6a] calcBaseline failed: ", conditionMessage(e))
      NULL
    })
  if (is.null(bl)) return(invisible(FALSE))

  # ---- groupBaseline -------------------------------------------------------
  grp <- tryCatch(
    shazam::groupBaseline(
      bl,
      groupBy = c("subject_id", "phenotype", "substate")),
    error = function(e) {
      log_message("[figS6a] groupBaseline failed: ", conditionMessage(e))
      NULL
    })
  if (is.null(grp)) return(invisible(FALSE))

  # ---- summarizeBaseline ---------------------------------------------------
  # shazam returns a "Baseline" S4. Try returnType="df" first; if the
  # version refuses, fall back to slot(grp, "stats").
  sum_df <- tryCatch(
    shazam::summarizeBaseline(grp, returnType = "df"),
    error = function(e) {
      log_message("[figS6a] summarizeBaseline df mode failed: ",
                  conditionMessage(e), "; trying slot(grp, 'stats').")
      NULL
    })
  if (is.null(sum_df)) {
    sum_df <- tryCatch(
      as.data.frame(methods::slot(grp, "stats")),
      error = function(e) {
        log_message("[figS6a] slot(grp,'stats') failed: ",
                    conditionMessage(e))
        NULL
      })
  }
  if (is.null(sum_df) || nrow(sum_df) == 0L) {
    log_message("[figS6a] empty BASELINe summary; skipping.")
    return(invisible(FALSE))
  }

  # Rename subject_id -> subject to match R/58's clone_architecture_metrics.csv
  # naming convention (consistent across architecture outputs).
  out <- sum_df |>
    dplyr::rename_with(~ "subject", dplyr::any_of("subject_id"))

  out_path <- file.path(out_dir, "baseline_selection.csv")
  readr::write_csv(out, out_path)
  log_message("[figS6a] wrote ", nrow(out), " rows to ", out_path)
  invisible(TRUE)
}
