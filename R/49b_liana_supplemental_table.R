# R/49b_liana_supplemental_table.R
# Combine the per-direction LIANA combined CSVs (NIU vs Viral, with disease
# bias columns) into a single tidy supplemental table for the manuscript.
#
# Reads (from outputs/tables/cross_compartment/):
#   liana_myeloid_to_tcell_combined.csv   (written by R/47)
#   liana_<dir>_combined.csv              (written by R/49 for tcell_to_bcell,
#                                          bcell_to_tcell, myeloid_to_bcell,
#                                          bcell_to_myeloid)
# Writes:
#   liana_supplemental_all_directions.csv
#
# Schema reconciliation:
#   All inputs share source, target, ligand_complex, receptor_complex,
#   ligand_family, per-condition rank/score columns, consensus_rank, and
#   disease_bias_logfc. Only the myeloid->tcell combined adds
#   disease_bias_rank; that column ends up NA for the B cell directions after
#   bind_rows, which is the expected (and explicit) behavior.
#
#   R/49 already stamps a `direction` column on its outputs; R/47 does not.
#   We overwrite with the canonical name from the file map so the column is
#   uniform across both sources.
#
# Missing inputs are logged and skipped rather than fatal, so partial reruns
# (e.g. only Phase 6 enabled) still produce a sensible combined table.
suppressPackageStartupMessages({
  library(dplyr)
})

run_liana_supplemental_table <- function(cfg) {
  cc_paths <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/",
    viz    = "outputs/viz/cross_compartment/")
  tbl_dir <- cc_paths$tables
  ensure_dir(tbl_dir)

  file_map <- c(
    myeloid_to_tcell = "liana_myeloid_to_tcell_combined.csv",
    bcell_to_tcell   = "liana_bcell_to_tcell_combined.csv",
    tcell_to_bcell   = "liana_tcell_to_bcell_combined.csv",
    bcell_to_myeloid = "liana_bcell_to_myeloid_combined.csv",
    myeloid_to_bcell = "liana_myeloid_to_bcell_combined.csv"
  )

  log_message("=== liana_supplemental_table ===")
  parts <- list()
  for (dir_name in names(file_map)) {
    path <- file.path(tbl_dir, file_map[[dir_name]])
    if (!file.exists(path)) {
      log_message("  missing: ", path, " (skipping ", dir_name, ")")
      next
    }
    df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(df) == 0L) {
      log_message("  empty: ", path, " (skipping ", dir_name, ")")
      next
    }
    df$direction <- dir_name
    parts[[dir_name]] <- df
    log_message("  ", dir_name, ": ", nrow(df), " LR pairs")
  }

  if (length(parts) == 0L) {
    log_message("liana_supplemental_table: no input CSVs found; nothing to write.")
    return(invisible(FALSE))
  }

  supp <- dplyr::bind_rows(parts) |>
    dplyr::relocate("direction", .before = 1) |>
    dplyr::arrange(.data$direction, .data$consensus_rank)

  out_path <- file.path(tbl_dir, "liana_supplemental_all_directions.csv")
  utils::write.csv(supp, out_path, row.names = FALSE)
  log_message("  wrote ", nrow(supp), " rows across ", length(parts),
              " directions to ", out_path)
  invisible(TRUE)
}
