# R/55_tcell_paired_eye_blood_metrics.R
# Compute per-(Subject, Tissue) TCR summary metrics from the FULL atlas so
# the Figure 6 T cell Gini / sharing panel works. The tcell compartment
# Seurat object is eye-only, so paired Eye+Blood comparisons must come from
# the full-atlas IntegratedSeuratObject.
#
# Mirrors R/54_bcell_paired_eye_blood_metrics.R. T cells have no SHM or
# class-switching, so only n_cells, shannon_evenness, and gini are written.
#
# Outputs:
#   outputs/tables/eye/tcell/tcell_paired_eye_blood_metrics.csv
#     Cols: Subject, Tissue_1, Phenotype_2, n_cells, shannon_evenness, gini

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

.gini_one_tcr <- function(x) {
  x <- x[x > 0]; n <- length(x)
  if (n < 2) return(NA_real_)
  x <- sort(x)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
}

.shannon_evenness_tcr <- function(counts) {
  f <- counts / sum(counts)
  H <- -sum(f * log(f))
  if (length(counts) > 1) H / log(length(counts)) else NA_real_
}

run_tcell_paired_eye_blood_metrics <- function(cfg) {
  paths_all <- get_target_paths(cfg, "all")
  obj_path  <- file.path(paths_all$results_objects,
                         "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("tcell paired metrics: full atlas missing at ", obj_path,
                "; skipping.")
    return(invisible(FALSE))
  }
  log_message("=== tcell paired eye<->blood metrics (full atlas) ===")
  obj <- readRDS(obj_path)
  meta <- obj@meta.data

  ctb <- resolve_celltype_broad(meta)
  if (is.null(ctb)) {
    log_message("  no celltype_broad column; aborting.")
    return(invisible(FALSE))
  }
  is_t <- meta[[ctb]] %in% c("T cell", "NK", "T/NK")
  meta <- meta[is_t, , drop = FALSE]
  log_message("  T/NK cells in full atlas: ", nrow(meta))

  clone_col <- intersect(c("CTstrict", "CTaa", "tcr_clone_id"),
                         colnames(meta))[1]
  if (is.na(clone_col)) {
    log_message("  no clone column (CTstrict / CTaa / tcr_clone_id); skipping.")
    return(invisible(FALSE))
  }

  d <- meta |>
    dplyr::filter(!is.na(.data[[clone_col]]),
                  Tissue_1 %in% c("Eye", "Blood"),
                  Phenotype_2 %in% c("NIU", "Viral"))
  healthy_set <- as.character(cfg$etiology_groups$healthy %||% "Healthy")
  if ("Etiology" %in% colnames(d)) {
    n_before <- nrow(d)
    d <- d |> dplyr::filter(!(Etiology %in% healthy_set))
    if (nrow(d) < n_before)
      log_message("  dropped ", n_before - nrow(d), " Healthy-etiology cells")
  }
  if (nrow(d) == 0L) {
    log_message("  no cells pass filter; aborting.")
    return(invisible(FALSE))
  }

  summ <- d |>
    dplyr::group_by(Subject, Tissue_1, Phenotype_2) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      shannon_evenness =
        .shannon_evenness_tcr(as.numeric(table(.data[[clone_col]]))),
      gini = .gini_one_tcr(as.numeric(table(.data[[clone_col]]))),
      .groups = "drop"
    )
  out_dir <- file.path("outputs/tables/eye/tcell")
  ensure_dir(out_dir)
  out_path <- file.path(out_dir, "tcell_paired_eye_blood_metrics.csv")
  utils::write.csv(summ, out_path, row.names = FALSE)
  log_message("  wrote: ", out_path, " (", nrow(summ), " rows, clone col = ",
              clone_col, ")")
  invisible(TRUE)
}
