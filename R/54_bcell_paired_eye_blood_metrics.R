# R/54_bcell_paired_eye_blood_metrics.R
# Compute per-(Subject, Tissue) BCR summary metrics from the FULL atlas so
# Figure 6 cover and Gini panels work. The bcell compartment Seurat object
# is eye-only (subset_compartments runs on the eye sub-atlas), so paired
# Eye+Blood comparisons must come from the full-atlas IntegratedSeuratObject.
#
# Outputs:
#   outputs/tables/eye/bcell/bcell_paired_eye_blood_metrics.csv
#     Cols: Subject, Tissue_1, Phenotype_2, n_cells, shannon_evenness,
#           class_switched_frac, mean_shm, gini

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

.gini_one <- function(x) {
  x <- x[x > 0]; n <- length(x)
  if (n < 2) return(NA_real_)
  x <- sort(x)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
}

.shannon_evenness <- function(counts) {
  f <- counts / sum(counts)
  H <- -sum(f * log(f))
  if (length(counts) > 1) H / log(length(counts)) else NA_real_
}

run_bcell_paired_eye_blood_metrics <- function(cfg) {
  paths_all <- get_target_paths(cfg, "all")
  obj_path  <- file.path(paths_all$results_objects,
                         "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("bcell paired metrics: full atlas missing at ", obj_path,
                "; skipping.")
    return(invisible(FALSE))
  }
  log_message("=== bcell paired eye<->blood metrics (full atlas) ===")
  obj <- readRDS(obj_path)
  meta <- obj@meta.data

  # B/Plasma cells only
  ctb <- resolve_celltype_broad(meta)
  if (is.null(ctb)) {
    log_message("  no celltype_broad column; aborting.")
    return(invisible(FALSE))
  }
  is_b <- meta[[ctb]] %in% c("B cell", "Plasma")
  meta <- meta[is_b, , drop = FALSE]
  log_message("  B/plasma cells in full atlas: ", nrow(meta))

  clone_col <- if ("CTstrict" %in% colnames(meta)) "CTstrict" else "bcr_clone_id"
  if (!clone_col %in% colnames(meta)) {
    log_message("  no clone column (CTstrict / bcr_clone_id); skipping.")
    return(invisible(FALSE))
  }
  iso_col <- intersect(c("c_call_heavy","c_call","isotype"), colnames(meta))[1]
  shm_total <- NULL
  shm_cols <- intersect(c("mu_freq_cdr_r_heavy","mu_freq_cdr_s_heavy",
                          "mu_freq_fwr_r_heavy","mu_freq_fwr_s_heavy"),
                        colnames(meta))
  if (length(shm_cols) > 0L)
    shm_total <- rowSums(as.matrix(meta[, shm_cols, drop = FALSE]),
                         na.rm = TRUE)

  d <- meta |>
    dplyr::filter(!is.na(.data[[clone_col]]),
                  Tissue_1 %in% c("Eye","Blood"),
                  Phenotype_2 %in% c("NIU","Viral"))
  # Drop Healthy etiology — Fig 6 thesis is NIU vs Viral; Healthy subjects
  # are PBMC-only controls with no paired eye data.
  healthy_set <- as.character(cfg$etiology_groups$healthy %||% "Healthy")
  if ("Etiology" %in% colnames(d)) {
    n_before <- nrow(d)
    d <- d |> dplyr::filter(!(Etiology %in% healthy_set))
    if (nrow(d) < n_before)
      log_message("  dropped ", n_before - nrow(d), " Healthy-etiology cells")
  }
  if (!is.null(shm_total)) d$SHM_total <- shm_total[is_b][match(rownames(d),
                                                                rownames(meta))]
  if (nrow(d) == 0L) {
    log_message("  no cells pass filter; aborting.")
    return(invisible(FALSE))
  }

  summ <- d |>
    dplyr::group_by(Subject, Tissue_1, Phenotype_2) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      shannon_evenness = .shannon_evenness(as.numeric(table(.data[[clone_col]]))),
      class_switched_frac = if (!is.na(iso_col)) {
        cc <- sub("\\*.*$", "", as.character(.data[[iso_col]]))
        mean(grepl("^IGH[GAE]", cc), na.rm = TRUE)
      } else NA_real_,
      mean_shm = if ("SHM_total" %in% colnames(d)) {
        mean(SHM_total, na.rm = TRUE)
      } else NA_real_,
      gini = .gini_one(as.numeric(table(.data[[clone_col]]))),
      .groups = "drop"
    )
  out_dir <- file.path("outputs/tables/eye/bcell")
  ensure_dir(out_dir)
  out_path <- file.path(out_dir, "bcell_paired_eye_blood_metrics.csv")
  utils::write.csv(summ, out_path, row.names = FALSE)
  log_message("  wrote: ", out_path, " (", nrow(summ), " rows)")
  invisible(TRUE)
}
