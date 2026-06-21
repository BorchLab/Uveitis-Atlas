# R/84_viz_eye.R
# Figure 2 (eye sub-atlas) visualizations.
# Uses shared viz_* helpers defined in 82_viz_dispatch.R.

suppressPackageStartupMessages({
  library(Seurat)
})

# run_visualizations_eye renders Figure 2 panels for the eye sub-atlas. It
# calls into the shared viz_* helpers (defined in 82_viz_dispatch.R) with
# target = "eye" so the dge/milo/composition/go panels look at eye-specific
# tables. Cross-tissue blocks (repertoire / BCR lineage) are F1-only and
# live in 83_viz_full.R.
run_visualizations_eye <- function(cfg) {
  paths <- get_target_paths(cfg, "eye")
  log_message("=== Figure 2 (eye sub-atlas) visualizations ===")
  ensure_dir(paths$viz_dir)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Eye IntegratedSeuratObject.rds not found. Skipping F2 viz.")
    return(invisible(TRUE))
  }
  obj <- readRDS(obj_path)

  viz_qc_summary(obj, cfg, paths)
  viz_integration(obj, cfg, paths)
  viz_celltypes(obj, cfg, paths)
  viz_markers(obj, cfg, paths)
  viz_dge(obj, cfg, paths, target = "eye")
  viz_milo(obj, cfg, paths, target = "eye")
  viz_composition(cfg, paths, target = "eye")

  log_message("=== Figure 2 visualizations complete ===")
  invisible(TRUE)
}
