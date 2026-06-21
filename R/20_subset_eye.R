# R/20_subset_eye.R
# Subset the integrated full object to eye-only cells, preserve carry-over
# columns under *_full suffixes, save to outputs/objects/eye/.
suppressPackageStartupMessages({
  library(Seurat)
})

subset_eye <- function(cfg) {
  if (!isTRUE(cfg$eye_focus$enable)) {
    log_message("eye_focus disabled. Skipping eye subset.")
    return(invisible(TRUE))
  }

  src_path <- file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(src_path)) {
    log_message("Full integrated object not found at ", src_path, ". Skipping.")
    return(invisible(FALSE))
  }

  log_message("Loading full integrated object for eye subset...")
  obj <- readRDS(src_path)

  tcol <- cfg$eye_focus$tissue_col   %||% "Tissue_1"
  tval <- cfg$eye_focus$tissue_value %||% "Eye"

  if (!tcol %in% colnames(obj[[]])) {
    stop("subset_eye: column ", tcol, " not in object metadata")
  }

  eye_cells <- colnames(obj)[obj[[tcol, drop = TRUE]] == tval]
  log_message(sprintf("Subsetting to %d %s cells (of %d total)",
                      length(eye_cells), tval, ncol(obj)))
  if (length(eye_cells) < 1000L) {
    stop("subset_eye: only ", length(eye_cells), " eye cells. Aborting.")
  }

  eye_obj <- subset(obj, cells = eye_cells)

  # ---- Lens fiber contamination filter ---------------------------------------
  # Intraocular fluid samples can carry lens fiber cells (crystallin/MIP/BFSP+).
  # They form a contaminating Leiden cluster on re-integration. Drop them here
  # so they never enter HVG selection, fastMNN, or clustering.
  lens_cfg <- cfg$eye_qc$lens_filter %||% list()
  if (isTRUE(lens_cfg$enable)) {
    lens_genes <- intersect(lens_cfg$genes %||% character(0), rownames(eye_obj))
    lens_thr   <- as.numeric(lens_cfg$threshold %||% 3.0)
    if (length(lens_genes) >= 3L) {
      DefaultAssay(eye_obj) <- "RNA"
      eye_obj <- JoinLayers(eye_obj)
      expr <- tryCatch(
        GetAssayData(eye_obj, assay = "RNA", layer = "data"),
        error = function(e) GetAssayData(eye_obj, assay = "RNA", slot = "data")
      )
      lens_score <- Matrix::colSums(expr[lens_genes, , drop = FALSE])
      drop_mask  <- lens_score > lens_thr
      n_drop     <- sum(drop_mask)
      log_message(sprintf(
        "Lens contamination filter: %d/%d genes present, threshold=%.2f, dropping %d cells (%.2f%%)",
        length(lens_genes), length(lens_cfg$genes), lens_thr,
        n_drop, 100 * n_drop / ncol(eye_obj)))
      if (n_drop > 0 && n_drop < ncol(eye_obj)) {
        eye_obj <- subset(eye_obj, cells = colnames(eye_obj)[!drop_mask])
      }
    } else {
      log_message("Lens filter enabled but <3 lens genes present in object. Skipping.")
    }
  }

  # Rename inherited cluster + celltype labels so re-integration / re-annotation
  # can write fresh ones into the canonical column names.
  if (isTRUE(cfg$eye_focus$rename_inherited_labels)) {
    rename_pairs <- c(
      "knn.leiden.cluster"      = "knn.leiden.cluster_full",
      "celltype"                = "celltype_full",
      "celltype_broad"          = "celltype_broad_full",
      "merged.celltype.cluster" = "merged.celltype.cluster_full"
    )
    md <- eye_obj[[]]
    for (old in names(rename_pairs)) {
      new <- rename_pairs[[old]]
      if (old %in% colnames(md)) {
        eye_obj[[new]] <- md[[old]]
        eye_obj[[old]] <- NULL
        log_message("  Renamed ", old, " -> ", new)
      }
    }
  }

  # Drop full-object reductions; they will be replaced by the eye-only
  # re-integration in the next step.
  for (rd in c("fastMNN", "UMAP", "pca")) {
    if (rd %in% Reductions(eye_obj)) {
      eye_obj[[rd]] <- NULL
    }
  }

  out_paths <- get_target_paths(cfg, "eye")
  ensure_dir(out_paths$results_objects)
  out_path <- file.path(out_paths$results_objects, "IntegratedSeuratObject.rds")
  saveRDS(eye_obj, out_path)
  log_message("Saved eye subset to ", out_path)

  invisible(TRUE)
}
