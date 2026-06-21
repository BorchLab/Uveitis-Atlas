# R/23_substate_labels.R
# Helpers that propagate compartment-substate keys across the pipeline so that
# downstream modules (PCA per substate, cross-compartment LR, panel labelling)
# read from a single consistent column rather than each rebuilding its own
# label from cluster IDs.
#
# Conventions (locked 2026-05-19):
#   substate_key      Stable machine-readable key per cell, of the form
#                     "<compartment>_<cluster_id>" (e.g. "myeloid_0",
#                     "tcell_3"). Same scheme `substate_joint` uses on the
#                     eye object so cross-compartment merges hit a single key.
#   substate_display  Human-readable label per cell, of the form
#                     "<id>: <curated label>" (e.g. "0: Classical
#                     inflammatory mono"). Falls back to "<id>: cluster_<id>"
#                     when cfg$compartment_substate_labels has no entry.
#   substate_joint    Eye-object column built by `build_joint_substate_labels`:
#                     every cell on the eye object that lives in a compartment
#                     gets the matching `substate_key`; cells outside any
#                     compartment get "other_<eye_cluster>" so cross-cell
#                     analyses see a single coherent grouping.
#
# This file replaces the substate-labelling half of R/64_cellchat.R (the
# CellChat runner itself is retired in favour of LIANA — see R/47).
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Pull the curated human label for one compartment + cluster id out of
# cfg$compartment_substate_labels (config/config.yml:601-628). Returns
# "cluster_<id>" when no curated label exists so plots stay legible.
get_substate_display <- function(cfg, target, id) {
  labels <- cfg$compartment_substate_labels[[target]]
  if (!is.null(labels) && !is.null(labels[[as.character(id)]])) {
    paste0(as.character(id), ": ", labels[[as.character(id)]])
  } else {
    paste0(as.character(id), ": cluster_", as.character(id))
  }
}

# Resolve a substate_key vector for a compartment object that may or may not
# have been stamped by apply_substate_keys (older objects don't carry the
# column yet). Reconstructs from knn.leiden.cluster as fallback. Used by
# R/47 and R/48 so the cross-compartment modules don't need both objects to
# be re-saved before they can run.
get_substate_key_vector <- function(obj, target) {
  if ("substate_key" %in% colnames(obj[[]])) return(as.character(obj$substate_key))
  cluster_col <- if ("knn.leiden.cluster" %in% colnames(obj[[]]))
                   "knn.leiden.cluster" else "seurat_clusters"
  paste0(target, "_", as.character(obj[[cluster_col, drop = TRUE]]))
}

# Stamp `substate_key` and `substate_display` columns onto a compartment
# Seurat object. Reads `knn.leiden.cluster` (the compartment's Leiden labels)
# and writes the two columns in place. Called from R/22 right before saving
# the compartment object so every downstream reader sees the same keys.
apply_substate_keys <- function(obj, target, cfg) {
  cluster_col <- if ("knn.leiden.cluster" %in% colnames(obj[[]]))
                   "knn.leiden.cluster" else "seurat_clusters"
  if (!cluster_col %in% colnames(obj[[]])) {
    log_message("  apply_substate_keys: no cluster column on ", target,
                " object; skipping.")
    return(obj)
  }
  ids <- as.character(obj[[cluster_col, drop = TRUE]])
  obj$substate_key     <- paste0(target, "_", ids)
  obj$substate_display <- vapply(ids,
                                 function(id) get_substate_display(cfg, target, id),
                                 character(1), USE.NAMES = FALSE)
  obj
}

# Map per-compartment Leiden cluster IDs onto the eye sub-atlas as a unified
# `substate_joint` column. Cells in the eye object that did not survive any
# compartment subset (lens-contaminated, NK outside the three F3-F5 lineages)
# get a generic "other_<eye_cluster>" label so downstream code still has a
# coherent grouping rather than NA groups.
#
# Relocated from R/64_cellchat.R:32-81. Reads compartment objects from disk
# via get_target_paths(); writes EyeJointSubstate.rds.
build_joint_substate_labels <- function(cfg) {
  eye_paths <- get_target_paths(cfg, "eye")
  eye_path  <- file.path(eye_paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(eye_path)) {
    log_message("Joint substates: eye object missing at ", eye_path,
                ". Run eye_subset / eye_reintegrate first.")
    return(invisible(FALSE))
  }
  eye <- readRDS(eye_path)
  eye$substate_joint <- NA_character_

  cluster_col <- "knn.leiden.cluster"
  for (cmp in c("myeloid", "bcell", "tcell")) {
    cmp_path <- file.path(get_target_paths(cfg, cmp)$results_objects,
                          "IntegratedSeuratObject.rds")
    if (!file.exists(cmp_path)) {
      log_message("  Joint substates: missing ", cmp, " object; skipping ", cmp)
      next
    }
    cmp_obj <- readRDS(cmp_path)
    if ("substate_key" %in% colnames(cmp_obj[[]])) {
      labels <- as.character(cmp_obj$substate_key)
    } else if (cluster_col %in% colnames(cmp_obj[[]])) {
      labels <- paste0(cmp, "_", as.character(cmp_obj[[cluster_col, drop = TRUE]]))
    } else {
      log_message("  Joint substates: ", cmp, " object has neither substate_key ",
                  "nor ", cluster_col, "; skipping.")
      next
    }
    idx <- match(colnames(cmp_obj), colnames(eye))
    keep <- !is.na(idx)
    eye$substate_joint[idx[keep]] <- labels[keep]
    log_message("  Joint substates: assigned ", sum(keep), " ", cmp, " cells")
  }

  eye_cluster_col <- if ("knn.leiden.cluster" %in% colnames(eye[[]]))
                       "knn.leiden.cluster" else "seurat_clusters"
  missing <- is.na(eye$substate_joint)
  if (any(missing) && eye_cluster_col %in% colnames(eye[[]])) {
    eye$substate_joint[missing] <- paste0(
      "other_", as.character(eye[[eye_cluster_col, drop = TRUE]][missing]))
  }
  log_message("Joint substates: ", length(unique(eye$substate_joint)),
              " group labels (", sum(grepl("^other_", eye$substate_joint)),
              " unassigned cells parked under 'other_*').")

  out_path <- file.path(eye_paths$results_objects, "EyeJointSubstate.rds")
  saveRDS(eye, out_path)
  log_message("Joint substates: saved to ", out_path)
  invisible(TRUE)
}
