# R/22_subset_compartments.R
# Subset the eye sub-atlas into myeloid / bcell / tcell compartment objects,
# re-integrate each, re-tune clustering. Mirrors R/21_integrate_eye.R but
# loops across the three compartments defined in cfg$compartments. A
# UCell-based per-cell lens filter (R/15_lens_qc.R) is applied once to the
# eye object before the split — this replaces the previous cluster-11
# mregDC hygiene rescue, which only operated on one cluster.
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(batchelor)
  library(igraph)
  library(uwot)
  library(Matrix)
  library(dplyr)
})

# Carry forward the eye-level metadata under *_eye suffixes so downstream
# panels can show inheritance. The compartment object then re-derives its own
# seurat_clusters / knn.leiden.cluster.
.carry_forward_eye_metadata <- function(obj) {
  if ("knn.leiden.cluster" %in% colnames(obj@meta.data)) {
    obj$knn.leiden.cluster_eye <- obj$knn.leiden.cluster
  }
  for (col in c("celltype", "celltype_broad", "seurat_clusters")) {
    if (col %in% colnames(obj@meta.data)) {
      obj[[paste0(col, "_eye")]] <- obj[[col]]
    }
  }
  if ("escape.UCell" %in% Seurat::Assays(obj)) {
    obj[["escape.UCell_eye"]] <- obj[["escape.UCell"]]
  }
  obj
}

# T2.5: drop later-visit cells per subject when cfg$compartments$timepoint_policy
# == "earliest". For each Subject with multiple orig.ident values, keep
# only the lexically smallest (typically encodes visit order). Writes an
# audit CSV listing which Subject/orig.ident were kept vs dropped.
.apply_compartment_timepoint_policy <- function(eye_obj, cfg) {
  policy <- cfg$compartments$timepoint_policy %||% "all"
  policy <- match.arg(as.character(policy), choices = c("all", "earliest"))
  if (policy == "all") return(eye_obj)
  meta <- eye_obj@meta.data
  if (!all(c("Subject", "orig.ident") %in% colnames(meta))) {
    log_message("Timepoint policy: Subject or orig.ident missing; leaving object as-is.")
    return(eye_obj)
  }

  per_subj <- split(as.character(meta$orig.ident), as.character(meta$Subject))
  keep_ident_by_subj <- vapply(per_subj, function(idents) {
    ux <- sort(unique(idents))
    ux[1]
  }, character(1))

  keep_mask <- meta$orig.ident == keep_ident_by_subj[as.character(meta$Subject)]
  keep_mask[is.na(keep_mask)] <- TRUE
  n_drop <- sum(!keep_mask)
  log_message(sprintf(
    "Timepoint policy (%s): keeping %d / %d cells (dropping %d later-visit cells)",
    policy, sum(keep_mask), length(keep_mask), n_drop))

  # Audit trail: subjects with >1 visit + which orig.ident was retained.
  multi <- per_subj[vapply(per_subj, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(multi) > 0L) {
    audit <- do.call(rbind, lapply(names(multi), function(s) {
      data.frame(Subject = s,
                 orig.ident = sort(unique(multi[[s]])),
                 kept = sort(unique(multi[[s]])) == keep_ident_by_subj[[s]],
                 policy = policy,
                 stringsAsFactors = FALSE)
    }))
    paths_eye <- get_target_paths(cfg, "eye")
    ensure_dir(paths_eye$results_tables)
    out_path <- file.path(paths_eye$results_tables, "timepoint_policy_audit.csv")
    write.csv(audit, out_path, row.names = FALSE)
    log_message("Wrote ", out_path)
  }

  if (n_drop == 0L) return(eye_obj)
  if (n_drop == ncol(eye_obj))
    stop("Timepoint policy dropped every cell — bug in subject/orig.ident mapping.",
         call. = FALSE)
  subset(eye_obj, cells = colnames(eye_obj)[keep_mask])
}

# Write a per-cell mapping of eye cells to their compartment label. Cells
# unassigned to any compartment appear with compartment == "none"; cells
# matching multiple celltype_keys (shouldn't happen by construction but
# possible if config is mis-edited) appear once per compartment.
#
# Outputs:
#   outputs/tables/eye/compartment_assignment.csv
.write_compartment_assignment <- function(eye_obj, cfg, cmp_names) {
  meta <- eye_obj@meta.data
  cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                   "knn.leiden.cluster" else "seurat_clusters"
  ctb_col <- resolve_celltype_broad(meta)

  rows <- lapply(cmp_names, function(cmp) {
    cmp_cfg <- cfg$compartments[[cmp]]
    if (!is.null(cmp_cfg$celltype_keys) && !is.null(ctb_col)) {
      keys <- as.character(cmp_cfg$celltype_keys)
      hit  <- as.character(meta[[ctb_col]]) %in% keys
      method <- paste0("celltype_keys(", ctb_col, ")")
    } else {
      parents <- as.character(cmp_cfg$parent_clusters)
      hit     <- as.character(meta[[cluster_col]]) %in% parents
      method  <- "parent_clusters"
    }
    cells <- colnames(eye_obj)[hit]
    if (length(cells) == 0L) return(NULL)
    data.frame(
      cell_id          = cells,
      compartment      = cmp,
      selection_method = method,
      celltype_broad   = if (!is.null(ctb_col)) as.character(meta[cells, ctb_col])
                         else                    NA_character_,
      eye_cluster      = as.character(meta[cells, cluster_col]),
      stringsAsFactors = FALSE
    )
  })
  assigned <- do.call(rbind, rows[!sapply(rows, is.null)])

  # Cells not landing in any compartment (NK / Eryth / Platelet / lens-
  # filtered / orphan celltypes) → compartment = "none".
  unassigned_cells <- setdiff(colnames(eye_obj), assigned$cell_id %||% character(0))
  if (length(unassigned_cells) > 0L) {
    unassigned <- data.frame(
      cell_id          = unassigned_cells,
      compartment      = "none",
      selection_method = NA_character_,
      celltype_broad   = if (!is.null(ctb_col)) as.character(meta[unassigned_cells, ctb_col])
                         else                    NA_character_,
      eye_cluster      = as.character(meta[unassigned_cells, cluster_col]),
      stringsAsFactors = FALSE
    )
    assigned <- rbind(assigned, unassigned)
  }

  paths_eye <- get_target_paths(cfg, "eye")
  ensure_dir(paths_eye$results_tables)
  out_path <- file.path(paths_eye$results_tables, "compartment_assignment.csv")
  write.csv(assigned, out_path, row.names = FALSE)
  log_message(sprintf("Wrote %s (%d rows)", out_path, nrow(assigned)))

  # Summary log: compartment x celltype_broad cross-tab.
  if (!is.null(ctb_col)) {
    summary_tbl <- as.data.frame.matrix(
      table(assigned$compartment, assigned$celltype_broad, useNA = "ifany"))
    sum_path <- file.path(paths_eye$results_tables,
                          "compartment_assignment_summary.csv")
    write.csv(summary_tbl, sum_path)
    log_message("Wrote ", sum_path)
  }
  invisible(assigned)
}

# Per-sample fastMNN integration. Mirrors the pattern used in R/21_integrate_eye.R
# (per-sample SCE list + batchelor::fastMNN with cfg$compartment_integration$k_mnn).
# Compartment subsets often have very few cells from some samples (e.g. an
# eye sample may have only 1-4 myeloid cells); these tiny per-sample subsets
# break SelectIntegrationFeatures (FindVariableFeatures cannot fit a mean-var
# trend on <30 cells), so they get dropped before integration.
.compartment_fastmnn <- function(obj, cfg) {
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  min_cells <- cfg$compartment_integration$min_cells_per_sample %||% 30L

  sample_counts <- table(obj$orig.ident)
  small <- names(sample_counts)[sample_counts < min_cells]
  if (length(small) > 0) {
    n_drop <- sum(sample_counts[small])
    log_message(sprintf("  Dropping %d samples with <%d cells (%d cells total)",
                        length(small), min_cells, n_drop))
    keep_cells <- colnames(obj)[!obj$orig.ident %in% small]
    obj <- subset(obj, cells = keep_cells)
  }
  samples <- unique(obj$orig.ident)
  log_message(sprintf("  fastMNN on %d samples (%d cells)",
                      length(samples), ncol(obj)))

  SeuratList <- lapply(samples, function(s) {
    cells <- colnames(obj)[obj$orig.ident == s]
    subset(obj, cells = cells)
  })
  names(SeuratList) <- samples

  features <- SelectIntegrationFeatures(SeuratList,
                                        nfeatures = cfg$compartment_integration$nfeatures %||% 2000,
                                        verbose = FALSE)
  features <- scRepertoire::quietVDJgenes(features)
  features <- features[!grepl("^MT-|^RPL|^RPS|^HSP", features, ignore.case = TRUE)]
  log_message(sprintf("  features retained: %d", length(features)))

  sce_list <- lapply(SeuratList, function(x) {
    x <- subset(x, features = features)
    x <- as.SingleCellExperiment(x, assay = "RNA")
    scuttle::logNormCounts(x)
  })
  rm(SeuratList); gc()

  set_global_seed(cfg$seed)
  mnn <- batchelor::fastMNN(sce_list,
                            d = cfg$integration$pca_dims %||% 50,
                            k = cfg$compartment_integration$k_mnn %||% 20)
  emb <- as.matrix(SingleCellExperiment::reducedDim(mnn, "corrected"))
  obj[["fastMNN"]] <- CreateDimReducObject(embeddings = emb, key = "mnn")
  list(obj = obj, emb = emb)
}

# PCA + Harmony fallback for compartments below the harmony_fallback_threshold.
# Harmony works with fewer cells than fastMNN's per-batch requirement.
.compartment_harmony <- function(obj, cfg) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("harmony package not installed; cannot fall back from fastMNN. ",
         "Install harmony or raise compartment_integration$harmony_fallback_threshold.")
  }
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = cfg$compartment_integration$nfeatures %||% 2000,
                              verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = cfg$integration$pca_dims %||% 50, verbose = FALSE)
  set_global_seed(cfg$seed)
  obj <- harmony::RunHarmony(obj, group.by.vars = "orig.ident", verbose = FALSE)
  emb <- Embeddings(obj, "harmony")
  list(obj = obj, emb = emb)
}

# Re-integrate one compartment, choosing fastMNN or Harmony based on cell count.
.integrate_compartment <- function(obj, cfg, cmp) {
  thresh <- cfg$compartment_integration$harmony_fallback_threshold %||% 3000
  if (ncol(obj) >= thresh) {
    log_message("  ", cmp, " has ", ncol(obj), " cells; using fastMNN.")
    .compartment_fastmnn(obj, cfg)
  } else {
    log_message("  ", cmp, " has ", ncol(obj), " cells (below ", thresh,
                "); falling back to PCA + Harmony.")
    .compartment_harmony(obj, cfg)
  }
}

# Resolve clustering config for a compartment: shallow-merge
# `cfg$compartment_clustering` <- `cfg$<cmp>_clustering`. Mirrors the
# eye_clustering pattern in 21_integrate_eye. Score weights live under
# cfg$clustering$score_weights and can be overridden in <cmp>_clustering.
.resolve_compartment_clustering <- function(cfg, cmp) {
  base <- cfg$compartment_clustering %||% list()
  override <- cfg[[paste0(cmp, "_clustering")]] %||% list()
  merged <- utils::modifyList(base, override)
  # Score weights: own merge so user can override individual weight keys.
  base_w <- cfg$clustering$score_weights %||% list(
    stability = 0.25, modularity = 0.25,
    silhouette = 0.30, singletons = 0.20)
  merged$score_weights <- utils::modifyList(base_w, merged$score_weights %||% list())
  merged
}

# Tune Leiden clustering with the compartment-scoped grid and write the
# tuning summary. Returns the labels and the chosen k.
.cluster_compartment <- function(obj, emb, cfg, paths, cmp) {
  cl_cfg <- .resolve_compartment_clustering(cfg, cmp)
  log_message(sprintf(
    "  tuning Leiden (%s) over k=[%s], res=[%.2f, %.2f], clusters=[%d, %d], silhouette_w=%.2f",
    cmp,
    paste(cl_cfg$k_grid, collapse = ","),
    cl_cfg$res_grid_min, cl_cfg$res_grid_max,
    as.integer(cl_cfg$min_clusters), as.integer(cl_cfg$max_clusters),
    cl_cfg$score_weights$silhouette))

  knn_leiden <- tune_knn_leiden(
    emb           = emb,
    k_grid        = as.integer(cl_cfg$k_grid),
    res_grid      = seq(cl_cfg$res_grid_min, cl_cfg$res_grid_max, by = cl_cfg$res_grid_step),
    repeats       = as.integer(cfg$clustering$repeats %||% 3L),
    min_clusters  = as.integer(cl_cfg$min_clusters),
    max_clusters  = as.integer(cl_cfg$max_clusters),
    score_weights = cl_cfg$score_weights,
    approximate   = TRUE
  )

  ensure_dir(paths$results_tables)
  write.csv(knn_leiden$summary,
            file.path(paths$results_tables, "knn_leiden_tuning.csv"),
            row.names = FALSE)

  best_lab <- knn_leiden$best$labels
  na_idx <- which(is.na(best_lab))
  if (length(na_idx) > 0) {
    valid_idx <- which(!is.na(best_lab))
    cls_valid <- sort(unique(best_lab[valid_idx]))
    centroids <- t(sapply(cls_valid, function(cl)
      colMeans(emb[best_lab == cl & !is.na(best_lab), , drop = FALSE])))
    for (ci in na_idx) {
      d <- apply(centroids, 1, function(ctr) sqrt(sum((emb[ci, ] - ctr)^2)))
      best_lab[ci] <- cls_valid[which.min(d)]
    }
    log_message(sprintf("  Rescued %d NA cells", length(na_idx)))
  }

  obj$knn.leiden.cluster <- as.factor(best_lab)
  obj$seurat_clusters    <- as.factor(best_lab)   # alias for downstream code
  Idents(obj) <- "knn.leiden.cluster"

  umap_tuned <- tune_umap(
    emb          = emb,
    labels       = best_lab,
    nn_grid      = c(knn_leiden$best$k,
                     max(10L, round(knn_leiden$best$k * 0.5)),
                     round(knn_leiden$best$k * 1.5)),
    mindist_grid = c(0.05, 0.1, 0.2, 0.4)
  )
  write.csv(umap_tuned$grid,
            file.path(paths$results_tables, "umap_tuning.csv"),
            row.names = FALSE)

  U <- uwot::umap(emb,
                  n_neighbors = umap_tuned$best$n_neighbors %||% 15L,
                  min_dist    = umap_tuned$best$min_dist    %||% 0.1,
                  metric      = "euclidean",
                  fast_sgd    = TRUE,
                  ret_model   = FALSE,
                  verbose     = FALSE)
  colnames(U) <- paste0("UMAP_", 1:2)
  obj[["UMAP"]] <- CreateDimReducObject(
    embeddings = U,
    assay      = "RNA",
    stdev      = numeric(),
    key        = "UMAP_"
  )

  obj@misc$knn_leiden <- list(
    best_k          = knn_leiden$best$k,
    best_resolution = knn_leiden$best$resolution,
    summary         = knn_leiden$summary,
    snn_graph       = knn_leiden$best$graph
  )
  obj@misc$umap_tuned <- umap_tuned$grid

  obj
}

# Post-hoc cluster-level contamination filter (T1.4). After initial Leiden
# clustering on a compartment, scores each cluster's per-gene detection rate
# on a panel of expected (canonical lineage) and off-target (contamination)
# marker sets. Flags a cluster as contamination when:
#   max(off-target panel score) > expected panel score, AND
#   max(off-target panel score) >= cfg$compartment_contamination_filter$min_offtarget
# Flagged cells are dropped from the compartment object; the caller then
# re-clusters on the existing embedding restricted to kept cells.
#
# Config under cfg$compartment_contamination_filter:
#   enable: bool (default FALSE)
#   min_offtarget: numeric (default 0.20) — only flag clusters where the
#       dominant off-target panel reaches this fraction of genes detected.
#   <cmp>:                                  # per-compartment override block
#     enable: bool (default TRUE if parent enable is TRUE)
#     expected_panel: [genes…]              # canonical lineage panel
#     offtarget_panels: { name: [genes…] }  # one or more contamination panels
#
# Audit columns written under outputs/tables/qc/compartment_contamination_<cmp>.csv:
#   cluster | n_cells | expected_score | <off>_score (one column per off panel) |
#   off_max | off_max_panel | drop | reason
.identify_contamination_clusters <- function(obj, cmp, cfg) {
  fcfg <- cfg$compartment_contamination_filter
  if (is.null(fcfg) || !isTRUE(fcfg$enable)) {
    return(list(drop_clusters = character(0), audit = NULL))
  }
  cmp_cfg <- fcfg[[cmp]] %||% list()
  if (isFALSE(cmp_cfg$enable %||% TRUE)) {
    return(list(drop_clusters = character(0), audit = NULL))
  }
  expected_panel  <- as.character(cmp_cfg$expected_panel %||% character())
  offtarget_panels <- cmp_cfg$offtarget_panels %||% list()
  if (length(expected_panel) == 0L || length(offtarget_panels) == 0L) {
    log_message(sprintf(
      "Contamination filter %s: missing expected_panel or offtarget_panels; skipping.",
      cmp))
    return(list(drop_clusters = character(0), audit = NULL))
  }
  min_offtarget <- as.numeric(fcfg$min_offtarget %||% 0.20)

  DefaultAssay(obj) <- "RNA"
  e <- GetAssayData(obj, assay = "RNA", layer = "data")
  cl <- as.character(obj$knn.leiden.cluster)
  clusters <- sort(unique(cl))

  score_panel <- function(panel) {
    g <- intersect(as.character(panel), rownames(e))
    if (length(g) == 0L) return(setNames(rep(NA_real_, length(clusters)), clusters))
    vapply(clusters, function(k) {
      cells <- which(cl == k)
      mean(Matrix::rowMeans(e[g, cells, drop = FALSE] > 0))
    }, numeric(1))
  }

  exp_scores <- score_panel(expected_panel)
  off_scores <- sapply(offtarget_panels, score_panel)
  if (!is.matrix(off_scores)) {
    off_scores <- matrix(off_scores, ncol = 1,
                         dimnames = list(clusters, names(offtarget_panels)))
  }

  off_max_idx   <- max.col(off_scores, ties.method = "first")
  off_max_score <- off_scores[cbind(seq_along(off_max_idx), off_max_idx)]
  off_max_panel <- colnames(off_scores)[off_max_idx]

  drop <- (off_max_score > exp_scores) & (off_max_score >= min_offtarget)
  reason <- ifelse(drop,
                   sprintf("off-target %s=%.2f > expected=%.2f",
                           off_max_panel, off_max_score, exp_scores),
                   "")

  audit <- data.frame(
    cluster        = clusters,
    n_cells        = as.integer(table(cl)[clusters]),
    expected_score = round(exp_scores, 3),
    off_scores,
    off_max        = round(off_max_score, 3),
    off_max_panel  = off_max_panel,
    drop           = drop,
    reason         = reason,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  drop_clusters <- clusters[drop]
  if (length(drop_clusters) > 0L) {
    log_message(sprintf(
      "Contamination filter %s: dropping cluster(s) %s — %s",
      cmp,
      paste(drop_clusters, collapse = ","),
      paste(reason[drop], collapse = "; ")))
  } else {
    log_message(sprintf("Contamination filter %s: no clusters flagged.", cmp))
  }

  list(drop_clusters = drop_clusters, audit = audit)
}

# Subset, integrate, cluster a single compartment. Writes the resulting object
# under outputs/objects/eye/<cmp>/IntegratedSeuratObject.rds and a marker table
# at outputs/tables/eye/<cmp>/<cmp>_substate_markers.csv (markers are run by
# 30_markers.R; subset_one_compartment does not run markers itself).
subset_one_compartment <- function(eye_obj, cmp, cfg) {
  set_global_seed(cfg$seed)
  cmp_cfg <- cfg$compartments[[cmp]]
  meta    <- eye_obj@meta.data

  # Prefer celltype_keys (stable across re-clusterings). Resolve the
  # celltype_broad source column with the _full carry-over fallback so
  # eye objects built without eye_annotate still work.
  ctb_col <- resolve_celltype_broad(meta)
  use_celltype <- !is.null(cmp_cfg$celltype_keys) && !is.null(ctb_col)

  if (use_celltype) {
    keys  <- as.character(cmp_cfg$celltype_keys)
    cells <- colnames(eye_obj)[as.character(meta[[ctb_col]]) %in% keys]
    log_message(sprintf("Compartment %s: %d cells via celltype_keys [%s] on `%s`",
                        cmp, length(cells), paste(keys, collapse = ","), ctb_col))
  } else {
    parents <- as.character(cmp_cfg$parent_clusters)
    cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                     "knn.leiden.cluster" else "seurat_clusters"
    cells <- colnames(eye_obj)[as.character(meta[[cluster_col]]) %in% parents]
    log_message(sprintf("Compartment %s: %d cells via parent_clusters [%s] (LEGACY)",
                        cmp, length(cells), paste(parents, collapse = ",")))
  }

  # Canonical-marker lineage gate (R/16_compartment_lineage_gate.R). Drops
  # cells whose top scoring panel is outside the compartment's expected
  # lineages or that co-express any non-immune panel. No-op when
  # cfg$compartment_lineage_gate$enable is FALSE.
  gate_res <- apply_compartment_lineage_gate(eye_obj, cells, cmp, cfg)
  cells <- gate_res$cells
  if (!is.null(gate_res$audit)) {
    .gate_audit_buffer[[cmp]] <<- gate_res$audit
  }

  # Myeloid CTaa filter — config-driven dispatch via R/14_repertoire_qc.R.
  # Selects one of options (a/b/c/d) from cfg$compartments$myeloid_ctaa_filter.
  # Default "a" preserves the legacy behavior (drop any CTaa+).
  ctaa_meta <- NULL
  if (cmp == "myeloid") {
    res <- apply_myeloid_ctaa_filter(eye_obj, cells, cfg)
    cells <- res$cells
    ctaa_meta <- res$tagged   # non-NULL only for option "c"
  }

  if (length(cells) < 100) {
    log_message("  Compartment too small to integrate; skipping.")
    return(invisible(NULL))
  }

  sub <- subset(eye_obj, cells = cells)
  sub <- .carry_forward_eye_metadata(sub)

  paths <- get_target_paths(cfg, cmp)
  ensure_dir(paths$results_objects)
  ensure_dir(paths$results_tables)

  integrated <- .integrate_compartment(sub, cfg, cmp)
  sub <- integrated$obj
  emb <- integrated$emb

  sub <- .cluster_compartment(sub, emb, cfg, paths, cmp)

  # Post-hoc contamination filter (T1.4). When enabled in
  # cfg$compartment_contamination_filter, drop entire clusters whose top
  # off-target marker panel beats the canonical expected panel. After
  # dropping, re-cluster on the existing embedding restricted to kept cells
  # so the surviving cluster IDs renumber cleanly (no gaps).
  contam <- .identify_contamination_clusters(sub, cmp, cfg)
  if (!is.null(contam$audit)) {
    qc_dir <- "outputs/tables/qc"
    ensure_dir(qc_dir)
    audit_path <- file.path(qc_dir,
                            sprintf("compartment_contamination_%s.csv", cmp))
    write.csv(contam$audit, audit_path, row.names = FALSE)
    log_message("Wrote ", audit_path)
  }
  if (length(contam$drop_clusters) > 0L) {
    keep_mask <- !as.character(sub$knn.leiden.cluster) %in% contam$drop_clusters
    keep_cells <- colnames(sub)[keep_mask]
    n_before <- ncol(sub)
    sub <- subset(sub, cells = keep_cells)
    emb <- emb[keep_cells, , drop = FALSE]
    log_message(sprintf(
      "  Contamination filter: %d cells dropped (cluster(s) %s); re-clustering %d cells.",
      n_before - ncol(sub),
      paste(contam$drop_clusters, collapse = ","),
      ncol(sub)))
    sub <- .cluster_compartment(sub, emb, cfg, paths, cmp)
  }

  # Stamp compartment-aware substate keys (R/23_substate_labels.R). Downstream
  # modules (R/45_compartment_pca, R/47_liana, R/65 joint-substate) read these
  # columns rather than reconstructing labels from cluster IDs.
  if (exists("apply_substate_keys", mode = "function")) {
    sub <- apply_substate_keys(sub, cmp, cfg)
  }

  saveRDS(sub, file.path(paths$results_objects, "IntegratedSeuratObject.rds"))
  log_message(sprintf("  Wrote %s with %d cells, %d substates",
                      file.path(paths$results_objects, "IntegratedSeuratObject.rds"),
                      ncol(sub),
                      length(unique(sub$knn.leiden.cluster))))
  invisible(sub)
}

# Per-compartment audit rows accumulated by apply_compartment_lineage_gate()
# during subset_one_compartment() and flushed at the end of subset_compartments().
# Module-level list so the gate writer doesn't need to thread the buffer
# through .integrate_compartment / .cluster_compartment.
.gate_audit_buffer <- list()

# Top-level entry point called from run_pipeline.R Phase 1c.
subset_compartments <- function(cfg) {
  log_message("=== Phase 1c: compartment subsetting ===")
  .gate_audit_buffer <<- list()
  paths_eye <- get_target_paths(cfg, "eye")
  eye_path  <- file.path(paths_eye$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(eye_path)) {
    log_message("Eye IntegratedSeuratObject.rds not found at ", eye_path,
                ". Run subset_eye and reintegrate_eye first.")
    return(invisible(FALSE))
  }
  eye_obj <- readRDS(eye_path)
  log_message("Loaded eye object: ", ncol(eye_obj), " cells.")

  # Per-cell UCell lens filter (R/15_lens_qc.R). No-op when
  # cfg$compartments$lens_ucell_filter$enable is FALSE; errors when enabled
  # without a threshold so the user is forced to run the diagnostic first.
  eye_obj <- .apply_lens_ucell_filter(eye_obj, cfg)

  # T2.5: timepoint policy. Scoped to the compartment branch so eye/full
  # analyses keep their full per-cell view. Default "all"; "earliest"
  # keeps the lexically lowest orig.ident per Subject (works for sample
  # IDs encoding visit number, e.g. UV215_003_AC < UV215_004_AC).
  eye_obj <- .apply_compartment_timepoint_policy(eye_obj, cfg)

  # Iterate compartment names only. cfg$compartments also holds non-
  # compartment knobs (lens_ucell_filter, myeloid_ctaa_filter,
  # timepoint_policy); skip anything that isn't a list with at least one
  # selection key. Keeps the loop robust as new top-level knobs are added.
  is_compartment <- vapply(cfg$compartments, function(x) {
    is.list(x) && (!is.null(x$celltype_keys) || !is.null(x$parent_clusters))
  }, logical(1))
  cmp_names <- names(cfg$compartments)[is_compartment]

  # Audit trail: write compartment_assignment.csv mapping every (post-lens-
  # filter) eye cell to its compartment label. This makes overlaps and
  # orphan cells discoverable without rereading the compartment objects.
  .write_compartment_assignment(eye_obj, cfg, cmp_names)

  for (cmp in cmp_names) {
    subset_one_compartment(eye_obj, cmp, cfg)
  }

  # Flush the consolidated lineage-gate audit. Non-empty only when
  # cfg$compartment_lineage_gate$enable == TRUE.
  if (length(.gate_audit_buffer) > 0L) {
    audit <- do.call(rbind, .gate_audit_buffer)
    rownames(audit) <- NULL
    qc_dir <- "outputs/tables/qc"
    ensure_dir(qc_dir)
    out_path <- file.path(qc_dir, "compartment_gate_report.csv")
    write.csv(audit, out_path, row.names = FALSE)
    log_message(sprintf(
      "Wrote %s (%d rows; %d dropped overall)",
      out_path, nrow(audit), sum(!audit$pass)))

    summary_tbl <- audit %>%
      dplyr::group_by(compartment) %>%
      dplyr::summarise(
        n_in     = dplyr::n(),
        n_pass   = sum(pass),
        n_drop   = sum(!pass),
        pct_pass = round(100 * mean(pass), 1),
        .groups  = "drop"
      )
    summary_path <- file.path(qc_dir, "compartment_gate_summary.csv")
    write.csv(summary_tbl, summary_path, row.names = FALSE)
    log_message("Wrote ", summary_path)
  }

  log_message("=== Compartment subsetting complete ===")
  invisible(TRUE)
}

# subset_compartments reads the eye sub-atlas, applies the UCell-based per-
# cell lens contamination filter from R/15_lens_qc.R (no-op when disabled),
# then for each compartment subsets the parent clusters, re-integrates with
# fastMNN (or Harmony for compartments below the cell-count threshold),
# tunes Leiden clustering on the compartment-scoped grid, computes a UMAP,
# and writes the Seurat object under
# outputs/objects/eye/<cmp>/IntegratedSeuratObject.rds. For the myeloid
# compartment, any cell with a non-NA CTaa is dropped before integration to
# strip lymphocyte doublets / mis-clustered T/B cells. The function assumes
# the eye object already exists and treats it as read-only. Eye-level
# metadata (cluster IDs, celltype labels, escape.UCell assay) is carried
# forward on each compartment object under *_eye suffixes so downstream
# panels can show inheritance.
