# R/21_integrate_eye.R
# Re-derive HVG -> fastMNN -> kNN-Leiden -> UMAP on the eye-only subset.
# Reuses tune_knn_leiden, tune_umap, validate_clusters_by_lineage from
# 10_integrate_full.R. Does not re-run scDblFinder or reload per-sample
# rds files.
suppressPackageStartupMessages({
  library(Seurat)
  library(scRepertoire)
  library(SingleCellExperiment)
  library(batchelor)
  library(igraph)
  library(uwot)
  library(Matrix)
})

reintegrate_eye <- function(cfg) {
  if (!isTRUE(cfg$eye_focus$enable)) {
    log_message("eye_focus disabled. Skipping eye reintegration.")
    return(invisible(TRUE))
  }

  paths <- get_target_paths(cfg, "eye")
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Eye subset not found at ", obj_path, ". Run subset_eye() first.")
    return(invisible(FALSE))
  }

  log_message("Loading eye subset for re-integration...")
  obj <- readRDS(obj_path)

  # Eye-specific overrides shallow-merged on top of shared config
  int_cfg <- utils::modifyList(cfg$integration %||% list(),
                               cfg$eye_integration %||% list())
  cl_cfg  <- utils::modifyList(cfg$clustering  %||% list(),
                               cfg$eye_clustering  %||% list())
  log_message(sprintf(
    "Eye integration params: k_mnn=%s, pca_dims=%s | clustering: k_grid=[%s], min/max=%s/%s, res=[%s, %s]",
    int_cfg$k_mnn, int_cfg$pca_dims,
    paste(cl_cfg$k_grid, collapse = ","),
    cl_cfg$min_clusters, cl_cfg$max_clusters,
    cl_cfg$res_grid_min, cl_cfg$res_grid_max))

  # --- Per-sample split for fastMNN ---
  # Use orig.ident as the batch (matches 10_integrate_full.R behavior).
  log_message("Splitting eye object by orig.ident for fastMNN...")
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  samples <- unique(obj$orig.ident)
  log_message(sprintf("  %d samples in eye object", length(samples)))

  SeuratList <- lapply(samples, function(s) {
    cells <- colnames(obj)[obj$orig.ident == s]
    subset(obj, cells = cells)
  })
  names(SeuratList) <- samples

  # --- HVG selection ---
  features <- SelectIntegrationFeatures(SeuratList,
                                        nfeatures = cfg$integration$nfeatures,
                                        verbose = FALSE)
  features <- scRepertoire::quietVDJgenes(features)
  exclude_patterns <- "^MT-|^RPL|^RPS|^HSP"
  features <- features[!grepl(exclude_patterns, features, ignore.case = TRUE)]
  n_keep <- min(2000L, length(features))
  features <- features[seq_len(n_keep)]
  log_message(sprintf("Eye integration features: %d retained", length(features)))

  sce_list <- lapply(SeuratList, function(x) {
    x <- subset(x, features = features)
    x <- as.SingleCellExperiment(x, assay = "RNA")
    scuttle::logNormCounts(x)
  })
  rm(SeuratList); gc()

  # --- fastMNN ---
  log_message("Running fastMNN on eye subset...")
  set_global_seed(cfg$seed)
  mnn <- batchelor::fastMNN(sce_list,
                            d = int_cfg$pca_dims,
                            k = int_cfg$k_mnn)
  rm(sce_list); gc()
  emb <- as.matrix(SingleCellExperiment::reducedDim(mnn, "corrected"))
  obj[["fastMNN"]] <- CreateDimReducObject(embeddings = emb, key = "mnn")

  # --- kNN-Leiden tuning (reuses helper from 10_integrate_full.R) ---
  log_message("Tuning kNN-Leiden clustering on eye subset...")
  knn_leiden <- tune_knn_leiden(
    emb           = emb,
    k_grid        = as.integer(cl_cfg$k_grid %||% c(15L, 20L, 30L)),
    res_grid      = seq(cl_cfg$res_grid_min %||% 0.1,
                        cl_cfg$res_grid_max %||% 1.0,
                        by = cl_cfg$res_grid_step %||% 0.05),
    repeats       = as.integer(cl_cfg$repeats %||% 3L),
    min_clusters  = as.integer(cl_cfg$min_clusters %||% 12L),
    max_clusters  = as.integer(cl_cfg$max_clusters %||% 20L),
    score_weights = cl_cfg$score_weights,
    approximate   = TRUE
  )

  best_k   <- knn_leiden$best$k
  best_res <- knn_leiden$best$resolution
  best_lab <- knn_leiden$best$labels

  ensure_dir(paths$results_tables)
  write.csv(knn_leiden$summary,
            file.path(paths$results_tables, "knn_leiden_tuning.csv"),
            row.names = FALSE)

  obj <- NormalizeData(obj, verbose = FALSE)

  # --- Lineage validation (reassign cross-lineage cells) ---
  log_message("Validating eye clusters for cross-lineage mixing...")
  best_lab <- validate_clusters_by_lineage(
    obj           = obj,
    emb           = emb,
    labels        = best_lab,
    mix_threshold = 0.10
  )

  # NA rescue
  na_idx <- which(is.na(best_lab))
  if (length(na_idx) > 0L) {
    valid_idx <- which(!is.na(best_lab))
    cls_valid <- sort(unique(best_lab[valid_idx]))
    na_centroids <- t(sapply(cls_valid, function(cl)
      colMeans(emb[best_lab == cl & !is.na(best_lab), , drop = FALSE])))
    for (ci in na_idx) {
      d <- apply(na_centroids, 1, function(ctr) sqrt(sum((emb[ci, ] - ctr)^2)))
      best_match <- which.min(d)
      if (length(best_match) > 0L) best_lab[ci] <- cls_valid[best_match]
    }
    log_message(sprintf("Rescued %d NA cells", length(na_idx)))
  }

  obj$knn.leiden.cluster <- as.factor(best_lab)
  Seurat::Idents(obj) <- "knn.leiden.cluster"

  # --- UMAP tuning ---
  log_message("Tuning UMAP for eye subset...")
  umap_tuned <- tune_umap(
    emb          = emb,
    labels       = best_lab,
    nn_grid      = c(best_k, max(10L, round(best_k * 0.5)), round(best_k * 1.5)),
    mindist_grid = c(0.05, 0.1, 0.2, 0.4)
  )

  write.csv(umap_tuned$grid,
            file.path(paths$results_tables, "umap_tuning.csv"),
            row.names = FALSE)

  U <- uwot::umap(
    emb,
    n_neighbors = umap_tuned$best$n_neighbors %||% 15L,
    min_dist    = umap_tuned$best$min_dist    %||% 0.1,
    metric      = "euclidean",
    fast_sgd    = TRUE,
    ret_model   = FALSE,
    verbose     = FALSE
  )
  colnames(U) <- paste0("UMAP_", 1:2)
  obj[["UMAP"]] <- CreateDimReducObject(
    embeddings = U,
    assay      = "RNA",
    stdev      = numeric(),
    key        = "UMAP_"
  )

  obj@misc$knn_leiden <- list(
    best_k          = best_k,
    best_resolution = best_res,
    summary         = knn_leiden$summary,
    snn_graph       = knn_leiden$best$graph
  )
  obj@misc$umap_tuned <- umap_tuned$grid

  saveRDS(obj, file.path(paths$results_objects, "IntegratedSeuratObject.rds"))
  log_message("Saved re-integrated eye object")
  invisible(TRUE)
}
