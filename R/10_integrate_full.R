# R/10_integrate_full.R
suppressPackageStartupMessages({
  library(Seurat)
  library(scRepertoire)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scDblFinder)
  library(batchelor)
  library(igraph)
  library(batchelor)
  library(bluster)
  library(leidenAlg)
  library(dplyr)
  library(stringr)
  library(mclust)         
  library(cluster)       
  library(BiocNeighbors)  
  library(Matrix)
  library(uwot)           
})

integrate_fastmnn <- function(cfg) {
  
  #Isolating files with meta data associated
  files <- list.files(cfg$paths$processed_dir, pattern="\\.rds$", full.names=TRUE)
  meta <- read.csv("inputs/data/metadata.csv", row.names = 1, stringsAsFactors = FALSE)
  # Run-folder name = two dirname levels up from the CellRanger output path
  # (matches 10_ingest_data.R's `basename(dirname(dirname(...)))` used when
  # building the QC metadata join). Previously we extracted the 10th path
  # component, which breaks whenever absolute-path depth changes.
  meta$run <- basename(dirname(dirname(meta$RNA_CRoutput)))
  files <- files[grep(paste0(meta$run, collapse = "|"), files)]
  
  log_message("Loading and merging data...")
  SeuratList <- lapply(files, function(x) {
    tmp <- readRDS(x) 
    DefaultAssay(tmp) <- "RNA"
    tmp[["SCT"]] <- NULL
    tmp
  })
  SeuratObj <- merge(SeuratList[[1]], SeuratList[-1])
  SeuratObj <- JoinLayers(SeuratObj)
  
  #Filtering Doublets One Last Time
  sce_merged <- as.SingleCellExperiment(SeuratObj)
  sce_merged <- scDblFinder(sce_merged, 
                            samples = "orig.ident",  # one model per capture
                            dims = 50,
                            dbr.sd = 1,
                            BPPARAM = BiocParallel::MulticoreParam(4))
  
  cells <- colnames(sce_merged)[sce_merged$scDblFinder.class == "singlet"]
  SeuratObj <- subset(SeuratObj, cells = cells)
  rm(sce_merged)

  # Subset SeuratList to match filtered cells (remove post-merge doublets)
  SeuratList <- lapply(SeuratList, function(x) {
    keep <- intersect(Cells(x), cells)
    subset(x, cells = keep)
  })

  # Extracting Stable Variable Features
  features <- SelectIntegrationFeatures(SeuratList,
                                        nfeatures = cfg$integration$nfeatures,
                                        verbose = FALSE)
  
  # Removing Complicating Variable Features
  features <- scRepertoire::quietVDJgenes(features)
  exclude_patterns <- "^MT-|^RPL|^RPS|^HSP"
  features <- features[!grepl(exclude_patterns, features, ignore.case = TRUE)]
  # Take up to 2000; avoid NAs when fewer remain after exclusion filters.
  n_keep   <- min(2000L, length(features))
  features <- features[seq_len(n_keep)]
  log_message(sprintf("Integration features: %d retained after filtering",
                      length(features)))

  sce_list <- lapply(SeuratList, function(x) {
    x <- subset(x, features = features)
    x <- as.SingleCellExperiment(x, assay = "RNA")
    scuttle::logNormCounts(x)
  })
  rm(SeuratList); gc()
  
  log_message("Running fastMNN on variable features...")
  set_global_seed(cfg$seed)
  mnn <- batchelor::fastMNN(sce_list,
                            d = cfg$integration$pca_dims,
                            k = cfg$integration$k_mnn)
  rm(sce_list); gc()
  emb <- as.matrix(SingleCellExperiment::reducedDim(mnn, "corrected"))
  SeuratObj[["fastMNN"]] <- CreateDimReducObject(embeddings = emb, 
                                                 key = "mnn")
  
  log_message("Calculating Optimal Leiden Clustering...")
  cl_cfg <- cfg$clustering %||% list()
  knn_leiden <- tune_knn_leiden(
    emb          = emb,
    k_grid       = as.integer(cl_cfg$k_grid %||% c(15L, 20L, 30L)),
    res_grid     = seq(cl_cfg$res_grid_min %||% 0.1,
                       cl_cfg$res_grid_max %||% 1.0,
                       by = cl_cfg$res_grid_step %||% 0.05),
    repeats      = as.integer(cl_cfg$repeats %||% 3L),
    min_clusters = as.integer(cl_cfg$min_clusters %||% 12L),
    max_clusters = as.integer(cl_cfg$max_clusters %||% 20L),
    score_weights = cl_cfg$score_weights,
    approximate  = TRUE
  )
  
  best_k   <- knn_leiden$best$k
  best_res <- knn_leiden$best$resolution
  best_lab <- knn_leiden$best$labels
  
  write.csv(knn_leiden$summary,
            file.path(cfg$paths$results_tables, "knn_leiden_tuning.csv"),
            row.names = FALSE)

  SeuratObj <- NormalizeData(SeuratObj, verbose = FALSE)

  log_message("Validating clusters for cross-lineage mixing...")
  best_lab <- validate_clusters_by_lineage(
    obj            = SeuratObj,
    emb            = emb,
    labels         = best_lab,
    mix_threshold  = 0.10
  )

  # Rescue any cells with NA cluster labels (e.g. degenerate embeddings)
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
    log_message(sprintf("Rescued %d NA cells by nearest-centroid assignment", length(na_idx)))
  }

  SeuratObj$knn.leiden.cluster <- as.factor(best_lab)
  Seurat::Idents(SeuratObj) <- "knn.leiden.cluster"
  
  log_message("Calculating Optimal UMAP Embeddings...")
  umap_tuned <- tune_umap(
    emb         = emb,
    labels      = best_lab,
    nn_grid     = c(best_k, max(10L, round(best_k*0.5)), round(best_k*1.5)),
    mindist_grid= c(0.05, 0.1, 0.2, 0.4)
  )
  
  write.csv(umap_tuned$grid, 
            file.path(cfg$paths$results_tables, "umap_tuning.csv"), 
            row.names = FALSE)
  
  # Use the tuned UMAP params; fall back to uwot defaults (15 / 0.1) if
  # tuning produced NA. Previously this referenced cfg$umap, which doesn't
  # exist in config.yml — the fallback was a no-op and there was also a
  # stray double-comma making the call syntactically suspect.
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
  
  SeuratObj[["UMAP"]] <- CreateDimReducObject(
                              embeddings = U,
                              assay = "RNA",
                              stdev = numeric(),
                              key = "UMAP_")
  
  #Store in misc (including SNN graph for downstream refinement)
  SeuratObj@misc$knn_leiden <- list(
    best_k = best_k,
    best_resolution = best_res,
    summary = knn_leiden$summary,
    snn_graph = knn_leiden$best$graph
  )
  SeuratObj@misc$umap_tuned <- umap_tuned$grid

  ensure_dir(cfg$paths$results_objects)
  saveRDS(SeuratObj, file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds"))
  invisible(TRUE)
}

# -------- Silhouette helper (with automatic subsampling for speed) --------
.silhouette_mean <- function(labels, X, metric = "euclidean", max_n = 5000L, seed = 1L) {
  labs <- as.integer(as.factor(labels))
  if (length(unique(labs)) < 2L || nrow(X) < 3L) return(NA_real_)
  set.seed(seed)
  if (nrow(X) > max_n) {
    keep <- sample.int(nrow(X), max_n)
    Xsub  <- X[keep, , drop = FALSE]
    labs  <- labs[keep]
  } else {
    Xsub <- X
  }
  # Drop clusters that became singletons after subsampling —
  # silhouette() errors on clusters with < 2 members
  tab <- table(labs)
  valid_cls <- as.integer(names(tab[tab >= 2L]))
  keep2 <- labs %in% valid_cls
  if (sum(keep2) < 3L || length(valid_cls) < 2L) return(NA_real_)
  Xsub <- Xsub[keep2, , drop = FALSE]
  labs  <- as.integer(as.factor(labs[keep2]))   # re-factor to drop empty levels
  d <- dist(Xsub, method = metric)
  sil <- try(silhouette(labs, d), silent = TRUE)
  if (inherits(sil, "try-error")) return(NA_real_)
  mean(sil[, "sil_width"], na.rm = TRUE)
}

# -------- Build a fast SNN graph with bluster + edge pruning
.make_snn_graph <- function(emb, k = 30L, approximate = TRUE, metric = "cosine",
                            prune_snn = 1/15) {
  bnparam <- if (approximate) {
    BiocNeighbors::AnnoyParam(ntrees = 50, search.mult = 2, distance = metric)
  } else {
    BiocNeighbors::KmknnParam(distance = metric)
  }
  # SNN weights edges by shared-neighbor overlap (rank-based Jaccard),
  # reducing spurious bridges between lineages that share gene programs (e.g. ISG)
  g <- bluster::makeSNNGraph(emb, k = k, BNPARAM = bnparam, type = "jaccard")
  if (is.directed(g)) g <- as.undirected(g, mode = "collapse")
  # Prune weak SNN edges — analogous to Seurat FindNeighbors(prune.SNN = 1/15)
  if (prune_snn > 0 && !is.null(igraph::E(g)$weight)) {
    cutoff <- prune_snn * max(igraph::E(g)$weight)
    g <- igraph::delete_edges(g, which(igraph::E(g)$weight < cutoff))
  }
  g
}

# -------- One Leiden run at a resolution on graph g ------------------------
.one_leiden_run <- function(g, resolution, embedding = NULL, seed = 1L) {
  set.seed(seed)
  memb <- leidenAlg::leiden.community(g, resolution = resolution)
  labels <- memb$membership
  # Metrics (use weights if present)
  w <- E(g)$weight
  mod <- if (!is.null(w)) modularity(g, labels, weights = w) else modularity(g, labels)
  k    <- length(unique(labels))
  pct_singletons <- mean(table(labels) == 1) * 100
  sil  <- if (is.null(embedding)) NA_real_ else .silhouette_mean(labels, embedding)
  list(
    labels = labels,
    metrics = data.frame(
      kNN_k = attr(g, "knn_k") %||% NA_integer_,
      resolution = resolution,
      clusters = k,
      modularity = mod,
      mean_silhouette = sil,
      pct_singletons = pct_singletons,
      seed = seed,
      stringsAsFactors = FALSE
    )
  )
}

# -------- Grid search over k and resolution with repeats & stability -------
tune_knn_leiden <- function(
    emb,
    k_grid        = c(15L, 20L, 30L),
    res_grid      = seq(0.1, 1.0, by = 0.05),
    repeats       = 3L,
    seeds         = NULL,
    min_clusters  = 12L,
    max_clusters  = 20L,
    score_weights = NULL,
    metric_knn    = "Cosine",
    metric_sil    = "Cosine",
    approximate   = TRUE
) {
  if (is.null(seeds)) seeds <- seq_len(repeats)

  n <- nrow(emb)
  dropped <- k_grid[k_grid >= n]
  k_grid <- k_grid[k_grid < n]
  if (length(dropped) > 0L)
    message("tune_knn_leiden: dropping k = ", paste(dropped, collapse = ", "),
            " (>= ", n, " data points)")
  if (length(k_grid) == 0L)
    stop("All k_grid values >= number of data points (", n, "). ",
         "Reduce k_grid or check that 'emb' is correct.")

  all_rows <- list()
  per_setting <- list()

  for (k in k_grid) {
    g <- .make_snn_graph(emb, k = k, approximate = approximate, metric = metric_knn)
    attr(g, "knn_k") <- k
    
    for (r in res_grid) {
      run_labels <- vector("list", length(seeds))
      metrics_df <- vector("list", length(seeds))
      
      for (i in seq_along(seeds)) {
        out <- .one_leiden_run(g, r, embedding = emb, seed = seeds[i])
        run_labels[[i]] <- out$labels
        metrics_df[[i]] <- out$metrics
      }
      metrics_df <- do.call(rbind, metrics_df)
      
      # Stability across repeats: mean pairwise ARI
      aris <- c()
      if (length(run_labels) >= 2L) {
        for (i in 1:(length(run_labels) - 1L)) {
          for (j in (i + 1L):length(run_labels)) {
            aris <- c(aris, mclust::adjustedRandIndex(run_labels[[i]], run_labels[[j]]))
          }
        }
      }
      mean_ari <- if (length(aris)) mean(aris) else NA_real_
      
      # Aggregate per (k, res)
      agg <- aggregate(
        . ~ kNN_k + resolution,
        data = metrics_df[, c("kNN_k","resolution","clusters","modularity","mean_silhouette","pct_singletons")],
        FUN = function(x) mean(x, na.rm = TRUE),
        na.action = na.pass
      )
      agg$repeats <- length(seeds)
      agg$stability_ARI <- mean_ari
      
      key <- sprintf("k=%d|res=%.4f", k, r)
      per_setting[[key]] <- list(labels = run_labels, metrics = metrics_df)
      all_rows[[length(all_rows) + 1L]] <- agg
    }
  }
  
  summary_df <- do.call(rbind, all_rows)
  rownames(summary_df) <- NULL

  # --- Filter to target cluster range [min_clusters, max_clusters] ---
  in_range <- summary_df$clusters >= min_clusters & summary_df$clusters <= max_clusters
  if (sum(in_range) == 0L) {
    # Fallback: pick solutions closest to the target range
    dist_to_range <- pmin(abs(summary_df$clusters - min_clusters),
                          abs(summary_df$clusters - max_clusters))
    closest <- min(dist_to_range, na.rm = TRUE)
    in_range <- dist_to_range <= closest + 1L
    warning(sprintf(
      "No solutions in [%d, %d] clusters. Relaxing to closest (%d-%d clusters).",
      as.integer(min_clusters), as.integer(max_clusters),
      as.integer(min(summary_df$clusters[in_range])),
      as.integer(max(summary_df$clusters[in_range]))
    ))
  }
  scored_df <- summary_df[in_range, , drop = FALSE]

  # Composite score on filtered solutions (no cluster-count penalty needed)
  norm <- function(x) if (all(is.na(x))) rep(NA_real_, length(x)) else (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE) + 1e-12)

  w <- score_weights %||% list(stability = 0.25, modularity = 0.25,
                                silhouette = 0.30, singletons = 0.20)

  s_ari  <- norm(scored_df$stability_ARI)
  s_mod  <- norm(scored_df$modularity)
  s_sil  <- norm(scored_df$mean_silhouette)
  s_sing <- 1 - norm(scored_df$pct_singletons)

  scored_df$score <- w$stability * s_ari + w$modularity * s_mod +
                     w$silhouette * s_sil + w$singletons * s_sing

  # Propagate scores back to full summary (unscored rows get NA)
  summary_df$score <- NA_real_
  summary_df$score[in_range] <- scored_df$score

  best_row <- which.max(scored_df$score)
  best_k   <- scored_df$kNN_k[best_row]
  best_res <- scored_df$resolution[best_row]

  log_message(sprintf("Best Leiden: k=%d, res=%.2f, %.0f clusters (score=%.3f)",
                      as.integer(best_k), best_res, scored_df$clusters[best_row],
                      scored_df$score[best_row]))
  
  # Recompute final graph and a single labeling at best settings (seed=1)
  best_k <- min(best_k, n - 1L)
  g_best <- .make_snn_graph(emb, k = best_k, approximate = TRUE, metric = metric_knn)
  attr(g_best, "knn_k") <- best_k
  memb_best <- leidenAlg::leiden.community(g_best, resolution = best_res)$membership
  
  list(
    summary         = summary_df[order(summary_df$kNN_k, summary_df$resolution), ],
    per_setting     = per_setting,
    best = list(k = best_k, resolution = best_res, labels = memb_best, graph = g_best)
  )
}

# -------- UMAP tuning given fixed labels (optimize separation on UMAP) -----
tune_umap <- function(
    emb,
    labels,
    nn_grid      = c(10L, 15L, 30L, 50L),
    mindist_grid = c(0.05, 0.1, 0.2, 0.4),
    metric       = "euclidean",
    n_threads    = max(1L, parallel::detectCores() - 1L),
    seed         = 1L
) {
  rows <- list()
  umaps <- list()
  idx <- 1L
  for (nn in nn_grid) {
    for (md in mindist_grid) {
      set.seed(seed)
      U <- uwot::umap(
        emb,
        n_neighbors = nn,
        min_dist    = md,
        metric      = metric,
        n_threads   = n_threads,
        fast_sgd    = TRUE,
        ret_model   = FALSE,
        verbose     = FALSE
      )
      sil <- .silhouette_mean(labels, U, metric = "euclidean", max_n = 10000L, seed = seed)
      rows[[idx]] <- data.frame(n_neighbors = nn, min_dist = md, umap_silhouette = sil)
      umaps[[idx]] <- U
      idx <- idx + 1L
    }
  }
  df <- do.call(rbind, rows)
  best_idx <- which.max(df$umap_silhouette)
  if (length(best_idx) == 0L) {
    warning("All UMAP silhouette scores were NA; falling back to first parameter set")
    best_idx <- 1L
  }
  list(
    grid   = df[order(df$n_neighbors, df$min_dist), ],
    best   = list(n_neighbors = df$n_neighbors[best_idx], min_dist = df$min_dist[best_idx]),
    coords = umaps[[best_idx]]
  )
}

# -------- Post-clustering lineage validation: reassign cross-lineage cells ----
# Detects clusters contaminated with cells from a DIFFERENT major lineage
# (T, NK, B, Myeloid) and moves those minority cells to the nearest pure
# cluster of their own lineage. Cluster count is preserved (no splitting).
#
# Previously this only split Lymphoid vs Myeloid — a T+B or NK+B mix was
# missed because all three are "Lymphoid". Now every lineage is a separate
# bucket so any cross-lineage contamination is flagged.
validate_clusters_by_lineage <- function(
    obj,
    emb,
    labels,
    t_genes       = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
    nk_genes      = c("KLRF1", "SH2D1B", "NCR3", "SPON2"),
    b_genes       = c("CD19", "MS4A1", "CD79A", "CD79B", "BANK1", "MZB1", "JCHAIN"),
    myeloid_genes = c("CD14", "LYZ", "S100A8", "S100A9", "FCGR3A"),
    mix_threshold = 0.15
) {

  # 1. Ensure genes exist in the object to prevent AddModuleScore errors
  counts <- Seurat::GetAssayData(obj, layer = "data")
  .avail <- function(genes) genes[genes %in% rownames(counts)]

  features_list <- list(
    T_cell  = .avail(t_genes),
    NK_cell = .avail(nk_genes),
    B_cell  = .avail(b_genes),
    Myeloid = .avail(myeloid_genes)
  )
  # Require at least two non-empty lineage marker sets to attempt validation.
  if (sum(vapply(features_list, length, integer(1)) > 0L) < 2L) {
    message("Lineage validation skipped: insufficient marker genes found in assay.")
    return(labels)
  }

  # 2. Calculate Module Scores (background-subtracted expression)
  obj <- Seurat::AddModuleScore(
    object   = obj,
    features = features_list,
    name     = "LineageScore_",
    nbin     = 24,
    ctrl     = 100
  )

  # 3. Extract scores
  scores <- obj@meta.data[, paste0("LineageScore_", seq_along(features_list))]
  colnames(scores) <- names(features_list)

  # 4. Determine per-cell lineage call (four buckets, not two)
  max_score   <- apply(scores, 1, max)
  max_lineage <- apply(scores, 1, function(x) names(which.max(x)))
  cell_calls  <- rep("Unknown", nrow(scores))
  cell_calls[max_score > 0] <- max_lineage[max_score > 0]
  # cell_calls ∈ {T_cell, NK_cell, B_cell, Myeloid, Unknown}

  lineage_names <- names(features_list)  # target lineages

  # 5. Per-cluster lineage composition. A cluster is "mixed" if two or more
  #    lineages each exceed mix_threshold of its cells. Pure clusters carry
  #    the name of their dominant lineage.
  clusters <- sort(unique(labels))
  cluster_info <- lapply(clusters, function(cl) {
    idx <- which(labels == cl)
    fracs <- vapply(lineage_names, function(ln)
      sum(cell_calls[idx] == ln) / length(idx), numeric(1))
    names(fracs) <- lineage_names
    list(fracs = fracs,
         dominant = lineage_names[which.max(fracs)],
         n_above = sum(fracs >= mix_threshold))
  })
  names(cluster_info) <- as.character(clusters)

  is_mixed    <- vapply(cluster_info, function(x) x$n_above >= 2L, logical(1))
  dominant_of <- vapply(cluster_info, function(x) x$dominant, character(1))

  # Pure clusters per lineage (reassignment targets)
  pure_by_lineage <- split(clusters[!is_mixed], dominant_of[!is_mixed])

  if (sum(!is_mixed) == 0L) {
    message("Lineage validation: no pure clusters; skipping reassignment.")
    return(labels)
  }

  # 6. Cluster centroids in the MNN embedding space
  centroids <- t(sapply(clusters, function(cl)
    colMeans(emb[labels == cl, , drop = FALSE])))
  rownames(centroids) <- as.character(clusters)

  # 7. For each mixed cluster, reassign minority-lineage cells to the nearest
  #    pure cluster of their own lineage.
  new_labels <- labels
  n_moved <- 0L
  mixed_cls <- clusters[is_mixed]
  for (cl in mixed_cls) {
    idx      <- which(labels == cl)
    dom_lin  <- cluster_info[[as.character(cl)]]$dominant
    fracs    <- cluster_info[[as.character(cl)]]$fracs
    minority_lineages <- setdiff(lineage_names, dom_lin)

    frac_str <- paste(sprintf("%s=%.0f%%", lineage_names, fracs * 100),
                      collapse = ", ")
    message(sprintf("Cluster %s: mixed (%s, n=%d, dominant=%s)",
                    cl, frac_str, length(idx), dom_lin))

    for (lin in minority_lineages) {
      # Only move cells whose lineage has pure-cluster targets and whose
      # fraction exceeds the mix threshold
      if (fracs[lin] < mix_threshold) next
      target_cls <- pure_by_lineage[[lin]]
      if (is.null(target_cls) || length(target_cls) == 0L) next

      minority_idx <- idx[cell_calls[idx] == lin]
      if (length(minority_idx) == 0L) next

      target_centroids <- centroids[as.character(target_cls), , drop = FALSE]
      for (ci in minority_idx) {
        dists <- apply(target_centroids, 1,
                       function(ctr) sqrt(sum((emb[ci, ] - ctr)^2)))
        best_idx <- which.min(dists)
        if (length(best_idx) == 0L) next
        new_labels[ci] <- target_cls[best_idx]
        n_moved <- n_moved + 1L
      }
    }
  }

  if (n_moved > 0L) {
    log_message(sprintf(
      "Reassigned %d minority-lineage cells from %d mixed clusters",
      n_moved, length(mixed_cls)))
  } else {
    log_message("No cross-lineage mixed clusters detected")
  }

  new_labels
}

# -------- (removed: reassign_misplaced_cells — logic folded into validate_clusters_by_lineage)
