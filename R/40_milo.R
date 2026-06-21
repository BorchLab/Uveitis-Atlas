# R/40_milo.R
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(miloR)
  library(scater)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# T2.4 helper. For each nhood in `res`, count Subjects contributing ≥1 cell
# within each contrast level, then flag nhoods where any arm has fewer than
# `min_subjects`. Adds n_subjects_<level> columns and an `underpowered`
# boolean. miloR stores per-sample counts in `nhoodCounts(milo)` — we use
# that matrix instead of recomputing cell→nhood membership.
.annotate_milo_subject_n <- function(res, milo, da_meta, contrast_col,
                                     sample_col, min_subjects = 4L) {
  if (!("Nhood" %in% colnames(res))) {
    log_message("  Milo N gate: Nhood column missing; skipping.")
    return(res)
  }
  nc <- tryCatch(miloR::nhoodCounts(milo), error = function(e) NULL)
  if (is.null(nc) || nrow(nc) == 0L) {
    log_message("  Milo N gate: nhoodCounts unavailable; skipping.")
    return(res)
  }
  has_subject <- "Subject" %in% colnames(da_meta)
  subject_for <- if (has_subject)
                   stats::setNames(as.character(da_meta$Subject),
                                   as.character(da_meta[[sample_col]]))
                 else stats::setNames(colnames(nc), colnames(nc))
  group_for   <- stats::setNames(as.character(da_meta[[contrast_col]]),
                                 as.character(da_meta[[sample_col]]))
  levels_ct   <- sort(unique(stats::na.omit(group_for)))

  count_subjects <- function(nhood_row, level) {
    samples_in <- colnames(nc)[as.numeric(nhood_row) > 0]
    samples_in <- samples_in[group_for[samples_in] == level]
    length(unique(stats::na.omit(subject_for[samples_in])))
  }

  ns_mat <- matrix(0L, nrow = nrow(nc), ncol = length(levels_ct),
                   dimnames = list(NULL, paste0("n_subjects_", levels_ct)))
  for (i in seq_len(nrow(nc))) {
    row_vec <- as.numeric(nc[i, ])
    for (j in seq_along(levels_ct)) {
      ns_mat[i, j] <- count_subjects(row_vec, levels_ct[j])
    }
  }
  ns_df <- as.data.frame(ns_mat)
  ns_df$Nhood <- seq_len(nrow(nc))
  res <- dplyr::left_join(res, ns_df, by = "Nhood")
  res$underpowered <- apply(ns_mat, 1, function(r) any(r < min_subjects))[res$Nhood]
  n_under <- sum(res$underpowered, na.rm = TRUE)
  log_message(sprintf("  Milo N gate: %d / %d neighbourhoods underpowered (any arm <%d subjects)",
                      n_under, nrow(res), min_subjects))
  res
}

# Refresh a Milo's colData and reducedDims from the current Seurat object.
# - Adds/overwrites celltype_broad, merged.celltype.cluster (and anything else
#   in `meta_cols`) in colData, aligned by cell barcode.
# - Copies the Seurat "UMAP" reduction into reducedDim(milo, "UMAP") so that
#   plotNhoodGraphDA(..., layout = "UMAP") lays nhoods out on the same
#   coordinates as the Seurat UMAPs elsewhere in the pipeline.
# Safe to call on both freshly built and cached Milo objects.
sync_milo_from_seurat <- function(milo, obj,
                                  meta_cols = c("celltype_broad",
                                                "merged.celltype.cluster",
                                                "knn.leiden.cluster",
                                                "Tissue_1", "Phenotype_2",
                                                "Etiology", "Subject",
                                                "Subject_Timepoint"),
                                  umap_name = "UMAP") {
  milo_cells <- colnames(milo)
  obj_cells  <- colnames(obj)
  common     <- intersect(milo_cells, obj_cells)
  if (length(common) < length(milo_cells)) {
    log_message("  sync_milo_from_seurat: ", length(milo_cells) - length(common),
                " milo cells missing from Seurat object; those rows stay NA.")
  }

  cd <- as.data.frame(colData(milo))
  meta <- obj@meta.data
  for (mc in meta_cols) {
    if (!mc %in% colnames(meta)) next
    # Create/overwrite column, aligned by barcode
    vals <- rep(NA, length(milo_cells))
    names(vals) <- milo_cells
    vals[common] <- as.character(meta[common, mc])
    # Preserve factor levels if the Seurat col was a factor
    if (is.factor(meta[[mc]])) {
      cd[[mc]] <- factor(vals, levels = levels(meta[[mc]]))
    } else {
      cd[[mc]] <- vals
    }
  }
  # Write back, preserving rownames
  rownames(cd) <- milo_cells
  colData(milo) <- S4Vectors::DataFrame(cd, check.names = FALSE)

  # Copy Seurat UMAP → reducedDim(milo, "UMAP"), aligned by barcode.
  if (umap_name %in% names(obj@reductions)) {
    emb <- obj@reductions[[umap_name]]@cell.embeddings
    umap_mat <- matrix(NA_real_, nrow = length(milo_cells), ncol = ncol(emb),
                       dimnames = list(milo_cells, colnames(emb)))
    keep <- intersect(milo_cells, rownames(emb))
    umap_mat[keep, ] <- emb[keep, , drop = FALSE]
    reducedDim(milo, "UMAP") <- umap_mat
  } else {
    log_message("  sync_milo_from_seurat: no '", umap_name,
                "' reduction in Seurat object; nhood graph will fall back.")
  }

  milo
}

run_milo_da <- function(cfg, target = c("all", "eye", "myeloid", "bcell", "tcell")) {
  target <- match.arg(target)
  paths  <- get_target_paths(cfg, target)
  log_message("Starting MiloR differential abundance analysis (target=", target, ")...")

  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Integrated object not found. Skipping MiloR.")
    return(invisible(TRUE))
  }

  obj <- readRDS(obj_path)
  ensure_dir(paths$results_tables)
  ensure_dir(paths$viz_dir)

  # --- Build or load Milo object ---
  milo_path <- file.path(paths$results_objects, "MiloObject.rds")
  if (file.exists(milo_path)) {
    log_message("Loading existing MiloObject...")
    milo <- readRDS(milo_path)
  } else {
    log_message("Building Milo object from scratch...")
    sce <- as.SingleCellExperiment(obj)
    milo <- Milo(sce)
    reducedDim(milo, "PCA") <- obj@reductions$fastMNN@cell.embeddings

    milo <- buildGraph(milo, k = cfg$milo$k, d = cfg$milo$d, reduced.dim = "PCA")
    milo <- makeNhoods(milo, prop = cfg$milo$prop, k = cfg$milo$k,
                       d = cfg$milo$d, refined = TRUE, reduced_dims = "PCA")
    saveRDS(milo, milo_path)
    log_message("Saved MiloObject.rds")
  }

  # Always sync colData + UMAP from the CURRENT Seurat object. This ensures
  # cached MiloObjects pick up celltype_broad / merged.celltype.cluster (added
  # downstream of the original build) and carry the Seurat UMAP for plotting.
  milo <- sync_milo_from_seurat(milo, obj)

  if (target == "all") {
    # --- Cross-tissue DA: Eye vs Blood ---
    log_message("MiloR DA (target=all): Tissue Eye vs Blood")
    tryCatch({
      run_milo_contrast(
        milo, obj,
        contrast_col  = "Tissue_1",
        sample_col    = "Subject_Timepoint",
        covariates    = NULL,
        contrast_name = "tissue",
        cfg           = cfg,
        paths         = paths,
        target        = target
      )
    }, error = function(e) {
      log_message("MiloR tissue contrast failed: ", conditionMessage(e))
    })
  } else {
    # --- Within-eye / within-compartment DA across phenotype contrasts ---
    # Compartment targets share the same etiology contrast as the eye branch.
    # Granulomatous contrast is eye-only (compartments may have too few cells
    # split across NIU/Viral/Gran/Non_Gran for stable Milo testing).
    log_message("MiloR DA (target=", target, "): phenotype contrasts")
    eye_contrasts <- if (target == "eye") {
      list(
        list(name = "etiology",  col = "Phenotype_2", groups = c("NIU",  "Viral")),
        list(name = "granulom",  col = "Phenotype_2", groups = c("Gran", "Non_Gran"))
      )
    } else {
      list(
        list(name = paste0("etiology_", target), col = "Phenotype_2", groups = c("NIU", "Viral"))
      )
    }
    for (ec in eye_contrasts) {
      tryCatch({
        cells_keep <- colnames(obj)[obj[[ec$col, drop = TRUE]] %in% ec$groups]
        if (length(cells_keep) < 100) {
          log_message("  Too few cells for ", ec$name, ", skipping.")
          next
        }
        obj_sub <- subset(obj, cells = cells_keep)
        sce_sub <- as.SingleCellExperiment(obj_sub)
        milo_sub <- Milo(sce_sub)
        reducedDim(milo_sub, "PCA") <- obj_sub@reductions$fastMNN@cell.embeddings
        milo_sub <- buildGraph(milo_sub, k = cfg$milo$k, d = cfg$milo$d,
                               reduced.dim = "PCA")
        milo_sub <- makeNhoods(milo_sub, prop = cfg$milo$prop, k = cfg$milo$k,
                               d = cfg$milo$d, refined = TRUE, reduced_dims = "PCA")
        milo_sub <- sync_milo_from_seurat(milo_sub, obj_sub)
        run_milo_contrast(
          milo_sub, obj_sub,
          contrast_col  = ec$col,
          sample_col    = "Subject_Timepoint",
          covariates    = NULL,
          contrast_name = paste0(ec$name, "_eye"),
          cfg           = cfg,
          paths         = paths,
          target        = target
        )
      }, error = function(e) {
        log_message("MiloR ", ec$name, " contrast failed: ", conditionMessage(e))
      })
    }
  }

  rm(milo, obj)
  invisible(TRUE)
}

# --- Helper: run a single MiloR contrast ---
# try_subject_block: when TRUE (default), auto-add `+ Subject` to the design
# if the contrast varies within Subject (paired eye/blood) and the resulting
# design is rank-sufficient. Set FALSE to force an unpaired design.
run_milo_contrast <- function(milo, obj, contrast_col, sample_col,
                               covariates = NULL, contrast_name, cfg,
                               paths = NULL,
                               try_subject_block = TRUE,
                               target = NULL) {
  if (is.null(paths)) paths <- cfg$paths
  # Resolve sample_col — prefer provided, then fall back
  if (!sample_col %in% colnames(colData(milo))) {
    if ("Subject" %in% colnames(colData(milo))) {
      sample_col <- "Subject"
    } else {
      sample_col <- "orig.ident"
    }
  }

  # Pre-flight: sample_col must uniquely identify each biological sample
  # *within* the contrast. If one sample_col value (e.g. Subject_Timepoint)
  # spans multiple levels of contrast_col (paired Eye/Blood on the same
  # subject+timepoint), countCells builds a malformed sparse matrix and
  # dies with "too many replacement values" out of replCmat4. Fall back to
  # orig.ident preemptively — Subject pairing is still carried through the
  # design matrix via the blocking-covariate logic below.
  cd_pre <- as.data.frame(colData(milo))
  if (all(c(sample_col, contrast_col) %in% colnames(cd_pre))) {
    pre_keys <- cd_pre[, c(sample_col, contrast_col), drop = FALSE]
    pre_keys <- pre_keys[stats::complete.cases(pre_keys), , drop = FALSE]
    sample_by_contrast <- unique(pre_keys)
    if (anyDuplicated(sample_by_contrast[[sample_col]])) {
      log_message("  sample_col '", sample_col, "' not 1:1 with '",
                  contrast_col, "'; falling back to orig.ident preemptively.")
      sample_col <- "orig.ident"
    }
  }

  # Count cells per neighbourhood per sample
  milo <- countCells(milo, meta.data = as.data.frame(colData(milo)),
                     samples = sample_col)

  # Build sample-level metadata; include Subject if available (for blocking).
  # Drop cells with NA in any of the key columns — these can't contribute
  # to the contrast and would otherwise split a single sample_col value
  # into multiple da_meta rows (e.g. cells lingering in the Milo graph
  # after being removed from the current Seurat object get NA from
  # sync_milo_from_seurat and create spurious (sample, NA) rows alongside
  # the real (sample, Eye) / (sample, Blood) rows).
  base_cols <- c(sample_col, contrast_col)
  cd_all <- as.data.frame(colData(milo))
  has_subject <- "Subject" %in% colnames(cd_all)
  if (has_subject) base_cols <- unique(c(base_cols, "Subject"))

  cd_nonmiss <- cd_all[stats::complete.cases(cd_all[, base_cols, drop = FALSE]), ,
                       drop = FALSE]
  da_meta <- unique(cd_nonmiss[, base_cols, drop = FALSE])

  # Second-chance duplicate check.
  # After NA filtering the pre-flight should have already resolved the
  # tissue-style ambiguity, so this branch only fires in edge cases
  # (e.g. sample_col was never the pre-flight candidate). Short-circuit
  # if we're already on orig.ident so we don't re-run countCells pointlessly.
  if (anyDuplicated(da_meta[[sample_col]])) {
    if (sample_col != "orig.ident") {
      log_message("  sample_col '", sample_col,
                  "' has duplicates; falling back to orig.ident")
      sample_col <- "orig.ident"
      milo <- countCells(milo, meta.data = as.data.frame(colData(milo)),
                         samples = sample_col)
      base_cols <- unique(c(sample_col, contrast_col,
                            if (has_subject) "Subject" else character(0)))
      cd_all <- as.data.frame(colData(milo))
      cd_nonmiss <- cd_all[stats::complete.cases(cd_all[, base_cols, drop = FALSE]), ,
                           drop = FALSE]
      da_meta <- unique(cd_nonmiss[, base_cols, drop = FALSE])
    }
    if (anyDuplicated(da_meta[[sample_col]])) {
      stop("sample_col '", sample_col, "' still has duplicate values after ",
           "NA filtering and orig.ident fallback. This means one sample is ",
           "genuinely aliased across contrast or Subject levels — check the ",
           "colData for cells with inconsistent annotations.")
    }
  }

  rownames(da_meta) <- da_meta[[sample_col]]
  da_meta[[contrast_col]] <- factor(da_meta[[contrast_col]])

  if (length(levels(da_meta[[contrast_col]])) < 2) {
    log_message("Only one level in ", contrast_col, ". Skipping.")
    return(invisible(NULL))
  }

  # Auto-detect whether to block by Subject (paired design).
  # Condition: Subject must vary within the contrast AND the augmented
  # design must be rank-sufficient (no aliasing).
  if (try_subject_block && has_subject && !"Subject" %in% covariates) {
    per_sbj_grp <- tapply(as.character(da_meta[[contrast_col]]),
                          as.character(da_meta$Subject),
                          function(x) length(unique(x)))
    within_variation <- any(per_sbj_grp > 1, na.rm = TRUE)
    if (within_variation) {
      test_f <- as.formula(paste("~ Subject +", contrast_col))
      da_meta$Subject <- factor(as.character(da_meta$Subject))
      mm_try <- try(model.matrix(test_f, data = da_meta), silent = TRUE)
      if (!inherits(mm_try, "try-error") &&
          qr(mm_try)$rank == ncol(mm_try) &&
          nrow(mm_try) > ncol(mm_try)) {
        covariates <- unique(c(covariates, "Subject"))
        log_message("  Paired design: adding Subject as blocking covariate.")
      } else {
        log_message("  Subject blocking skipped (rank-deficient design).")
      }
    }
  }

  formula_str <- paste("~", paste(c(covariates, contrast_col), collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = da_meta)

  # Calculate nhood distances for spatial FDR
  milo <- calcNhoodDistance(milo, d = cfg$milo$d, reduced.dim = "PCA")

  # Test
  res <- testNhoods(milo, design = design, design.df = da_meta,
                    reduced.dim = "PCA")

  # Annotate nhoods with all grouping variables present in colData.
  # knn.leiden.cluster is the fine-grained Leiden labelling; for the eye
  # target it is the only cluster column (no merge_clusters step).
  annot_cols <- c("celltype_broad", "merged.celltype.cluster",
                  "knn.leiden.cluster")
  for (ac in annot_cols) {
    if (ac %in% colnames(colData(milo))) {
      res <- annotateNhoods(milo, res, coldata_col = ac)
    } else {
      log_message("Column '", ac, "' not found in milo colData; skipping annotation.")
    }
  }

  # T2.4: pseudobulk-N gate for compartment targets. For each nhood, count
  # the number of distinct Subjects contributing >=1 cell within each
  # contrast level, and flag nhoods where any arm has <min_subjects.
  # Eye / full-object milo outputs are unchanged.
  if (!is.null(target) && target %in% c("myeloid", "bcell", "tcell")) {
    res <- .annotate_milo_subject_n(res, milo, da_meta, contrast_col,
                                    sample_col, min_subjects = 4L)
  }

  # Save full nhood-level results
  write.csv(res,
            file.path(paths$results_tables, paste0("Milo_", contrast_name, "_DA.csv")),
            row.names = FALSE)
  log_message("Saved MiloR results: Milo_", contrast_name, "_DA.csv")

  # Per-group summary tables (one per annotation)
  summarise_by <- function(res_df, group_col, purity_col = NULL,
                           purity_min = 0.7, fdr_thresh = 0.1) {
    if (!group_col %in% colnames(res_df)) return(NULL)
    df <- res_df %>% filter(!is.na(.data[[group_col]]))
    if (!is.null(purity_col) && purity_col %in% colnames(df)) {
      df <- df %>% filter(.data[[purity_col]] >= purity_min)
    }
    if (nrow(df) == 0) return(NULL)
    df %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        n_neighborhoods = n(),
        n_sig_DA     = sum(SpatialFDR < fdr_thresh, na.rm = TRUE),
        n_enriched   = sum(SpatialFDR < fdr_thresh & logFC > 0, na.rm = TRUE),
        n_depleted   = sum(SpatialFDR < fdr_thresh & logFC < 0, na.rm = TRUE),
        mean_logFC   = mean(logFC, na.rm = TRUE),
        median_logFC = median(logFC, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(contrast = contrast_name, grouping = group_col)
  }

  # celltype_broad summary
  s_broad <- summarise_by(res, "celltype_broad",
                          purity_col = "celltype_broad_fraction")
  if (!is.null(s_broad)) {
    write.csv(s_broad,
              file.path(paths$results_tables,
                        paste0("milo_summary_", contrast_name, "_celltype_broad.csv")),
              row.names = FALSE)
    log_message("Saved: milo_summary_", contrast_name, "_celltype_broad.csv")
  }

  # merged.celltype.cluster summary
  s_cluster <- summarise_by(res, "merged.celltype.cluster",
                            purity_col = "merged.celltype.cluster_fraction")
  if (!is.null(s_cluster)) {
    write.csv(s_cluster,
              file.path(paths$results_tables,
                        paste0("milo_summary_", contrast_name, "_merged_cluster.csv")),
              row.names = FALSE)
    log_message("Saved: milo_summary_", contrast_name, "_merged_cluster.csv")
  }

  # knn.leiden.cluster summary
  s_leiden <- summarise_by(res, "knn.leiden.cluster",
                           purity_col = "knn.leiden.cluster_fraction")
  if (!is.null(s_leiden)) {
    write.csv(s_leiden,
              file.path(paths$results_tables,
                        paste0("milo_summary_", contrast_name, "_knn_leiden.csv")),
              row.names = FALSE)
    log_message("Saved: milo_summary_", contrast_name, "_knn_leiden.csv")
  }

  # Build nhood graph for plotting
  milo <- buildNhoodGraph(milo)

  # Persist the per-contrast Milo object (post-test, graph included) so
  # downstream viz (plotNhoodGraphDA) can reuse it without rebuilding.
  milo_contrast_path <- file.path(paths$results_objects,
                                  paste0("MiloObject_", contrast_name, ".rds"))
  tryCatch({
    saveRDS(milo, milo_contrast_path)
    log_message("Saved: MiloObject_", contrast_name, ".rds")
  }, error = function(e) {
    log_message("WARN: failed to save MiloObject_", contrast_name,
                ".rds: ", conditionMessage(e))
  })

  # Helper: strip edge layers from miloR nhood graphs so only the nhood
  # points (vertices) remain. Edges between overlapping neighbourhoods are
  # visually noisy and obscure logFC on dense graphs.
  drop_edge_layers <- function(p) {
    keep <- vapply(p$layers, function(l) {
      !any(grepl("edge", class(l$geom), ignore.case = TRUE))
    }, logical(1))
    p$layers <- p$layers[keep]
    p
  }

  # Nhood graph DA plot — viridis fill/color, edges suppressed.
  tryCatch({
    p <- plotNhoodGraphDA(milo, res, alpha = 0.05)
    p <- drop_edge_layers(p)
    p <- p +
      viridis::scale_fill_viridis(option = "viridis",
                                  na.value = "grey85", name = "logFC") +
      viridis::scale_color_viridis(option = "viridis",
                                   na.value = "grey85", name = "logFC") +
      ggtitle(paste("MiloR DA:", contrast_name))
    ggsave(file.path(paths$viz_dir, paste0("milo_", contrast_name, "_nhood.pdf")),
           p, width = 8, height = 6, dpi = 300)
  }, error = function(e) {
    log_message("MiloR nhood plot failed: ", conditionMessage(e))
  })

  # Beeswarm plots for each annotation (logFC distribution per group). Color
  # by logFC on a viridis scale so direction is encoded in the same palette
  # as the nhood graph.
  for (ac in c("celltype_broad", "merged.celltype.cluster",
               "knn.leiden.cluster")) {
    if (!ac %in% colnames(res)) next
    tryCatch({
      pb <- plotDAbeeswarm(res, group.by = ac, alpha = 0.1)
      pb <- pb +
        viridis::scale_color_viridis(option = "viridis", name = "logFC") +
        ggtitle(paste("MiloR DA:", contrast_name, "by", ac))
      ggsave(file.path(paths$viz_dir,
                       paste0("milo_", contrast_name, "_beeswarm_", ac, ".pdf")),
             pb, width = 8, height = 6, dpi = 300)
    }, error = function(e) {
      log_message("MiloR beeswarm (", ac, ") failed: ", conditionMessage(e))
    })
  }

  invisible(res)
}
