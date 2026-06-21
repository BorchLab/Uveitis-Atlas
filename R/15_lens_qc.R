# R/15_lens_qc.R
#
# Pre-compartment lens contamination filter (UCell-based).
#
# Why this exists: the eye-stage lens filter in R/20_subset_eye.R uses
# sum-of-log-norm expression on a fixed gene panel. That catches the heaviest
# contamination before re-integration but leaves a residual tail of lens-
# adjacent cells (most visibly in cluster 11, where mregDC and lens fibers
# co-cluster after re-integration). Rather than reinstate the cluster-11
# hygiene split, we score every eye cell with UCell on the same lens gene
# panel and apply a per-cell threshold before the T/B/Myeloid compartment
# split.
#
# UCell is rank-based, so the score is robust to per-cell sequencing depth
# differences that the sum-of-log-norm filter is sensitive to. Same package
# already used by R/32_escape.R (cfg$escape$method = "UCell"), so no new
# install requirement.
#
# Two entry points:
#   * run_lens_ucell_diagnostic(cfg) — compute scores on the eye object,
#     write per-cell distribution + per-cluster summary + suggested
#     threshold via density-valley detection. Re-runnable at zero cost
#     while tuning. Does NOT modify the eye object.
#   * .apply_lens_ucell_filter(eye_obj, cfg) — helper called from
#     subset_compartments at the start of Phase 1c. Reads
#     cfg$compartments$lens_ucell_filter$threshold (must be set after the
#     diagnostic), drops cells above the threshold, writes the
#     lens_ucell_filter_report.csv audit trail, and returns the filtered
#     Seurat object.
#
# Outputs:
#   outputs/tables/eye/lens_ucell_distribution.csv
#   outputs/tables/eye/lens_ucell_by_cluster.csv
#   outputs/tables/eye/lens_ucell_threshold_suggestion.csv
#   outputs/tables/eye/lens_ucell_filter_report.csv  (written at apply time)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Resolve the lens gene panel. Prefer cfg$compartments$lens_ucell_filter$genes
# when set; fall back to cfg$eye_qc$lens_filter$genes so the two filters
# share a panel by default and edits only need to land in one place.
.lens_ucell_genes <- function(cfg) {
  g <- cfg$compartments$lens_ucell_filter$genes %||%
       cfg$eye_qc$lens_filter$genes %||%
       character(0)
  unique(as.character(g))
}

# Run UCell on the lens panel. Returns a numeric vector named by cell barcode.
.compute_lens_ucell_scores <- function(eye_obj, lens_genes, seed = 42L) {
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("UCell is required for the lens UCell filter. Install with ",
         "BiocManager::install('UCell').")
  }
  present <- intersect(lens_genes, rownames(eye_obj))
  if (length(present) < 3L) {
    stop("Lens UCell filter: only ", length(present),
         " of ", length(lens_genes), " genes present in the eye object. ",
         "Need at least 3.")
  }
  log_message(sprintf("Lens UCell: %d / %d genes present", length(present), length(lens_genes)))

  DefaultAssay(eye_obj) <- "RNA"
  eye_obj <- JoinLayers(eye_obj)
  features <- list(Lens = present)

  set_global_seed(seed)
  scored <- UCell::AddModuleScore_UCell(eye_obj, features = features,
                                        name = "_ucell", assay = "RNA")
  scores <- scored[["Lens_ucell", drop = TRUE]]
  names(scores) <- colnames(scored)
  scores
}

# Suggest a threshold from one score distribution. Tries density-valley
# bimodal detection first (the lens contamination tail typically sits as a
# distinct mode at the upper end of the score range); falls back to a high
# percentile when no clean valley exists.
#
# Returns a list with $method ("density_valley" or "percentile_<q>"),
# $threshold, and $diagnostic (a single-row data.frame for the suggestion
# CSV).
.suggest_lens_threshold_one <- function(scores, fallback_quantile = 0.99) {
  scores <- scores[is.finite(scores)]
  out <- list(method = NA_character_, threshold = NA_real_,
              diagnostic = data.frame())

  if (length(scores) < 100L) {
    return(out)
  }

  d <- stats::density(scores, n = 2048)
  y <- d$y; z <- d$x

  dy    <- diff(y)
  s     <- sign(dy)
  curv  <- diff(s)
  peaks <- which(curv == -2) + 1L
  pits  <- which(curv ==  2) + 1L

  peaks <- peaks[y[peaks] > 0.02 * max(y)]
  thr_valley <- NA_real_

  if (length(peaks) >= 2L && length(pits) > 0L) {
    top2 <- peaks[order(y[peaks], decreasing = TRUE)][seq_len(2L)]
    rng  <- sort(range(top2))
    between <- pits[pits > rng[1] & pits < rng[2]]
    if (length(between) > 0L) {
      thr_valley <- z[between[which.min(y[between])]]
    }
  }

  thr_quant <- stats::quantile(scores, fallback_quantile, names = FALSE)

  if (is.finite(thr_valley)) {
    out$method    <- "density_valley"
    out$threshold <- thr_valley
  } else {
    out$method    <- paste0("percentile_", round(100 * fallback_quantile))
    out$threshold <- thr_quant
  }

  out$diagnostic <- data.frame(
    method                = out$method,
    threshold             = round(out$threshold, 4),
    threshold_valley      = round(thr_valley, 4),
    threshold_percentile  = round(thr_quant, 4),
    fallback_quantile     = fallback_quantile,
    n_scores              = length(scores),
    n_peaks_detected      = length(peaks),
    n_pits_detected       = length(pits),
    stringsAsFactors      = FALSE
  )
  out
}

# Global + per-batch threshold suggestion. When `batches` is NULL or has a
# single level, returns only the global suggestion in $global. Otherwise also
# returns $per_batch (data.frame keyed by batch level) and $batch_effect
# (single-row data.frame with summary metrics for divergence across batches).
.suggest_lens_threshold <- function(scores, batches = NULL, fallback_quantile = 0.99) {
  out <- list(global = .suggest_lens_threshold_one(scores, fallback_quantile),
              per_batch = NULL, batch_effect = NULL)

  if (is.null(batches)) return(out)
  batches <- as.character(batches)
  if (length(batches) != length(scores))
    stop("Lens UCell: `batches` length must match `scores` length.")

  levels_present <- unique(batches[!is.na(batches)])
  if (length(levels_present) < 2L) return(out)

  per_rows <- lapply(levels_present, function(lv) {
    idx <- which(batches == lv & is.finite(scores))
    s   <- scores[idx]
    sug <- .suggest_lens_threshold_one(s, fallback_quantile)
    if (nrow(sug$diagnostic) == 0L) {
      data.frame(
        batch                = lv,
        n_cells              = length(s),
        method               = NA_character_,
        threshold            = NA_real_,
        threshold_valley     = NA_real_,
        threshold_percentile = NA_real_,
        median               = if (length(s)) round(stats::median(s), 4) else NA_real_,
        q95                  = if (length(s)) round(stats::quantile(s, 0.95, names = FALSE), 4) else NA_real_,
        q99                  = if (length(s)) round(stats::quantile(s, 0.99, names = FALSE), 4) else NA_real_,
        n_above_global       = NA_integer_,
        stringsAsFactors     = FALSE
      )
    } else {
      d <- sug$diagnostic
      data.frame(
        batch                = lv,
        n_cells              = length(s),
        method               = d$method,
        threshold            = d$threshold,
        threshold_valley     = d$threshold_valley,
        threshold_percentile = d$threshold_percentile,
        median               = round(stats::median(s), 4),
        q95                  = round(stats::quantile(s, 0.95, names = FALSE), 4),
        q99                  = round(stats::quantile(s, 0.99, names = FALSE), 4),
        n_above_global       = sum(s > out$global$threshold, na.rm = TRUE),
        stringsAsFactors     = FALSE
      )
    }
  })
  per_batch <- dplyr::bind_rows(per_rows)
  per_batch$pct_above_global <- round(100 * per_batch$n_above_global / per_batch$n_cells, 2)

  # Cross-batch divergence summary. `median_range` and `threshold_range`
  # quantify how much the per-batch distributions shift; high values relative
  # to the global IQR indicate the global threshold will over-/under-filter
  # specific batches and per-batch thresholds should be used instead.
  finite_scores <- scores[is.finite(scores)]
  global_iqr <- if (length(finite_scores) >= 4L)
                  stats::IQR(finite_scores) else NA_real_
  out$batch_effect <- data.frame(
    n_batches            = length(levels_present),
    median_min           = round(min(per_batch$median, na.rm = TRUE), 4),
    median_max           = round(max(per_batch$median, na.rm = TRUE), 4),
    median_range         = round(diff(range(per_batch$median, na.rm = TRUE)), 4),
    threshold_min        = round(min(per_batch$threshold, na.rm = TRUE), 4),
    threshold_max        = round(max(per_batch$threshold, na.rm = TRUE), 4),
    threshold_range      = round(diff(range(per_batch$threshold, na.rm = TRUE)), 4),
    global_iqr           = round(global_iqr, 4),
    median_range_vs_iqr  = round(diff(range(per_batch$median, na.rm = TRUE)) / global_iqr, 3),
    stringsAsFactors     = FALSE
  )

  out$per_batch <- per_batch
  out
}

# Diagnostic entry point — call once (or repeatedly while tuning), inspect
# outputs, then commit a threshold to cfg$compartments$lens_ucell_filter$threshold.
run_lens_ucell_diagnostic <- function(cfg) {
  paths_eye <- get_target_paths(cfg, "eye")
  eye_path  <- file.path(paths_eye$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(eye_path)) {
    log_message("Lens UCell diagnostic: eye object not found at ", eye_path, ". Skipping.")
    return(invisible(FALSE))
  }
  log_message("Lens UCell diagnostic: loading eye object")
  eye_obj <- readRDS(eye_path)

  lens_genes <- .lens_ucell_genes(cfg)
  if (length(lens_genes) == 0L) {
    log_message("Lens UCell diagnostic: no genes configured; skipping.")
    return(invisible(FALSE))
  }

  scores <- .compute_lens_ucell_scores(eye_obj, lens_genes,
                                       seed = cfg$seed %||% 42L)

  meta <- eye_obj[[]]
  cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                   "knn.leiden.cluster" else "seurat_clusters"

  batch_key <- cfg$compartments$lens_ucell_filter$batch_key
  batch_vec <- NULL
  if (!is.null(batch_key)) {
    if (!(batch_key %in% colnames(meta))) {
      log_message(sprintf(
        "Lens UCell: batch_key '%s' not in metadata; falling back to global analysis.",
        batch_key))
      batch_key <- NULL
    } else {
      batch_vec <- as.character(meta[[batch_key]])
    }
  }

  per_cell <- data.frame(
    cell_id         = colnames(eye_obj),
    lens_ucell      = scores,
    eye_cluster     = as.character(meta[[cluster_col]]),
    celltype        = { c <- resolve_celltype(meta);       if (!is.null(c)) as.character(meta[[c]]) else NA_character_ },
    celltype_broad  = { c <- resolve_celltype_broad(meta); if (!is.null(c)) as.character(meta[[c]]) else NA_character_ },
    Tissue_1        = if ("Tissue_1" %in% colnames(meta)) as.character(meta$Tissue_1) else NA_character_,
    batch           = if (!is.null(batch_vec)) batch_vec else NA_character_,
    stringsAsFactors = FALSE
  )

  ensure_dir(paths_eye$results_tables)
  ensure_dir(paths_eye$viz_dir)

  per_cell_path <- file.path(paths_eye$results_tables, "lens_ucell_distribution.csv")
  write.csv(per_cell, per_cell_path, row.names = FALSE)
  log_message("Wrote ", per_cell_path)

  per_cluster <- per_cell %>%
    dplyr::group_by(eye_cluster) %>%
    dplyr::summarise(
      n_cells   = dplyr::n(),
      mean      = round(mean(lens_ucell, na.rm = TRUE), 4),
      median    = round(stats::median(lens_ucell, na.rm = TRUE), 4),
      q90       = round(stats::quantile(lens_ucell, 0.90, na.rm = TRUE, names = FALSE), 4),
      q95       = round(stats::quantile(lens_ucell, 0.95, na.rm = TRUE, names = FALSE), 4),
      q99       = round(stats::quantile(lens_ucell, 0.99, na.rm = TRUE, names = FALSE), 4),
      max       = round(max(lens_ucell, na.rm = TRUE), 4),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(median))
  per_cluster_path <- file.path(paths_eye$results_tables, "lens_ucell_by_cluster.csv")
  write.csv(per_cluster, per_cluster_path, row.names = FALSE)
  log_message("Wrote ", per_cluster_path)

  fb_q <- cfg$compartments$lens_ucell_filter$fallback_quantile %||% 0.99
  sugg <- .suggest_lens_threshold(scores, batches = batch_vec,
                                  fallback_quantile = fb_q)

  global_thr <- sugg$global$threshold
  if (nrow(sugg$global$diagnostic) > 0L) {
    sugg$global$diagnostic$n_above_suggested  <- sum(scores > global_thr, na.rm = TRUE)
    sugg$global$diagnostic$pct_above_suggested <- round(
      100 * sugg$global$diagnostic$n_above_suggested / length(scores), 2)
    sugg_path <- file.path(paths_eye$results_tables, "lens_ucell_threshold_suggestion.csv")
    write.csv(sugg$global$diagnostic, sugg_path, row.names = FALSE)
    log_message(sprintf(
      "Suggested GLOBAL lens threshold (%s): %.4f — drops %d cells (%.2f%%)",
      sugg$global$method, global_thr,
      sugg$global$diagnostic$n_above_suggested,
      sugg$global$diagnostic$pct_above_suggested))
    log_message("Wrote ", sugg_path)
  }

  # Per-batch outputs when batch_key is set and at least 2 batches exist.
  if (!is.null(sugg$per_batch)) {
    pb_path <- file.path(paths_eye$results_tables,
                         "lens_ucell_threshold_suggestion_per_batch.csv")
    write.csv(sugg$per_batch, pb_path, row.names = FALSE)
    log_message("Wrote ", pb_path)

    be_path <- file.path(paths_eye$results_tables,
                         "lens_ucell_batch_effect_summary.csv")
    write.csv(sugg$batch_effect, be_path, row.names = FALSE)
    log_message("Wrote ", be_path)

    be <- sugg$batch_effect
    log_message(sprintf(
      "Batch (%s): %d levels, median range %.4f (%.2fx global IQR), threshold range %.4f",
      batch_key, be$n_batches, be$median_range,
      be$median_range_vs_iqr, be$threshold_range))
    if (isTRUE(be$median_range_vs_iqr > 0.5)) {
      log_message(
        "WARN: per-batch median range exceeds 50% of global IQR. ",
        "Per-batch thresholds recommended over a single global threshold.")
    }
  }

  invisible(list(scores = scores, suggestion = sugg, batch_key = batch_key))
}

# Resolve the threshold config into a per-cell numeric vector. Supports
# three modes:
#   * scalar (numeric)         — global threshold applied to every cell.
#   * named list (per-batch)   — keyed by levels of cfg$compartments$
#                                 lens_ucell_filter$batch_key. Cells whose
#                                 batch is not in the list keep a +Inf
#                                 threshold (i.e. never dropped) and the
#                                 caller is warned.
# Returns a numeric vector of length ncol(eye_obj).
.resolve_lens_thresholds <- function(eye_obj, flt) {
  thr <- flt$threshold
  if (is.null(thr)) {
    stop("Lens UCell filter enabled but cfg$compartments$lens_ucell_filter$threshold ",
         "is not set. Run run_lens_ucell_diagnostic(cfg) and commit a threshold ",
         "(see outputs/tables/eye/lens_ucell_threshold_suggestion.csv ",
         "or _per_batch.csv).", call. = FALSE)
  }
  n <- ncol(eye_obj)

  if (is.numeric(thr) && length(thr) == 1L && is.finite(thr)) {
    return(rep(as.numeric(thr), n))
  }

  if (is.list(thr)) {
    batch_key <- flt$batch_key
    if (is.null(batch_key))
      stop("Per-batch threshold list provided but cfg$compartments$",
           "lens_ucell_filter$batch_key is not set.", call. = FALSE)
    meta <- eye_obj[[]]
    if (!(batch_key %in% colnames(meta)))
      stop("Lens UCell filter: batch_key '", batch_key,
           "' not present in metadata.", call. = FALSE)
    batches <- as.character(meta[[batch_key]])

    missing_levels <- setdiff(unique(batches[!is.na(batches)]), names(thr))
    if (length(missing_levels) > 0L) {
      log_message(sprintf(
        "WARN: per-batch threshold missing for batches: %s. These cells are kept.",
        paste(missing_levels, collapse = ", ")))
    }

    per_cell_thr <- vapply(batches, function(b) {
      if (is.na(b)) return(Inf)
      v <- thr[[b]]
      if (is.null(v) || !is.finite(as.numeric(v))) return(Inf)
      as.numeric(v)
    }, numeric(1), USE.NAMES = FALSE)
    return(per_cell_thr)
  }

  stop("Lens UCell filter: threshold must be a scalar numeric or a named ",
       "list of per-batch thresholds; got class '", class(thr)[1], "'.",
       call. = FALSE)
}

# Apply-time helper called by subset_compartments. Reads the user-committed
# threshold from cfg$compartments$lens_ucell_filter$threshold, drops cells
# above it (with per-batch resolution when configured), writes the audit
# trail, and returns the filtered eye object.
#
# If cfg$compartments$lens_ucell_filter$enable is FALSE (or missing), this
# is a no-op and returns the input unchanged.
.apply_lens_ucell_filter <- function(eye_obj, cfg) {
  flt <- cfg$compartments$lens_ucell_filter %||% list()
  if (!isTRUE(flt$enable)) {
    log_message("Lens UCell filter disabled; passing eye object through unchanged.")
    return(eye_obj)
  }
  per_cell_thr <- .resolve_lens_thresholds(eye_obj, flt)

  lens_genes <- .lens_ucell_genes(cfg)
  scores <- .compute_lens_ucell_scores(eye_obj, lens_genes,
                                       seed = cfg$seed %||% 42L)

  keep_mask <- scores <= per_cell_thr
  keep_mask[is.na(keep_mask)] <- TRUE   # don't drop cells with missing score
  n_drop <- sum(!keep_mask)

  mode_lbl <- if (length(unique(per_cell_thr[is.finite(per_cell_thr)])) > 1L)
                "per-batch" else "global"
  log_message(sprintf(
    "Lens UCell filter (%s): dropping %d / %d cells (%.2f%%)",
    mode_lbl, n_drop, length(scores), 100 * n_drop / length(scores)))

  meta <- eye_obj[[]]
  cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                   "knn.leiden.cluster" else "seurat_clusters"
  batch_key <- flt$batch_key
  report <- data.frame(
    cell_id         = colnames(eye_obj),
    lens_ucell      = scores,
    threshold       = per_cell_thr,
    kept            = keep_mask,
    eye_cluster     = as.character(meta[[cluster_col]]),
    celltype        = { c <- resolve_celltype(meta);       if (!is.null(c)) as.character(meta[[c]]) else NA_character_ },
    celltype_broad  = { c <- resolve_celltype_broad(meta); if (!is.null(c)) as.character(meta[[c]]) else NA_character_ },
    stringsAsFactors = FALSE
  )
  if (!is.null(batch_key) && batch_key %in% colnames(meta)) {
    report$batch <- as.character(meta[[batch_key]])
  }
  paths_eye <- get_target_paths(cfg, "eye")
  ensure_dir(paths_eye$results_tables)
  report_path <- file.path(paths_eye$results_tables, "lens_ucell_filter_report.csv")
  write.csv(report, report_path, row.names = FALSE)
  log_message("Wrote ", report_path)

  if (n_drop == 0L) return(eye_obj)
  if (n_drop == ncol(eye_obj)) {
    stop("Lens UCell filter: threshold drops every cell. Threshold likely too low.",
         call. = FALSE)
  }
  subset(eye_obj, cells = colnames(eye_obj)[keep_mask])
}
