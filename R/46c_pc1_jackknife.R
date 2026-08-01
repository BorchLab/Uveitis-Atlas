# R/46c_pc1_jackknife.R
#
# Subject-level leave-one-out 
#
# Entry point: run_pc1_jackknife(cfg), gated on cfg$steps$pc1_jackknife.
#
# Outputs under outputs/tables/stress_sensitivity/jackknife/:
#   pca_pc1_loo_jackknife.csv     one row per (fold, target, substate)
#   pca_pc1_loo_summary.csv       one row per (target, substate)
#   pc1_bridge_loo_jackknife.csv  one row per (fold, weighting)
#   pc1_bridge_loo_summary.csv    one row per weighting
# Figures under outputs/viz/stress_sensitivity/.
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

.loo_cohens_d <- function(v, g) {
  x <- v[g == "NIU"]; y <- v[g == "Viral"]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  sp <- sqrt(((length(x) - 1) * stats::var(x) + (length(y) - 1) * stats::var(y)) /
               (length(x) + length(y) - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(y) - mean(x)) / sp
}

# Drop one subject's columns from a cached pseudobulk list. Substates that fall
# below 4 columns or lose a group entirely are dropped for that fold, which
# compute_per_substate_pca would do anyway.
.loo_subset_pbs <- function(pbs, drop_samples) {
  out <- lapply(pbs, function(pb) {
    if (inherits(pb, "SummarizedExperiment")) {
      k <- !(as.character(SummarizedExperiment::colData(pb)$sample) %in% drop_samples)
      if (sum(k) < 4L) return(NULL)
      pb <- pb[, k, drop = FALSE]
      if (length(unique(as.character(SummarizedExperiment::colData(pb)$group))) < 2L)
        return(NULL)
      pb
    } else {
      k <- !(as.character(pb$coldata$sample) %in% drop_samples)
      if (sum(k) < 4L) return(NULL)
      if (length(unique(as.character(pb$coldata$group[k]))) < 2L) return(NULL)
      list(counts = pb$counts[, k, drop = FALSE],
           coldata = pb$coldata[k, , drop = FALSE])
    }
  })
  Filter(Negate(is.null), out)
}

.loo_pca <- function(cfg, target, jk_dir) {
  p <- get_target_paths(cfg, target)
  op <- file.path(p$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(op)) {
    log_message("  jackknife[", target, "]: object missing; skipping.")
    return(NULL)
  }
  log_message("  jackknife[", target, "]: loading object")
  obj <- readRDS(op)
  cpcfg <- cfg$compartment_pca %||% list()
  floor <- .pca_min_cells(cfg, target)
  pca_args <- list(
    min_gene_count = as.integer(cpcfg$min_gene_count %||% 10L),
    hvg_n          = as.integer(cpcfg$hvg_n %||% 2000L),
    n_pcs          = as.integer(cpcfg$n_pcs %||% 5L),
    vst_blind      = isTRUE(cpcfg$vst_blind),
    pc1_split_fdr  = as.numeric(cpcfg$pc1_split_fdr %||% 0.05))

  log_message("  jackknife[", target, "]: building pseudobulks once (floor ",
              floor, ")")
  pbs <- build_per_substate_pseudobulks(obj, min_cells_per_pb = floor)
  if (!length(pbs)) { rm(obj); gc(verbose = FALSE); return(NULL) }

  full <- do.call(compute_per_substate_pca,
                  c(list(obj = obj, min_cells_per_pb = floor, pbs = pbs), pca_args))
  if (is.null(full)) { rm(obj, pbs); gc(verbose = FALSE); return(NULL) }
  # The object stays in memory across folds: compute_per_substate_pca still
  # reads obj[[]] for the sample -> subject / Etiology / Cohort lookups even
  # when pseudobulks are injected. Only AggregateExpression is skipped, which is
  # the expensive part.

  sample_to_subject <- stats::setNames(full$scores$subject, full$scores$sample)
  sample_to_group   <- stats::setNames(full$scores$Phenotype_2, full$scores$sample)
  subjects <- sort(unique(full$scores$subject))
  log_message("  jackknife[", target, "]: ", length(subjects), " folds")

  rows <- list()
  for (sbj in subjects) {
    drop <- names(sample_to_subject)[sample_to_subject == sbj]
    sub_pbs <- .loo_subset_pbs(pbs, drop)
    if (!length(sub_pbs)) next
    fold <- tryCatch(do.call(compute_per_substate_pca,
                             c(list(obj = obj, min_cells_per_pb = floor,
                                    pbs = sub_pbs), pca_args)),
                     error = function(e) NULL)
    if (is.null(fold)) next
    for (ss in unique(fold$scores$substate)) {
      a <- fold$scores[fold$scores$substate == ss, ]
      b <- full$scores[full$scores$substate == ss, ]
      j <- dplyr::inner_join(
        dplyr::select(a, sample = "sample", fold_pc1 = "PC1_oriented",
                      grp = "Phenotype_2"),
        dplyr::select(b, sample = "sample", full_pc1 = "PC1_oriented"),
        by = "sample")
      if (nrow(j) < 4L) next
      sa <- fold$significance[fold$significance$substate == ss, ]
      sf <- full$significance[full$significance$substate == ss, ]
      r <- suppressWarnings(stats::cor(j$fold_pc1, j$full_pc1))
      va <- fold$variance$var_explained[fold$variance$substate == ss &
                                          fold$variance$PC == "PC1"]
      d_fold <- .loo_cohens_d(j$fold_pc1, j$grp)
      d_full <- .loo_cohens_d(j$full_pc1, j$grp)
      rows[[length(rows) + 1]] <- data.frame(
        fold_subject = sbj,
        fold_group = unname(sample_to_group[drop[1]]),
        target = target, substate = as.character(ss),
        n_NIU = if (nrow(sa)) sa$n_NIU[1] else NA_integer_,
        n_Viral = if (nrow(sa)) sa$n_Viral[1] else NA_integer_,
        var_explained_PC1 = if (length(va)) va[1] else NA_real_,
        t_statistic = if (nrow(sa)) sa$t_statistic[1] else NA_real_,
        df = if (nrow(sa)) sa$df[1] else NA_real_,
        p_value = if (nrow(sa)) sa$p_value[1] else NA_real_,
        q_value = if (nrow(sa)) sa$q_value[1] else NA_real_,
        separating = if (nrow(sa)) isTRUE(sa$separating[1]) else NA,
        separating_full = if (nrow(sf)) isTRUE(sf$separating[1]) else NA,
        cohens_d = d_fold, cohens_d_full = d_full,
        delta_d = if (is.finite(d_full)) abs(d_fold - d_full) else NA_real_,
        # Sign flip means the substate's dominant axis reversed when one subject
        # was removed. Far more serious than a p-value moving.
        pc1_cor_to_full = r, pc1_sign_flipped = isTRUE(r < 0),
        stringsAsFactors = FALSE)
    }
  }
  rm(obj, pbs); gc(verbose = FALSE)
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

.loo_bridge <- function(cfg, jk_dir) {
  jcfg <- cfg$pc1_jackknife %||% list()
  n_perm <- as.integer(jcfg$n_permutation %||% 1000L)
  pm <- get_target_paths(cfg, "myeloid")$results_tables
  pt <- get_target_paths(cfg, "tcell")$results_tables
  fm <- file.path(pm, "pca_subject_scores.csv")
  ft <- file.path(pt, "pca_subject_scores.csv")
  if (!file.exists(fm) || !file.exists(ft)) {
    log_message("  jackknife[bridge]: PCA score CSVs missing; skipping.")
    return(NULL)
  }
  myel <- utils::read.csv(fm, stringsAsFactors = FALSE)
  tcel <- utils::read.csv(ft, stringsAsFactors = FALSE)

  rows <- list()
  for (w in c("weighted", "unweighted")) {
    is_w <- identical(w, "weighted")
    m <- .aggregate_subject_pc1(myel, weighted = is_w) |>
      dplyr::rename(myeloid_pc1 = "PC1")
    t <- .aggregate_subject_pc1(tcel, weighted = is_w) |>
      dplyr::rename(tcell_pc1 = "PC1")
    j <- dplyr::inner_join(dplyr::select(m, "subject", "Phenotype_2", "myeloid_pc1"),
                           dplyr::select(t, "subject", "Phenotype_2", "tcell_pc1"),
                           by = c("subject", "Phenotype_2"))
    if (nrow(j) < 6L) next
    full <- .partial_correlation_battery(j$myeloid_pc1, j$tcell_pc1,
                                         j$Phenotype_2, 0L, n_perm, seed = 42L)
    for (i in seq_len(nrow(j))) {
      k <- j[-i, , drop = FALSE]
      if (length(unique(k$Phenotype_2)) < 2L) next
      # Vary the seed by fold: reusing seed 42 across all folds would make the
      # permutation nulls share a stream and understate fold-to-fold spread.
      f <- .partial_correlation_battery(k$myeloid_pc1, k$tcell_pc1,
                                        k$Phenotype_2, 0L, n_perm,
                                        seed = 42L + i)
      rows[[length(rows) + 1]] <- data.frame(
        fold_subject = j$subject[i], fold_group = j$Phenotype_2[i],
        weighting = w, n = f$n,
        pearson_r = f$partial_pearson_r, pearson_p = f$partial_pearson_p,
        permutation_p = f$partial_permutation_p,
        full_r = full$partial_pearson_r,
        delta_r_vs_full = f$partial_pearson_r - full$partial_pearson_r,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

.loo_plots <- function(pca_loo, br_loo, cfg) {
  vd <- (cfg$paths_stress_sensitivity %||%
           list(viz = "outputs/viz/stress_sensitivity"))$viz
  ensure_dir(vd)
  if (!is.null(br_loo) && nrow(br_loo)) {
    d <- br_loo[br_loo$weighting == "unweighted", , drop = FALSE]
    if (nrow(d))
      save_pdf_png(
        ggplot(d, aes(stats::reorder(.data$fold_subject, .data$delta_r_vs_full),
                      .data$delta_r_vs_full, fill = .data$fold_group)) +
          geom_col() + coord_flip() +
          geom_hline(yintercept = 0, colour = "grey40") +
          scale_fill_manual(values = c(NIU = "#E21F26", Viral = "#397FB9"),
                            name = NULL) +
          labs(title = "Figure 4E leave-one-subject-out influence",
               subtitle = sprintf("Full-data partial r = %.3f (unweighted). Bars are the change when that subject is dropped.",
                                  d$full_r[1]),
               x = NULL, y = "change in partial r") +
          theme_bw(base_size = 10),
        file.path(vd, "pc1_bridge_loo_tornado"),
        w = 7, h = max(4, 0.22 * nrow(d) + 2))
  }
  if (!is.null(pca_loo) && nrow(pca_loo))
    save_pdf_png(
      ggplot(pca_loo[pca_loo$separating_full %in% TRUE, , drop = FALSE],
             aes(.data$substate, -log10(.data$q_value))) +
        geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.7) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                   colour = "#E21F26") +
        facet_wrap(~ target, scales = "free_x") +
        labs(title = "Per-substate PC1 significance across leave-one-out folds",
             subtitle = "Restricted to substates that separate in the full data. Dashed line is q = 0.05.",
             x = "substate", y = expression(-log[10](q))) +
        theme_bw(base_size = 10),
      file.path(vd, "pca_pc1_loo_stability"), w = 9, h = 4.5)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Pseudobulk-floor sensitivity of the F4E coupling
# ---------------------------------------------------------------------------

# Disclosure analysis, added 2026-07-31. The myeloid pseudobulk floor materially
# changes the reported coupling:
#
#   floor 10 -> n = 24, partial r = 0.423, CI [0.050, 0.734], perm p = 0.055
#   floor 20 -> n = 22, partial r = 0.700, CI [0.439, 0.841], perm p = 0.001
#
# Both are defensible readings of the same data. The higher floor gives more
# cells per pseudobulk and therefore a less noisy per-subject PC1; because
# measurement error in a predictor attenuates a correlation toward zero, the
# floor-20 estimate is the less attenuated one rather than simply the larger
# one. That is the argument for using it in the bridge, and this sweep is the
# evidence for it — plotted, with CIs, across the whole range, so the dependence
# is something we report rather than something a reviewer discovers.
#
# Reviewer 1 Major 3 asked exactly this class of question about Figure 4E, so
# this belongs in the supplement regardless of which floor is primary.
run_pc1_floor_sensitivity <- function(cfg, floors = NULL) {
  jcfg   <- cfg$pc1_jackknife %||% list()
  floors <- as.integer(floors %||% jcfg$floor_sweep %||% c(5, 10, 15, 20, 25, 30))
  n_boot <- as.integer(cfg$cross_compartment_bridge$n_bootstrap %||% 2000L)
  n_perm <- as.integer(cfg$cross_compartment_bridge$n_permutation %||% 1000L)
  out_dir <- file.path((cfg$paths_stress_sensitivity %||%
                          list(tables = "outputs/tables/stress_sensitivity"))$tables,
                       "floor_sensitivity")
  ensure_dir(out_dir)
  log_message("=== pc1_floor_sensitivity (myeloid pseudobulk floor sweep) ===")

  pm <- get_target_paths(cfg, "myeloid")
  op <- file.path(pm$results_objects, "IntegratedSeuratObject.rds")
  ft <- file.path(get_target_paths(cfg, "tcell")$results_tables,
                  "pca_subject_scores.csv")
  if (!file.exists(op) || !file.exists(ft)) {
    log_message("  floor sweep: myeloid object or T cell scores missing.")
    return(invisible(NULL))
  }
  cpcfg <- cfg$compartment_pca %||% list()
  obj  <- readRDS(op)
  tcel <- utils::read.csv(ft, stringsAsFactors = FALSE)

  rows <- list()
  for (fl in floors) {
    log_message("  floor ", fl)
    pbs <- build_per_substate_pseudobulks(obj, min_cells_per_pb = fl)
    if (!length(pbs)) next
    res <- tryCatch(compute_per_substate_pca(
      obj, min_cells_per_pb = fl,
      min_gene_count = as.integer(cpcfg$min_gene_count %||% 10L),
      hvg_n          = as.integer(cpcfg$hvg_n %||% 2000L),
      n_pcs          = as.integer(cpcfg$n_pcs %||% 5L),
      vst_blind      = isTRUE(cpcfg$vst_blind),
      pc1_split_fdr  = as.numeric(cpcfg$pc1_split_fdr %||% 0.05),
      pbs = pbs), error = function(e) NULL)
    if (is.null(res)) next
    for (w in c("weighted", "unweighted")) {
      is_w <- identical(w, "weighted")
      m <- .aggregate_subject_pc1(res$scores, weighted = is_w) |>
        dplyr::rename(myeloid_pc1 = "PC1")
      t <- .aggregate_subject_pc1(tcel, weighted = is_w) |>
        dplyr::rename(tcell_pc1 = "PC1")
      j <- dplyr::inner_join(
        dplyr::select(m, "subject", "Phenotype_2", "myeloid_pc1"),
        dplyr::select(t, "subject", "Phenotype_2", "tcell_pc1"),
        by = c("subject", "Phenotype_2"))
      if (nrow(j) < 6L) next
      b <- .partial_correlation_battery(j$myeloid_pc1, j$tcell_pc1,
                                        j$Phenotype_2, n_boot, n_perm, seed = 42L)
      rows[[length(rows) + 1]] <- data.frame(
        min_cells_per_pb = fl, weighting = w, n_subjects = b$n,
        n_pseudobulks = nrow(res$scores),
        n_substates = length(unique(res$scores$substate)),
        n_separating = sum(res$significance$separating %in% TRUE),
        separating_substates = paste(
          res$significance$substate[res$significance$separating %in% TRUE],
          collapse = ";"),
        partial_r = b$partial_pearson_r,
        ci_lo = b$partial_pearson_ci_lo, ci_hi = b$partial_pearson_ci_hi,
        parametric_p = b$partial_pearson_p,
        permutation_p = b$partial_permutation_p,
        stringsAsFactors = FALSE)
    }
    rm(pbs, res); gc(verbose = FALSE)
  }
  rm(obj); gc(verbose = FALSE)
  if (!length(rows)) return(invisible(NULL))
  sweep <- dplyr::bind_rows(rows)
  .sens_write(sweep, out_dir, "pc1_bridge_floor_sensitivity")

  vd <- (cfg$paths_stress_sensitivity %||%
           list(viz = "outputs/viz/stress_sensitivity"))$viz
  ensure_dir(vd)
  tryCatch(save_pdf_png(
    ggplot(sweep, aes(.data$min_cells_per_pb, .data$partial_r)) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_ribbon(aes(ymin = .data$ci_lo, ymax = .data$ci_hi), alpha = 0.18,
                  fill = "#397FB9") +
      geom_line(colour = "#397FB9") +
      geom_point(aes(shape = .data$permutation_p < 0.05), size = 2.6,
                 colour = "#397FB9") +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                         name = "permutation p < 0.05") +
      facet_wrap(~ weighting) +
      scale_x_continuous(breaks = floors) +
      labs(title = "Figure 4E coupling versus the myeloid pseudobulk floor",
           subtitle = paste("Partial r controlling for Phenotype_2, with",
                            "bootstrap 95% CI. Open points are not significant",
                            "by the permutation null."),
           x = "minimum cells per (subject, substate) pseudobulk",
           y = "partial Pearson r") +
      theme_bw(base_size = 10),
    file.path(vd, "pc1_bridge_floor_sensitivity"), w = 9, h = 4.5),
    error = function(e) log_message("  floor sweep figure failed: ",
                                    conditionMessage(e)))

  for (i in seq_len(nrow(sweep)))
    log_message(sprintf(
      "  floor %2d [%s]: n = %d, r = %.3f, CI [%.3f, %.3f], perm p = %.4g, separating = {%s}",
      sweep$min_cells_per_pb[i], sweep$weighting[i], sweep$n_subjects[i],
      sweep$partial_r[i], sweep$ci_lo[i], sweep$ci_hi[i],
      sweep$permutation_p[i], sweep$separating_substates[i]))
  invisible(sweep)
}

# ---------------------------------------------------------------------------
# F4F per-pair permutation-resolution stability
# ---------------------------------------------------------------------------

# Added 2026-08-01. The submitted Figure 4F reported six significant
# (myeloid substate x T cell substate) pairs, computed with 500 permutations.
# One of them -- myeloid 3 x T cell 1, r = -0.619, n = 11 -- was significant at
# permutation p = 0.0499. With 500 permutations the smallest representable
# p-value step is 1/501, so 0.0499 sits one draw from the threshold: the call is
# a coin flip against the resolution of the null, not a stable result. At 1000
# permutations it no longer clears 0.05, and a different pair (myeloid 4 x
# T cell 1, r = 0.762) enters instead.
#
# Per the 2026-08-01 decision the pair is RETAINED and reported with an explicit
# stability caveat rather than dropped. This function produces the evidence for
# that caveat: every pair evaluated at a ladder of permutation counts, with the
# p-value trajectory and a flag for pairs whose significance call changes. That
# is a more honest object than either silently keeping or silently dropping it.
run_pc1_pair_stability <- function(cfg, n_perm_grid = NULL, seeds = NULL) {
  n_perm_grid <- as.integer(n_perm_grid %||%
                              cfg$pc1_jackknife$pair_perm_grid %||%
                              c(500, 1000, 2000, 5000))
  # Several seeds at the smallest grid point show how much of the original call
  # was seed luck as opposed to iteration count.
  seeds <- as.integer(seeds %||% c(42, 43, 44))
  out_dir <- file.path((cfg$paths_stress_sensitivity %||%
                          list(tables = "outputs/tables/stress_sensitivity"))$tables,
                       "pair_stability")
  ensure_dir(out_dir)
  log_message("=== pc1_pair_stability (F4F permutation resolution) ===")

  pm <- get_target_paths(cfg, "myeloid")$results_tables
  pt <- get_target_paths(cfg, "tcell")$results_tables
  fm <- file.path(pm, "pca_subject_scores.csv")
  ft <- file.path(pt, "pca_subject_scores.csv")
  if (!file.exists(fm) || !file.exists(ft)) {
    log_message("  pair stability: PCA score CSVs missing.")
    return(invisible(NULL))
  }
  myel <- utils::read.csv(fm, stringsAsFactors = FALSE)
  tcel <- utils::read.csv(ft, stringsAsFactors = FALSE)

  rows <- list()
  for (np in n_perm_grid) {
    sd_use <- if (np == min(n_perm_grid)) seeds else seeds[1]
    for (sd in sd_use) {
      ps <- .per_substate_bridge(myel, tcel, n_boot = 0L, n_perm = np,
                                 seed = sd)
      if (!nrow(ps)) next
      ps$n_permutation <- np
      ps$seed <- sd
      rows[[length(rows) + 1]] <- ps
    }
  }
  if (!length(rows)) return(invisible(NULL))
  st <- dplyr::bind_rows(rows)
  st$significant <- !is.na(st$partial_permutation_p) &
                      st$partial_permutation_p < 0.05
  st$pair <- paste0("myeloid ", st$myeloid_substate, " x tcell ",
                    st$tcell_substate)

  summ <- st |>
    dplyr::group_by(.data$pair, .data$myeloid_substate, .data$tcell_substate) |>
    dplyr::summarise(
      n = .data$n[1],
      r = .data$partial_pearson_r[1],
      n_settings = dplyr::n(),
      n_significant = sum(.data$significant),
      min_p = min(.data$partial_permutation_p, na.rm = TRUE),
      max_p = max(.data$partial_permutation_p, na.rm = TRUE),
      # A pair is stable only if every resolution agrees. The ones that flip are
      # exactly the ones a reviewer should be told about.
      call_stable = dplyr::n() == sum(.data$significant) |
                      sum(.data$significant) == 0L,
      .groups = "drop") |>
    dplyr::arrange(dplyr::desc(abs(.data$r)))

  .sens_write(st,   out_dir, "pc1_pair_permutation_ladder")
  .sens_write(summ, out_dir, "pc1_pair_stability_summary")

  unstable <- summ[!summ$call_stable, , drop = FALSE]
  if (nrow(unstable))
    for (i in seq_len(nrow(unstable)))
      log_message(sprintf(
        "  UNSTABLE CALL: %s (r = %.3f, n = %d) significant in %d of %d permutation settings, p range [%.4f, %.4f]",
        unstable$pair[i], unstable$r[i], unstable$n[i],
        unstable$n_significant[i], unstable$n_settings[i],
        unstable$min_p[i], unstable$max_p[i]))
  else log_message("  All pair calls agree across permutation resolutions.")
  invisible(summ)
}

run_pc1_jackknife <- function(cfg, targets = NULL) {
  jcfg <- cfg$pc1_jackknife %||% list()
  targets <- as.character(targets %||% jcfg$targets %||%
                            c("myeloid", "tcell", "bridge"))
  jk_dir <- file.path((cfg$paths_stress_sensitivity %||%
                         list(tables = "outputs/tables/stress_sensitivity"))$tables,
                      "jackknife")
  ensure_dir(jk_dir)
  log_message("=== pc1_jackknife (leave-one-subject-out) ===")

  pca_loo <- dplyr::bind_rows(Filter(Negate(is.null), lapply(
    intersect(targets, c("myeloid", "tcell")),
    function(tg) .loo_pca(cfg, tg, jk_dir))))

  if (nrow(pca_loo)) {
    .sens_write(pca_loo, jk_dir, "pca_pc1_loo_jackknife")
    summ <- pca_loo |>
      dplyr::group_by(.data$target, .data$substate) |>
      dplyr::summarise(
        n_folds = dplyr::n(),
        separating_full = any(.data$separating_full %in% TRUE),
        frac_folds_separating = mean(.data$separating %in% TRUE),
        min_q = suppressWarnings(min(.data$q_value, na.rm = TRUE)),
        max_q = suppressWarnings(max(.data$q_value, na.rm = TRUE)),
        min_abs_pc1_cor = suppressWarnings(min(abs(.data$pc1_cor_to_full),
                                               na.rm = TRUE)),
        n_sign_flips = sum(.data$pc1_sign_flipped %in% TRUE),
        max_abs_delta_d = suppressWarnings(max(.data$delta_d, na.rm = TRUE)),
        most_influential_subject =
          .data$fold_subject[which.max(replace(.data$delta_d,
                                               is.na(.data$delta_d), -Inf))],
        most_influential_delta_d = suppressWarnings(max(.data$delta_d,
                                                        na.rm = TRUE)),
        .groups = "drop")
    summ$robust <- summ$frac_folds_separating >=
        as.numeric(jcfg$min_fold_separating_frac %||% 0.90) &
      summ$n_sign_flips == 0L &
      summ$min_abs_pc1_cor >= as.numeric(jcfg$min_abs_pc1_cor %||% 0.80)
    .sens_write(summ, jk_dir, "pca_pc1_loo_summary")
    frail <- summ[summ$separating_full %in% TRUE & !summ$robust, , drop = FALSE]
    if (nrow(frail))
      for (i in seq_len(nrow(frail)))
        log_message(sprintf(
          "  FRAGILE: %s substate %s — %.0f%% of folds separating, %d sign flip(s), min |r| to full = %.2f",
          frail$target[i], frail$substate[i],
          100 * frail$frac_folds_separating[i], frail$n_sign_flips[i],
          frail$min_abs_pc1_cor[i]))
    else log_message("  All full-data separating substates passed the LOO criteria.")
  }

  br_loo <- if ("bridge" %in% targets) .loo_bridge(cfg, jk_dir) else NULL
  if (!is.null(br_loo) && nrow(br_loo)) {
    .sens_write(br_loo, jk_dir, "pc1_bridge_loo_jackknife")
    bs <- br_loo |>
      dplyr::group_by(.data$weighting) |>
      dplyr::summarise(
        n_folds = dplyr::n(), full_r = .data$full_r[1],
        min_r = min(.data$pearson_r, na.rm = TRUE),
        max_r = max(.data$pearson_r, na.rm = TRUE),
        frac_folds_perm_p_lt_05 = mean(.data$permutation_p < 0.05, na.rm = TRUE),
        max_abs_delta_r = max(abs(.data$delta_r_vs_full), na.rm = TRUE),
        most_influential_subject =
          .data$fold_subject[which.max(abs(.data$delta_r_vs_full))],
        .groups = "drop")
    # Dropping from the smaller Viral arm perturbs the residualization
    # asymmetrically, so report the two directions separately.
    bg <- br_loo |>
      dplyr::group_by(.data$weighting, .data$fold_group) |>
      dplyr::summarise(n_folds = dplyr::n(),
                       min_r = min(.data$pearson_r, na.rm = TRUE),
                       max_abs_delta_r = max(abs(.data$delta_r_vs_full),
                                             na.rm = TRUE),
                       .groups = "drop")
    bs$robust <- bs$min_r >= as.numeric(jcfg$bridge_min_r %||% 0.40) &
      bs$frac_folds_perm_p_lt_05 >= as.numeric(jcfg$bridge_min_frac_perm_p %||% 0.95)
    .sens_write(bs, jk_dir, "pc1_bridge_loo_summary")
    .sens_write(bg, jk_dir, "pc1_bridge_loo_by_group")
    for (i in seq_len(nrow(bs)))
      log_message(sprintf(
        "  BRIDGE LOO [%s]: full r = %.3f, fold range [%.3f, %.3f], %.0f%% of folds perm p < 0.05, most influential = %s. Robust: %s",
        bs$weighting[i], bs$full_r[i], bs$min_r[i], bs$max_r[i],
        100 * bs$frac_folds_perm_p_lt_05[i], bs$most_influential_subject[i],
        bs$robust[i]))
  }

  tryCatch(.loo_plots(pca_loo, br_loo, cfg),
           error = function(e)
             log_message("  jackknife figures failed: ", conditionMessage(e)))
  invisible(list(pca = pca_loo, bridge = br_loo))
}
