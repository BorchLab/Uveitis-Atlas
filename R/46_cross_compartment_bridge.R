# R/46_cross_compartment_bridge.R
# Per-subject cross-compartment PC1 bridge — Figure 4 panel F.
#
# Mechanistic hinge of F4: if the myeloid compartment is instructing T cells
# differently in autoimmune vs viral disease, the per-subject myeloid PC1
# score (autoimmune-pole vs viral-pole) should correlate with the per-subject
# T cell PC1 score. Pearson is reported pooled across all subjects and
# stratified within etiology, with 2,000-iter bootstrap 95% CI and a 1,000-
# permutation null so the small-n within-etiology calls (Viral n ~ 12) are
# defensible against reviewer pushback on correlation stability.
#
# Inputs:
#   outputs/tables/eye/myeloid/pca_subject_scores.csv
#   outputs/tables/eye/tcell/pca_subject_scores.csv
# Outputs:
#   outputs/tables/cross_compartment/pc1_bridge_scores.csv
#   outputs/tables/cross_compartment/pc1_bridge_correlation.csv
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Collapse per-(sample, substate) PC1_oriented down to a single subject-level
# score. Two strategies — both saved so the user can audit whether rare
# substates carry the disease signal (unweighted) or whether population
# averages dominate (weighted).
#
# Step 1: pseudobulks are per Subject_Timepoint × substate. Collapse multiple
# timepoints to one subject row by averaging PC1_oriented across timepoints
# (weights are n_cells when weighted, equal when not).
# Step 2: aggregate across substates within subject to one PC1 per subject.
.aggregate_subject_pc1 <- function(scores, weighted = TRUE) {
  if (nrow(scores) == 0L) return(tibble::tibble())
  # Step 1: timepoints -> subject_substate. Weight by n_cells when requested.
  # IMPORTANT: compute PC1 BEFORE the n_cells summary, because dplyr::summarise
  # evaluates expressions left-to-right and once `n_cells` is reassigned it
  # becomes a scalar of length 1 within the group — `weighted.mean` then sees
  # mismatched x/w lengths. Computing PC1 first keeps `.data$n_cells` as the
  # original per-row vector for the weight argument.
  st1 <- scores |>
    dplyr::group_by(.data$subject, .data$substate, .data$Phenotype_2) |>
    dplyr::summarise(
      PC1     = if (isTRUE(weighted))
                  stats::weighted.mean(.data$PC1_oriented,
                                       w = .data$n_cells, na.rm = TRUE)
                else mean(.data$PC1_oriented, na.rm = TRUE),
      n_cells = sum(.data$n_cells, na.rm = TRUE),
      .groups = "drop")
  # Step 2: substates -> subject. Same weighting rule, same evaluation order.
  st2 <- st1 |>
    dplyr::group_by(.data$subject, .data$Phenotype_2) |>
    dplyr::summarise(
      PC1     = if (isTRUE(weighted))
                  stats::weighted.mean(.data$PC1, w = .data$n_cells, na.rm = TRUE)
                else mean(.data$PC1, na.rm = TRUE),
      n_cells = sum(.data$n_cells, na.rm = TRUE),
      n_substates = dplyr::n(),
      .groups = "drop")
  st2
}

# Partial correlation of x and y controlling for a categorical covariate z.
# Regresses x ~ z and y ~ z separately (intercepts per level), then correlates
# the residuals. This is the within-group coupling: "after accounting for
# disease arm, does a subject's myeloid PC1 still predict their T cell PC1?"
# Same bootstrap-CI / permutation-null battery as the marginal correlation.
.partial_correlation_battery <- function(x, y, z, n_boot, n_perm, seed = 42L) {
  out <- tibble::tibble(
    n = length(x),
    partial_pearson_r = NA_real_,
    partial_pearson_p = NA_real_,
    partial_pearson_ci_lo = NA_real_,
    partial_pearson_ci_hi = NA_real_,
    partial_permutation_p = NA_real_
  )
  if (length(x) < 4 || length(unique(z)) < 2) return(out)
  rx <- tryCatch(stats::resid(stats::lm(x ~ factor(z))),
                 error = function(e) NULL)
  ry <- tryCatch(stats::resid(stats::lm(y ~ factor(z))),
                 error = function(e) NULL)
  if (is.null(rx) || is.null(ry)) return(out)
  ct <- tryCatch(stats::cor.test(rx, ry, method = "pearson"),
                 error = function(e) NULL)
  if (!is.null(ct)) {
    out$partial_pearson_r <- unname(ct$estimate)
    out$partial_pearson_p <- ct$p.value
  }
  set.seed(seed)
  if (n_boot >= 100L) {
    n <- length(rx)
    rs <- vapply(seq_len(n_boot), function(i) {
      idx <- sample.int(n, size = n, replace = TRUE)
      suppressWarnings(stats::cor(rx[idx], ry[idx],
                                  use = "pairwise.complete.obs"))
    }, numeric(1))
    rs <- rs[is.finite(rs)]
    if (length(rs) >= 10L) {
      qq <- stats::quantile(rs, probs = c(0.025, 0.975), na.rm = TRUE)
      out$partial_pearson_ci_lo <- unname(qq[1])
      out$partial_pearson_ci_hi <- unname(qq[2])
    }
  }
  if (n_perm >= 100L && !is.na(out$partial_pearson_r)) {
    # Permute within each level of z so the null preserves group means.
    obs <- abs(out$partial_pearson_r)
    z_chr <- as.character(z)
    null_rs <- vapply(seq_len(n_perm), function(i) {
      ry_perm <- ry
      for (lvl in unique(z_chr)) {
        m <- z_chr == lvl
        ry_perm[m] <- sample(ry_perm[m])
      }
      suppressWarnings(abs(stats::cor(rx, ry_perm,
                                      use = "pairwise.complete.obs")))
    }, numeric(1))
    null_rs <- null_rs[is.finite(null_rs)]
    if (length(null_rs) >= 10L) {
      out$partial_permutation_p <- (sum(null_rs >= obs) + 1L) /
                                    (length(null_rs) + 1L)
    }
  }
  out
}

# Per-(myeloid substate × T cell substate) bridge. For each pair, builds a
# subject-level join on the per-substate PC1_oriented values (averaged across
# timepoints) and computes the within-etiology + partial Pearson r. Returns
# a long tibble with one row per (myeloid_substate, tcell_substate, stratum).
.per_substate_bridge <- function(myel, tcel, n_boot, n_perm) {
  # Reduce timepoints to subject level by averaging PC1_oriented per
  # (subject, substate, Phenotype_2).
  collapse_tp <- function(df) {
    df |>
      dplyr::group_by(.data$subject, .data$substate, .data$Phenotype_2) |>
      dplyr::summarise(PC1 = mean(.data$PC1_oriented, na.rm = TRUE),
                       n_cells = sum(.data$n_cells, na.rm = TRUE),
                       .groups = "drop")
  }
  m <- collapse_tp(myel)
  t <- collapse_tp(tcel)
  out_rows <- list()
  for (ms in unique(m$substate)) {
    m_sub <- dplyr::filter(m, .data$substate == ms) |>
      dplyr::rename(myeloid_pc1 = "PC1", n_cells_m = "n_cells")
    for (ts in unique(t$substate)) {
      t_sub <- dplyr::filter(t, .data$substate == ts) |>
        dplyr::rename(tcell_pc1 = "PC1", n_cells_t = "n_cells")
      jn <- dplyr::inner_join(m_sub |> dplyr::select(-.data$substate),
                              t_sub |> dplyr::select(-.data$substate),
                              by = c("subject", "Phenotype_2"))
      if (nrow(jn) < 4L) next
      row <- .partial_correlation_battery(jn$myeloid_pc1, jn$tcell_pc1,
                                          jn$Phenotype_2, n_boot, n_perm,
                                          seed = 42L)
      row$myeloid_substate <- ms
      row$tcell_substate   <- ts
      out_rows[[paste(ms, ts, sep = "::")]] <- row
    }
  }
  if (length(out_rows) == 0L) return(tibble::tibble())
  dplyr::bind_rows(out_rows)
}

# Bootstrap + permutation in one helper so both strata reuse the same nulls.
# Returns a single-row tibble of {n, pearson_r, pearson_p, ci_lo, ci_hi,
# spearman_r, spearman_p, permutation_p, slope, intercept}.
.correlation_battery <- function(x, y, n_boot, n_perm, seed = 42L) {
  out <- tibble::tibble(
    n = length(x),
    pearson_r = NA_real_, pearson_p = NA_real_,
    pearson_ci_lo = NA_real_, pearson_ci_hi = NA_real_,
    spearman_r = NA_real_, spearman_p = NA_real_,
    permutation_p = NA_real_,
    slope = NA_real_, intercept = NA_real_
  )
  if (length(x) < 3) return(out)

  pt <- tryCatch(stats::cor.test(x, y, method = "pearson"),
                 error = function(e) NULL)
  st <- tryCatch(suppressWarnings(stats::cor.test(x, y, method = "spearman")),
                 error = function(e) NULL)
  if (!is.null(pt)) {
    out$pearson_r <- unname(pt$estimate)
    out$pearson_p <- pt$p.value
  }
  if (!is.null(st)) {
    out$spearman_r <- unname(st$estimate)
    out$spearman_p <- st$p.value
  }
  lm_fit <- tryCatch(stats::lm(y ~ x), error = function(e) NULL)
  if (!is.null(lm_fit)) {
    coefs <- stats::coef(lm_fit)
    out$intercept <- unname(coefs[1])
    out$slope     <- unname(coefs[2])
  }

  set.seed(seed)
  if (n_boot >= 100L) {
    n <- length(x)
    rs <- vapply(seq_len(n_boot), function(i) {
      idx <- sample.int(n, size = n, replace = TRUE)
      suppressWarnings(stats::cor(x[idx], y[idx],
                                  use = "pairwise.complete.obs"))
    }, numeric(1))
    rs <- rs[is.finite(rs)]
    if (length(rs) >= 10L) {
      qq <- stats::quantile(rs, probs = c(0.025, 0.975), na.rm = TRUE)
      out$pearson_ci_lo <- unname(qq[1])
      out$pearson_ci_hi <- unname(qq[2])
    }
  }
  if (n_perm >= 100L && !is.na(out$pearson_r)) {
    obs <- abs(out$pearson_r)
    null_rs <- vapply(seq_len(n_perm), function(i) {
      suppressWarnings(abs(stats::cor(x, sample(y),
                                      use = "pairwise.complete.obs")))
    }, numeric(1))
    null_rs <- null_rs[is.finite(null_rs)]
    if (length(null_rs) >= 10L) {
      out$permutation_p <- (sum(null_rs >= obs) + 1L) / (length(null_rs) + 1L)
    }
  }
  out
}

# Main entry called from run_pipeline.R Phase 1d.
run_cross_compartment_pc1_bridge <- function(cfg,
                                             n_boot = NULL,
                                             n_perm = NULL) {
  bcfg <- cfg$cross_compartment_bridge %||% list()
  n_boot <- as.integer(n_boot %||% bcfg$n_bootstrap %||% 2000L)
  n_perm <- as.integer(n_perm %||% bcfg$n_permutation %||% 1000L)
  weight_default <- isTRUE(bcfg$weight_by_cluster_size %||% TRUE)
  also_unweighted <- isTRUE(bcfg$also_save_unweighted %||% TRUE)

  paths_myel <- get_target_paths(cfg, "myeloid")
  paths_tcel <- get_target_paths(cfg, "tcell")
  cc_paths   <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/",
    viz    = "outputs/viz/cross_compartment/")
  ensure_dir(cc_paths$tables)

  myel_csv <- file.path(paths_myel$results_tables, "pca_subject_scores.csv")
  tcel_csv <- file.path(paths_tcel$results_tables, "pca_subject_scores.csv")
  if (!file.exists(myel_csv) || !file.exists(tcel_csv)) {
    log_message("cross_compartment_bridge: PCA score CSV missing (",
                myel_csv, " | ", tcel_csv,
                "). Run compartment_pca first.")
    return(invisible(FALSE))
  }
  myel <- utils::read.csv(myel_csv, stringsAsFactors = FALSE)
  tcel <- utils::read.csv(tcel_csv, stringsAsFactors = FALSE)

  weightings <- if (also_unweighted) c("weighted", "unweighted") else
                  if (weight_default) "weighted" else "unweighted"

  bridge_rows <- list()
  cor_rows    <- list()
  for (w in weightings) {
    is_w <- identical(w, "weighted")
    m_subj <- .aggregate_subject_pc1(myel, weighted = is_w) |>
      dplyr::rename(myeloid_pc1 = "PC1", n_cells_myeloid = "n_cells",
                    n_substates_myeloid = "n_substates")
    t_subj <- .aggregate_subject_pc1(tcel, weighted = is_w) |>
      dplyr::rename(tcell_pc1 = "PC1", n_cells_tcell = "n_cells",
                    n_substates_tcell = "n_substates")
    joined <- dplyr::inner_join(m_subj, t_subj,
                                by = c("subject", "Phenotype_2"))
    joined$weighting <- w
    bridge_rows[[w]] <- joined

    if (nrow(joined) < 3L) {
      log_message("  bridge[", w, "]: fewer than 3 paired subjects (",
                  nrow(joined), "); skipping correlations.")
      next
    }

    # Pooled and per-etiology strata. Granulomatous stratification is run
    # only when n >= 6 within (etiology, granulomatous).
    pool_row <- .correlation_battery(joined$myeloid_pc1, joined$tcell_pc1,
                                     n_boot, n_perm, seed = 42L)
    pool_row$stratum <- "pooled"
    pool_row$weighting <- w
    cor_rows[[paste(w, "pooled", sep = "::")]] <- pool_row

    # Partial correlation controlling for Phenotype_2 — this is the actual
    # within-patient mechanistic claim. The pooled raw r is inflated by the
    # NIU-Viral separation that PC1_orientation guarantees by construction;
    # partial r removes the between-group component.
    partial <- .partial_correlation_battery(joined$myeloid_pc1,
                                            joined$tcell_pc1,
                                            joined$Phenotype_2,
                                            n_boot, n_perm, seed = 42L)
    partial_row <- pool_row  # inherit the n
    partial_row$pearson_r     <- partial$partial_pearson_r
    partial_row$pearson_p     <- partial$partial_pearson_p
    partial_row$pearson_ci_lo <- partial$partial_pearson_ci_lo
    partial_row$pearson_ci_hi <- partial$partial_pearson_ci_hi
    partial_row$permutation_p <- partial$partial_permutation_p
    partial_row$spearman_r    <- NA_real_  # not computed for partial
    partial_row$spearman_p    <- NA_real_
    partial_row$slope         <- NA_real_
    partial_row$intercept     <- NA_real_
    partial_row$stratum       <- "partial_controlling_Phenotype_2"
    partial_row$weighting     <- w
    cor_rows[[paste(w, "partial", sep = "::")]] <- partial_row

    for (eti in unique(joined$Phenotype_2)) {
      idx <- joined$Phenotype_2 == eti
      if (sum(idx) < 3L) next
      row <- .correlation_battery(joined$myeloid_pc1[idx],
                                  joined$tcell_pc1[idx],
                                  n_boot, n_perm, seed = 42L)
      row$stratum   <- paste0("etiology:", eti)
      row$weighting <- w
      cor_rows[[paste(w, eti, sep = "::")]] <- row
    }
  }

  bridge <- if (length(bridge_rows)) dplyr::bind_rows(bridge_rows) else
              tibble::tibble()
  correlations <- if (length(cor_rows)) dplyr::bind_rows(cor_rows) else
              tibble::tibble()
  if (nrow(correlations) > 0L) {
    correlations <- correlations |>
      dplyr::select("stratum", "weighting", "n",
                    "pearson_r", "pearson_p",
                    "pearson_ci_lo", "pearson_ci_hi",
                    "spearman_r", "spearman_p",
                    "permutation_p", "slope", "intercept")
  }

  utils::write.csv(bridge,
                   file.path(cc_paths$tables, "pc1_bridge_scores.csv"),
                   row.names = FALSE)
  utils::write.csv(correlations,
                   file.path(cc_paths$tables, "pc1_bridge_correlation.csv"),
                   row.names = FALSE)
  log_message("cross_compartment_bridge: wrote ", nrow(bridge),
              " subject-weighting rows and ", nrow(correlations),
              " correlation rows.")

  # Per-(myeloid substate × T cell substate) bridge. Reads the same per-
  # substate PCA scores and computes a partial correlation per pair so the
  # F4_F_per_substate heatmap can show which sender/receiver pairings carry
  # the cross-compartment coupling signal at substate granularity.
  per_substate <- .per_substate_bridge(myel, tcel, n_boot = max(500L, n_boot %/% 4L),
                                       n_perm = max(500L, n_perm %/% 2L))
  if (nrow(per_substate) > 0L) {
    out_path <- file.path(cc_paths$tables, "pc1_bridge_per_substate.csv")
    utils::write.csv(per_substate, out_path, row.names = FALSE)
    log_message("cross_compartment_bridge: wrote ", nrow(per_substate),
                " per-(myeloid x tcell) substate rows to ", out_path)
  }
  invisible(TRUE)
}
