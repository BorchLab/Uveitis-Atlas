# R/57_bcell_cross_compartment_bridge.R
# B cell-centered cross-compartment PC1 bridges. Mirrors R/46's myeloid x T
# cell bridge for two new pairings:
#   bcell <-> myeloid
#   bcell <-> tcell
# so the B cell story can be tied to the myeloid + T cell PC1 axes the same
# way Figure 4 ties myeloid -> T cell.
#
# Reuses R/46's helpers verbatim (they are pair-agnostic at the math level):
#   .aggregate_subject_pc1     subject-level PC1 aggregation
#   .correlation_battery       pooled / per-arm Pearson r with boot + perm
#   .partial_correlation_battery  partial r controlling Phenotype_2
#   .per_substate_bridge       per-(source x target) substate partial r
#
# Inputs:
#   outputs/tables/eye/<bcell|myeloid|tcell>/pca_subject_scores.csv
# Outputs (under outputs/tables/cross_compartment/):
#   pc1_bridge_bcell_myeloid_scores.csv
#   pc1_bridge_bcell_myeloid_correlation.csv
#   pc1_bridge_bcell_myeloid_per_substate.csv
#   pc1_bridge_bcell_tcell_scores.csv
#   pc1_bridge_bcell_tcell_correlation.csv
#   pc1_bridge_bcell_tcell_per_substate.csv
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Core for one pair. partner = "myeloid" or "tcell"; B cell is always the
# "primary" side (column prefixes are bcell_/<partner>_).
.bcell_bridge_one_pair <- function(cfg, partner, n_boot, n_perm,
                                   weight_default, also_unweighted,
                                   cc_paths) {
  paths_b <- get_target_paths(cfg, "bcell")
  paths_p <- get_target_paths(cfg, partner)
  b_csv <- file.path(paths_b$results_tables, "pca_subject_scores.csv")
  p_csv <- file.path(paths_p$results_tables, "pca_subject_scores.csv")
  if (!file.exists(b_csv) || !file.exists(p_csv)) {
    log_message("bcell_bridge[", partner,
                "]: PCA score CSV missing (", b_csv, " | ", p_csv,
                "). Run compartment_pca first.")
    return(invisible(FALSE))
  }
  bcl <- utils::read.csv(b_csv, stringsAsFactors = FALSE)
  pnr <- utils::read.csv(p_csv, stringsAsFactors = FALSE)

  weightings <- if (also_unweighted) c("weighted", "unweighted") else
                  if (weight_default) "weighted" else "unweighted"

  bridge_rows <- list()
  cor_rows    <- list()
  for (w in weightings) {
    is_w <- identical(w, "weighted")
    b_subj <- .aggregate_subject_pc1(bcl, weighted = is_w) |>
      dplyr::rename(bcell_pc1 = "PC1", n_cells_bcell = "n_cells",
                    n_substates_bcell = "n_substates")
    p_subj <- .aggregate_subject_pc1(pnr, weighted = is_w) |>
      dplyr::rename(partner_pc1 = "PC1", n_cells_partner = "n_cells",
                    n_substates_partner = "n_substates")
    joined <- dplyr::inner_join(b_subj, p_subj,
                                by = c("subject", "Phenotype_2"))
    joined$weighting <- w
    joined$partner   <- partner
    bridge_rows[[w]] <- joined

    if (nrow(joined) < 3L) {
      log_message("  bcell_bridge[", partner, "][", w,
                  "]: fewer than 3 paired subjects (",
                  nrow(joined), "); skipping correlations.")
      next
    }
    # Pooled marginal r (kept as supplement — between-arm separation
    # contributes to the apparent coupling because PC1 is oriented).
    pool_row <- .correlation_battery(joined$bcell_pc1, joined$partner_pc1,
                                     n_boot, n_perm, seed = 42L)
    pool_row$stratum   <- "pooled"
    pool_row$weighting <- w
    pool_row$partner   <- partner
    cor_rows[[paste(w, "pooled", sep = "::")]] <- pool_row

    # Partial r controlling Phenotype_2 — the within-patient mechanistic claim.
    partial <- .partial_correlation_battery(joined$bcell_pc1,
                                            joined$partner_pc1,
                                            joined$Phenotype_2,
                                            n_boot, n_perm, seed = 42L)
    partial_row <- pool_row
    partial_row$pearson_r     <- partial$partial_pearson_r
    partial_row$pearson_p     <- partial$partial_pearson_p
    partial_row$pearson_ci_lo <- partial$partial_pearson_ci_lo
    partial_row$pearson_ci_hi <- partial$partial_pearson_ci_hi
    partial_row$permutation_p <- partial$partial_permutation_p
    partial_row$spearman_r    <- NA_real_
    partial_row$spearman_p    <- NA_real_
    partial_row$slope         <- NA_real_
    partial_row$intercept     <- NA_real_
    partial_row$stratum       <- "partial_controlling_Phenotype_2"
    partial_row$weighting     <- w
    partial_row$partner       <- partner
    cor_rows[[paste(w, "partial", sep = "::")]] <- partial_row

    for (eti in unique(joined$Phenotype_2)) {
      idx <- joined$Phenotype_2 == eti
      if (sum(idx) < 3L) next
      row <- .correlation_battery(joined$bcell_pc1[idx],
                                  joined$partner_pc1[idx],
                                  n_boot, n_perm, seed = 42L)
      row$stratum   <- paste0("etiology:", eti)
      row$weighting <- w
      row$partner   <- partner
      cor_rows[[paste(w, eti, sep = "::")]] <- row
    }
  }

  bridge <- if (length(bridge_rows)) dplyr::bind_rows(bridge_rows) else
              tibble::tibble()
  correlations <- if (length(cor_rows)) dplyr::bind_rows(cor_rows) else
              tibble::tibble()
  if (nrow(correlations) > 0L) {
    correlations <- correlations |>
      dplyr::select("stratum", "weighting", "partner", "n",
                    "pearson_r", "pearson_p",
                    "pearson_ci_lo", "pearson_ci_hi",
                    "spearman_r", "spearman_p",
                    "permutation_p", "slope", "intercept")
  }

  prefix <- paste0("pc1_bridge_bcell_", partner)
  utils::write.csv(bridge,
                   file.path(cc_paths$tables, paste0(prefix, "_scores.csv")),
                   row.names = FALSE)
  utils::write.csv(correlations,
                   file.path(cc_paths$tables,
                             paste0(prefix, "_correlation.csv")),
                   row.names = FALSE)
  log_message("bcell_bridge[", partner, "]: wrote ", nrow(bridge),
              " subject-weighting rows and ", nrow(correlations),
              " correlation rows.")

  # Per-(bcell substate x partner substate) partial r. .per_substate_bridge
  # is parameterized on column names "myeloid_substate" / "tcell_substate"
  # internally, so we rename its output to a pair-agnostic schema.
  per_substate <- .per_substate_bridge(bcl, pnr,
                                       n_boot = max(500L, n_boot %/% 4L),
                                       n_perm = max(500L, n_perm %/% 2L))
  if (nrow(per_substate) > 0L) {
    # .per_substate_bridge returns columns myeloid_substate / tcell_substate
    # by name — these are role labels (first arg / second arg). Rename so
    # the bcell role is explicit.
    per_substate <- per_substate |>
      dplyr::rename(bcell_substate   = "myeloid_substate",
                    partner_substate = "tcell_substate")
    per_substate$partner <- partner
    out_path <- file.path(cc_paths$tables,
                          paste0(prefix, "_per_substate.csv"))
    utils::write.csv(per_substate, out_path, row.names = FALSE)
    log_message("bcell_bridge[", partner, "]: wrote ", nrow(per_substate),
                " per-(bcell x ", partner, ") substate rows to ", out_path)
  }
  invisible(TRUE)
}

# Main entry — runs the two B cell pairings sequentially. Pipeline phase
# 1d (after compartment_pca + cross_compartment_bridge for myeloid x T cell).
run_bcell_cross_compartment_bridge <- function(cfg,
                                               n_boot = NULL,
                                               n_perm = NULL) {
  bcfg <- cfg$cross_compartment_bridge %||% list()
  n_boot <- as.integer(n_boot %||% bcfg$n_bootstrap %||% 2000L)
  n_perm <- as.integer(n_perm %||% bcfg$n_permutation %||% 1000L)
  weight_default  <- isTRUE(bcfg$weight_by_cluster_size %||% TRUE)
  also_unweighted <- isTRUE(bcfg$also_save_unweighted   %||% TRUE)

  cc_paths <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/",
    viz    = "outputs/viz/cross_compartment/")
  ensure_dir(cc_paths$tables)

  log_message("=== bcell_cross_compartment_bridge: bcell x {myeloid, tcell} ===")
  ok_m <- .bcell_bridge_one_pair(cfg, "myeloid", n_boot, n_perm,
                                 weight_default, also_unweighted, cc_paths)
  ok_t <- .bcell_bridge_one_pair(cfg, "tcell",   n_boot, n_perm,
                                 weight_default, also_unweighted, cc_paths)
  invisible(isTRUE(ok_m) || isTRUE(ok_t))
}
