# R/45_compartment_pca.R
# Per-compartment pseudobulk PCA. Refactored from the inline implementation
# in R/85_viz_myeloid.R:.myeloid_panelD_pca_facets so the same logic powers
# both F3 panel F (myeloid) and F4 panels D/E (T cell), and so subject-level
# PC1 scores live as a canonical on-disk CSV that the cross-compartment
# bridge (R/46) and LIANA gene-set construction (R/48) can read directly.
#
# Public entry points:
#   run_compartment_pca(cfg, target)         orchestrates one compartment
#   compute_per_substate_pca(obj, ...)       does the actual aggregation +
#                                             vst + prcomp + sign-orient,
#                                             returns scores / loadings
#                                             tibbles ready to write.
#   .pc1_loadings_by_program(target, programs, cfg)
#                                             pulls per-program PC1 loadings
#                                             restricted to disease-separating
#                                             substates and writes
#                                             pc1_loadings_by_program.csv.
#
# Outputs under outputs/tables/eye/<target>/:
#   pca_subject_scores.csv         one row per (substate, pseudobulk sample)
#   pca_gene_loadings.csv          one row per (substate, gene)
#   pca_variance_explained.csv     one row per (substate, PC)
#   pca_pc1_significance.csv       one row per substate (Welch t + BH q)
#   pc1_loadings_by_program.csv    one row per (substate, program, gene)
#
# Stress-sensitivity extension (2026-07-31). Three optional
# arguments make the disease axis re-derivable with an acute-stress /
# dissociation covariate removed:
#   covariate_col   name of a per-(substate, sample) column carried on the
#                   pseudobulk colData by build_per_substate_pseudobulks()
#   covariate_mode  none | rbe_group_protected | rbe_naive | hvg_drop
#   exclude_genes   gene vector for hvg_drop
# All three default to the published behaviour, and .pca_adjust_matrix() returns
# the input matrix untouched when covariate_mode == "none", so the published
# path is provably unchanged. See .pca_adjust_matrix for why a DESeq2 design of
# ~ covariate + group is NOT one of the offered arms.
#
# `pbs` lets a caller inject prebuilt pseudobulks so AggregateExpression is paid
# for once (R/85 panel D; the 34-fold jackknife in R/46c).
#
# Compatibility note: matches the F3 inline implementation's conventions
# exactly: DESeq2::vst(blind = FALSE, design = ~ group), HVG top-2000 by
# rowVars, prcomp(scale. = TRUE), PC1 sign-flipped so the Viral centroid is
# at positive PC1. This means the new module's PC1 scores reproduce the
# (formerly inline) F3 panel F values up to numerical noise.
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Asymmetric pseudobulk floor. Myeloid is sparser per (subject, substate) so
# it gets a lower floor than T / B cell. Sourced from cfg$compartment_pca to
# stay overrideable from config without code edits.
.pca_min_cells <- function(cfg, target) {
  cpcfg <- cfg$compartment_pca %||% list()
  key <- paste0("min_cells_per_subject_substate_", target)
  as.integer(cpcfg[[key]] %||% 20L)
}

# Rank guard shared by every model-fitting site in the stress-sensitivity work.
# In the eye arm Phenotype_2, Tissue_2 and Cohort are perfectly collinear
# (23 aqueous/US/NIU samples vs 13 vitreous/Japan/Viral), so any design that
# tries to hold two of them at once is unidentifiable. lm() would silently
# return NA for the aliased term; we would rather stop.
.assert_full_rank <- function(mm, what = "design") {
  r <- qr(mm)$rank
  if (r < ncol(mm))
    stop(what, ": design is rank-deficient (rank ", r, " < ", ncol(mm),
         " columns). In the eye arm Phenotype_2, Tissue_2 and Cohort are ",
         "perfectly collinear; no categorical site adjustment is identifiable.",
         call. = FALSE)
  invisible(TRUE)
}

# Variance inflation of a continuous covariate against the group term. A
# covariate that the group term already explains is a categorical adjustment
# wearing a numeric costume: it inflates the group SE and yields a stable-
# looking but meaningless coefficient. Callers abort above vif_abort.
.covariate_vif <- function(cov, grp) {
  ok <- is.finite(cov) & !is.na(grp)
  if (sum(ok) < 4L || length(unique(grp[ok])) < 2L) return(NA_real_)
  if (stats::sd(cov[ok]) <= .Machine$double.eps) return(Inf)
  r2 <- suppressWarnings(
    summary(stats::lm(cov[ok] ~ factor(as.character(grp[ok]))))$r.squared)
  if (!is.finite(r2)) return(NA_real_)
  if (r2 >= 1 - 1e-12) return(Inf)
  1 / (1 - r2)
}

# Variance-stabilise a pseudobulk count matrix, degrading gracefully instead of
# failing.
#
# DESeq2::vst() fits its dispersion trend on a subsample of `nsub` genes with
# mean normalized count > 5, and errors outright when a sparse substate cannot
# supply them ("less than 'nsub' rows with mean normalized count > 5"). For the
# well-populated substates that is never hit.
#
# The chain only ever engages AFTER the standard call has already errored, so
# every substate that currently transforms successfully is untouched and the
# published panels cannot move. Callers get `method` back so a fallback can be
# labelled in the figure rather than passed off as an ordinary vst.
.pca_vst_robust <- function(dds, blind = FALSE, min_gene_count = 10L) {
  out <- function(m, meth) list(mat = m, method = meth)
  r <- tryCatch(out(SummarizedExperiment::assay(DESeq2::vst(dds, blind = blind)),
                    "vst"),
                error = function(e) NULL)
  if (!is.null(r)) return(r)
  # Retry with a dispersion subsample small enough for a sparse substate.
  n_ok <- sum(rowMeans(DESeq2::counts(dds, normalized = FALSE)) > 5)
  r <- tryCatch(out(SummarizedExperiment::assay(
                      DESeq2::vst(dds, blind = blind,
                                  nsub = max(10L, min(1000L, n_ok - 1L)))),
                    "vst_reduced_nsub"), error = function(e) NULL)
  if (!is.null(r)) return(r)
  # Full parametric VST: no nsub subsampling, slower but robust at small n.
  r <- tryCatch(out(SummarizedExperiment::assay(
                      DESeq2::varianceStabilizingTransformation(dds, blind = blind)),
                    "vst_full"), error = function(e) NULL)
  if (!is.null(r)) return(r)
  # Last resort: log2 CPM. Not variance-stabilised, so label it in the figure.
  m <- DESeq2::counts(dds, normalized = FALSE)
  cs <- pmax(colSums(m), 1)
  out(log2(sweep(m, 2, cs, "/") * 1e6 + 1), "log2cpm_fallback")
}

# Remove the stress axis from a vst matrix before HVG selection and prcomp.
#
# Modes:
#   none                 return mat untouched (the published path)
#   rbe_group_protected  fit expr ~ group + cov, subtract only the cov term.
#                        Removes the WITHIN-group stress-associated variation,
#                        which is the only identifiable component. Primary.
#   rbe_naive            fit expr ~ 1 + cov, subtract the whole term including
#                        the part shared with group. Under perfect group/site
#                        confounding this deliberately deletes group signal, so
#                        it is a worst-case bound, never a primary result.
#   hvg_drop             drop the stress genes outright before variance ranking.
#                        No over-correction risk, but blind to indirect effects
#                        through correlated genes.
#
# Deliberately NOT offered: a DESeq2 design of ~ cov + group. vst(blind=FALSE)
# uses the design only for dispersion estimation and does not residualize the
# returned matrix, so that arm looks like an adjustment and does nothing.
#
# Returns list(mat=, status=, vif=, n_excluded=). status is "ok", "unadjusted"
# (mode none) or a reason string; callers record it as arm_failed.
.pca_adjust_matrix <- function(mat, grp, cov, mode,
                               exclude_genes = NULL, vif_abort = 10) {
  if (identical(mode, "none"))
    return(list(mat = mat, status = "unadjusted", vif = NA_real_,
                n_excluded = 0L))

  if (identical(mode, "hvg_drop")) {
    drop <- intersect(rownames(mat), unique(as.character(exclude_genes)))
    if (length(drop) >= nrow(mat) - 10L)
      return(list(mat = mat, status = "exclude_list_too_large",
                  vif = NA_real_, n_excluded = length(drop)))
    return(list(mat = mat[setdiff(rownames(mat), drop), , drop = FALSE],
                status = "ok", vif = NA_real_, n_excluded = length(drop)))
  }

  if (!requireNamespace("limma", quietly = TRUE))
    return(list(mat = mat, status = "limma_unavailable", vif = NA_real_,
                n_excluded = 0L))
  if (is.null(cov) || anyNA(cov) || length(cov) != ncol(mat))
    return(list(mat = mat, status = "covariate_missing_or_ragged",
                vif = NA_real_, n_excluded = 0L))
  if (stats::sd(cov) <= .Machine$double.eps)
    return(list(mat = mat, status = "covariate_zero_variance",
                vif = NA_real_, n_excluded = 0L))

  vif <- .covariate_vif(cov, grp)
  if (is.finite(vif) && vif > vif_abort)
    return(list(mat = mat, status = sprintf("near_aliased_vif_%.1f", vif),
                vif = vif, n_excluded = 0L))

  cov_z <- as.numeric(scale(cov))   # centered within substate
  adj <- tryCatch({
    if (identical(mode, "rbe_group_protected")) {
      dm <- stats::model.matrix(~ factor(as.character(grp)))
      .assert_full_rank(cbind(dm, cov_z), "pca rbe_group_protected")
      limma::removeBatchEffect(mat, covariates = cov_z, design = dm)
    } else {
      limma::removeBatchEffect(mat, covariates = cov_z)
    }
  }, error = function(e) NULL)
  if (is.null(adj))
    return(list(mat = mat, status = "removeBatchEffect_failed",
                vif = vif, n_excluded = 0L))
  list(mat = adj, status = "ok", vif = vif, n_excluded = 0L)
}

# Heart of the module. Given a Seurat compartment object, run per-substate
# pseudobulk PCA following the F3 convention. Returns a list with elements
#   scores      tibble (substate, sample, n_cells, group, PC1..PC5, PC1_oriented)
#   loadings    tibble (substate, gene, PC1..PC5, PC1_oriented, loading_rank_within_substate)
#   variance    tibble (substate, PC, var_explained, cum_var)
#   significance  tibble (substate, t_statistic, df, p_value, q_value, n_NIU, n_Viral, separating)
# Callers (run_compartment_pca, .myeloid_panelD_pca_facets) handle the writing
# and the plotting.
compute_per_substate_pca <- function(obj,
                                     cluster_col       = "knn.leiden.cluster",
                                     group_col         = "Phenotype_2",
                                     subject_col       = "Subject",
                                     sample_col        = NULL,
                                     groups            = c("NIU", "Viral"),
                                     min_cells_per_pb  = 30L,
                                     min_gene_count    = 10L,
                                     hvg_n             = 2000L,
                                     n_pcs             = 5L,
                                     vst_blind         = FALSE,
                                     pc1_split_fdr     = 0.05,
                                     covariate_col     = NULL,
                                     covariate_mode    = c("none",
                                                           "rbe_group_protected",
                                                           "rbe_naive",
                                                           "hvg_drop"),
                                     exclude_genes     = NULL,
                                     vif_abort         = 10,
                                     # Minimum pseudobulk columns for a substate
                                     # to be attempted. Default 4 preserves the
                                     # published behaviour; the completeness
                                     # supplement lowers it to 3 (prcomp needs
                                     # >= 2, PC2 needs >= 3) so no cluster is
                                     # dropped for being small.
                                     min_pb_cols       = 4L,
                                     robust_vst        = FALSE,
                                     pbs               = NULL) {
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("compute_per_substate_pca: DESeq2 required.")
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("compute_per_substate_pca: matrixStats required.")
  covariate_mode <- match.arg(covariate_mode)
  if (!identical(covariate_mode, "none") && is.null(covariate_col) &&
      !identical(covariate_mode, "hvg_drop"))
    stop("compute_per_substate_pca: covariate_mode='", covariate_mode,
         "' requires covariate_col.")

  if (is.null(sample_col)) {
    sample_col <- if ("Subject_Timepoint" %in% colnames(obj[[]]))
                    "Subject_Timepoint" else "orig.ident"
  }

  # pbs may be injected by a caller that already paid for AggregateExpression
  # (R/85 panel D, and the leave-one-out jackknife in R/46c, which subsets a
  # single cached build 34 times rather than re-aggregating per fold). When
  # injected, min_cells_per_pb has already been applied by the builder.
  if (is.null(pbs)) {
    pbs <- build_per_substate_pseudobulks(
      obj,
      cluster_col      = cluster_col,
      group_col        = group_col,
      groups           = groups,
      min_cells_per_pb = min_cells_per_pb,
      covariate_cols   = covariate_col
    )
  }
  if (length(pbs) == 0) {
    log_message("  compute_per_substate_pca: no pseudobulks after floor=",
                min_cells_per_pb, " filter; aborting.")
    return(NULL)
  }

  # Subject lookup keyed by the sample column used in the pseudobulks. When
  # sample_col is Subject_Timepoint, this collapses "<subject>_<visit>" to
  # "<subject>". Downstream bridge (R/46) further averages across timepoints.
  meta <- obj[[]]
  sample_to_subject <- tapply(as.character(meta[[subject_col]]),
                              as.character(meta[[sample_col]]),
                              function(x) x[1])
  sample_to_pheno   <- tapply(as.character(meta[[group_col]]),
                              as.character(meta[[sample_col]]),
                              function(x) x[1])
  sample_to_eti     <- if ("Etiology" %in% colnames(meta))
                         tapply(as.character(meta$Etiology),
                                as.character(meta[[sample_col]]),
                                function(x) x[1]) else NULL
  sample_to_gran    <- if ("Phenotype" %in% colnames(meta))
                         tapply(as.character(meta$Phenotype),
                                as.character(meta[[sample_col]]),
                                function(x) x[1]) else NULL
  sample_to_cohort  <- if ("Cohort" %in% colnames(meta))
                         tapply(as.character(meta$Cohort),
                                as.character(meta[[sample_col]]),
                                function(x) x[1]) else NULL

  scores_rows   <- list()
  loadings_rows <- list()
  variance_rows <- list()
  sig_rows      <- list()

  for (ck in names(pbs)) {
    pb <- pbs[[ck]]
    if (inherits(pb, "SummarizedExperiment")) {
      cd <- SummarizedExperiment::colData(pb)
      m  <- SummarizedExperiment::assay(pb, "counts")
    } else {
      cd <- pb$coldata; m <- pb$counts
    }
    grp <- factor(cd$group, levels = groups)
    if (ncol(m) < min_pb_cols || length(unique(grp)) < 2) {
      log_message("  PCA substate ", ck, ": <", min_pb_cols,
                  " columns or single group; skipping.")
      next
    }

    # Per-(substate, sample) covariate, computed by the pseudobulk builder from
    # exactly the cells that entered each column. Not a per-sample lookup.
    cov_vec <- if (!is.null(covariate_col) && covariate_col %in% colnames(cd))
                 as.numeric(cd[[covariate_col]]) else NULL

    fit <- tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(m),
        colData   = data.frame(group = grp),
        design    = ~ group)
      dds <- dds[rowSums(DESeq2::counts(dds)) > min_gene_count, ]
      tr <- if (isTRUE(robust_vst))
              .pca_vst_robust(dds, blind = vst_blind,
                              min_gene_count = min_gene_count)
            else list(mat = SummarizedExperiment::assay(
                        DESeq2::vst(dds, blind = vst_blind)), method = "vst")
      mat <- tr$mat
      vst_method <- tr$method
      # Adjust after the count filter and before variance ranking, so the HVG
      # set itself reflects the adjustment rather than being chosen on the
      # unadjusted matrix and only then de-stressed.
      adj <- .pca_adjust_matrix(mat, grp, cov_vec, covariate_mode,
                                exclude_genes = exclude_genes,
                                vif_abort = vif_abort)
      mat  <- adj$mat
      vars <- matrixStats::rowVars(mat)
      keep <- order(-vars)[seq_len(min(hvg_n, length(vars)))]
      mat  <- mat[keep, ]
      # scale. = TRUE fails on any zero-variance gene, which small substates can
      # produce after the fallback transform; drop those rather than lose the
      # whole substate.
      mat <- mat[matrixStats::rowVars(mat) > 0, , drop = FALSE]
      list(pca = stats::prcomp(t(mat), scale. = TRUE), adj = adj,
           vst_method = vst_method)
    }, error = function(e) {
      log_message("    PCA failed for substate ", ck, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(fit)) next
    pca_res    <- fit$pca
    adj_status <- fit$adj$status
    adj_vif    <- fit$adj$vif
    adj_nex    <- fit$adj$n_excluded
    if (!identical(covariate_mode, "none") && !identical(adj_status, "ok"))
      log_message("    PCA substate ", ck, ": arm '", covariate_mode,
                  "' not applied (", adj_status, "); reporting UNADJUSTED.")

    # Sign orient so positive PC1 is the Viral centroid (matches F3 convention).
    grp_chr   <- as.character(grp)
    viral_mean <- mean(pca_res$x[grp_chr == "Viral", 1], na.rm = TRUE)
    niu_mean   <- mean(pca_res$x[grp_chr == "NIU",   1], na.rm = TRUE)
    flip <- if (is.finite(viral_mean) && is.finite(niu_mean) &&
                viral_mean < niu_mean) -1 else 1
    pc1_oriented <- pca_res$x[, 1] * flip
    rotation_pc1 <- pca_res$rotation[, 1] * flip

    n_keep_pc <- min(n_pcs, ncol(pca_res$x))
    scores_mat <- pca_res$x[, seq_len(n_keep_pc), drop = FALSE]
    colnames(scores_mat) <- paste0("PC", seq_len(n_keep_pc))
    rot_mat <- pca_res$rotation[, seq_len(n_keep_pc), drop = FALSE]
    colnames(rot_mat) <- paste0("PC", seq_len(n_keep_pc))

    sample_id <- as.character(cd$sample)
    n_cells   <- as.integer(cd$n_cells)
    scores_df <- tibble::tibble(
      substate     = as.character(ck),
      sample       = sample_id,
      subject      = unname(sample_to_subject[sample_id]),
      Phenotype_2  = grp_chr,
      Etiology     = if (!is.null(sample_to_eti))    unname(sample_to_eti[sample_id])    else NA_character_,
      Phenotype    = if (!is.null(sample_to_gran))   unname(sample_to_gran[sample_id])   else NA_character_,
      Cohort       = if (!is.null(sample_to_cohort)) unname(sample_to_cohort[sample_id]) else NA_character_,
      n_cells      = n_cells
    )
    scores_df <- dplyr::bind_cols(scores_df,
                                  tibble::as_tibble(scores_mat),
                                  PC1_oriented = pc1_oriented)

    rot_df <- tibble::tibble(
      substate = as.character(ck),
      gene     = rownames(rot_mat)
    )
    rot_df <- dplyr::bind_cols(rot_df,
                               tibble::as_tibble(rot_mat),
                               PC1_oriented = rotation_pc1)
    rot_df$loading_rank_within_substate <- rank(-abs(rot_df$PC1_oriented),
                                                ties.method = "first")

    var_pct <- (pca_res$sdev[seq_len(n_keep_pc)]^2 / sum(pca_res$sdev^2)) * 100
    var_df <- tibble::tibble(
      substate      = as.character(ck),
      PC            = paste0("PC", seq_len(n_keep_pc)),
      var_explained = round(var_pct, 3),
      cum_var       = round(cumsum(var_pct), 3)
    )

    # Welch t on PC1_oriented across the two groups. Underpowered substates
    # (n < 3 per group) get NA p-values rather than t-test errors.
    pc1_niu   <- pc1_oriented[grp_chr == "NIU"]
    pc1_viral <- pc1_oriented[grp_chr == "Viral"]
    sig_row <- tibble::tibble(
      substate     = as.character(ck),
      n_NIU        = length(pc1_niu),
      n_Viral      = length(pc1_viral),
      mean_NIU     = if (length(pc1_niu)   > 0) mean(pc1_niu,   na.rm = TRUE) else NA_real_,
      mean_Viral   = if (length(pc1_viral) > 0) mean(pc1_viral, na.rm = TRUE) else NA_real_,
      t_statistic  = NA_real_,
      df           = NA_real_,
      p_value      = NA_real_,
      # Arm provenance. Carried on the significance table (not on scores) so
      # the published four-column readers in R/85 / R/86 / R/88 are unaffected
      # and every adjusted run is self-describing.
      covariate_mode = covariate_mode,
      covariate_col  = covariate_col %||% NA_character_,
      arm_status     = adj_status,
      arm_vif        = adj_vif,
      arm_n_excluded = adj_nex,
      # "vst" for every substate on the published path. Anything else means the
      # standard transform errored and a fallback was used, which only happens
      # for the sparse substates the completeness supplement exists to show.
      vst_method     = fit$vst_method %||% "vst",
      # Centroid-based sign flip actually applied. Under adjustment a substate
      # can lose its separation, at which point the flip is arbitrary; downstream
      # concordance aligns to the published PC1 by correlation instead.
      pc1_flip       = flip
    )
    if (length(pc1_niu) >= 3 && length(pc1_viral) >= 3) {
      tt <- tryCatch(stats::t.test(pc1_viral, pc1_niu, var.equal = FALSE),
                     error = function(e) NULL)
      if (!is.null(tt)) {
        sig_row$t_statistic <- unname(tt$statistic)
        sig_row$df          <- unname(tt$parameter)
        sig_row$p_value     <- tt$p.value
      }
    }

    scores_rows[[ck]]   <- scores_df
    loadings_rows[[ck]] <- rot_df
    variance_rows[[ck]] <- var_df
    sig_rows[[ck]]      <- sig_row
  }

  if (length(scores_rows) == 0) return(NULL)

  scores   <- dplyr::bind_rows(scores_rows)
  loadings <- dplyr::bind_rows(loadings_rows)
  variance <- dplyr::bind_rows(variance_rows)
  sig      <- dplyr::bind_rows(sig_rows)
  sig$q_value    <- stats::p.adjust(sig$p_value, method = "BH")
  sig$separating <- !is.na(sig$q_value) & sig$q_value < pc1_split_fdr

  list(scores = scores, loadings = loadings,
       variance = variance, significance = sig)
}

# Entry point called from run_pipeline.R Phase 1d. Reads the compartment
# Seurat object from disk, runs compute_per_substate_pca with the
# compartment-specific floor, writes four CSVs, then optionally calls
# .pc1_loadings_by_program when cfg$<target>_programs is defined.
run_compartment_pca <- function(cfg,
                                target = c("myeloid", "tcell", "bcell"),
                                covariate_col  = NULL,
                                covariate_mode = "none",
                                exclude_genes  = NULL,
                                out_dir        = NULL,
                                out_suffix     = "") {
  target <- match.arg(target)
  paths <- get_target_paths(cfg, target)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("compartment_pca[", target, "]: object missing at ",
                obj_path, "; skipping.")
    return(invisible(FALSE))
  }
  log_message("=== compartment_pca [", target, "] ===")
  obj <- readRDS(obj_path)

  min_cells <- .pca_min_cells(cfg, target)
  cpcfg     <- cfg$compartment_pca %||% list()
  res <- compute_per_substate_pca(
    obj,
    min_cells_per_pb = min_cells,
    min_gene_count   = as.integer(cpcfg$min_gene_count %||% 10L),
    hvg_n            = as.integer(cpcfg$hvg_n %||% 2000L),
    n_pcs            = as.integer(cpcfg$n_pcs %||% 5L),
    vst_blind        = isTRUE(cpcfg$vst_blind),
    pc1_split_fdr    = as.numeric(cpcfg$pc1_split_fdr %||% 0.05),
    covariate_col    = covariate_col,
    covariate_mode   = covariate_mode,
    exclude_genes    = exclude_genes,
    vif_abort        = as.numeric(
      (cfg$stress_sensitivity$collinearity$vif_abort) %||% 10)
  )
  if (is.null(res)) {
    log_message("compartment_pca[", target, "]: PCA returned nothing.")
    return(invisible(FALSE))
  }

  # out_dir / out_suffix let the sensitivity harness (R/46b) reuse this
  # orchestration while writing to outputs/tables/stress_sensitivity/ instead of
  # over the published CSVs. Defaults reproduce the published paths exactly.
  tab_dir <- out_dir %||% paths$results_tables
  ensure_dir(tab_dir)
  fp <- function(stem) file.path(tab_dir, paste0(stem, out_suffix, ".csv"))
  utils::write.csv(res$scores,       fp("pca_subject_scores"),      row.names = FALSE)
  utils::write.csv(res$loadings,     fp("pca_gene_loadings"),       row.names = FALSE)
  utils::write.csv(res$variance,     fp("pca_variance_explained"),  row.names = FALSE)
  utils::write.csv(res$significance, fp("pca_pc1_significance"),    row.names = FALSE)
  log_message("compartment_pca[", target, "]: wrote ",
              nrow(res$scores), " sample rows, ",
              nrow(res$loadings), " gene rows, ",
              sum(res$significance$separating, na.rm = TRUE),
              " separating substate(s).")

  # Per-program loadings (Figure 4 panel E). Skipped when the program list
  # is not defined in config — the bridge / LIANA / NicheNet steps don't
  # depend on it.
  programs <- cfg[[paste0(target, "_programs")]]
  if (!is.null(programs) && length(programs) > 0) {
    .pc1_loadings_by_program(target, programs, cfg,
                             scores = res$scores,
                             loadings = res$loadings,
                             significance = res$significance)
  }

  invisible(TRUE)
}

# Restrict PC1 loadings to a curated set of program gene panels. Used for
# F3 panel E (myeloid programs: HLA_I, HLA_II, CD1, ...) and F4 panel E
# (tcell programs: TCR_signal, Checkpoint, ...). One row per (substate,
# program, gene) keeping only separating substates and only genes present
# in the substate's loadings table.
.pc1_loadings_by_program <- function(target, programs, cfg,
                                     scores, loadings, significance) {
  paths <- get_target_paths(cfg, target)
  sep_subs <- significance$substate[significance$separating]
  if (length(sep_subs) == 0L) {
    log_message("  pc1_loadings_by_program[", target,
                "]: no separating substates; writing empty stub.")
    out <- tibble::tibble(substate = character(), program = character(),
                          gene = character(), PC1_oriented = double(),
                          loading_rank_within_program = integer())
    utils::write.csv(out,
                     file.path(paths$results_tables,
                               "pc1_loadings_by_program.csv"),
                     row.names = FALSE)
    return(invisible(out))
  }
  rows <- list()
  for (prog in names(programs)) {
    genes <- as.character(programs[[prog]])
    for (ck in sep_subs) {
      df <- dplyr::filter(loadings,
                          .data$substate == ck,
                          .data$gene %in% genes)
      if (nrow(df) == 0L) next
      df <- df |>
        dplyr::mutate(program = prog,
                      loading_rank_within_program =
                        rank(-abs(.data$PC1_oriented),
                             ties.method = "first")) |>
        dplyr::select(substate, program, gene, PC1_oriented,
                      loading_rank_within_program)
      rows[[paste(ck, prog, sep = "::")]] <- df
    }
  }
  out <- if (length(rows) == 0L) {
    tibble::tibble(substate = character(), program = character(),
                   gene = character(), PC1_oriented = double(),
                   loading_rank_within_program = integer())
  } else dplyr::bind_rows(rows)

  utils::write.csv(out,
                   file.path(paths$results_tables,
                             "pc1_loadings_by_program.csv"),
                   row.names = FALSE)
  log_message("  pc1_loadings_by_program[", target, "]: wrote ",
              nrow(out), " (substate, program, gene) rows across ",
              length(sep_subs), " separating substate(s).")
  invisible(out)
}
