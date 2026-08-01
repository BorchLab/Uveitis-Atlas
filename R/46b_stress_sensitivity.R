# R/46b_stress_sensitivity.R
#
# Evaluate-then-decide harness for Reviewer 1 Major 1. Recomputes the four core
# disease-axis results with the dissociation / acute-stress score removed, and
# scores the answer against criteria that were committed to config BEFORE the
# first run, so the outcome cannot be threshold-shopped afterwards.
#
# Targets:
#   eye_dge            eye pseudobulk DGE (Figure 2), design ~ stress + group
#   myeloid_pca        per-substate PC1 (Figure 3E)
#   tcell_pca          per-substate PC1 (Figure 4D)
#   bridge             myeloid -> T cell PC1 coupling (Figures 4E / 4F)
#   blood_replication  the falsification arm; see below
#
# Entry point: run_stress_sensitivity(cfg), gated on cfg$steps$stress_sensitivity.
# Prerequisite: run_stress_qc(cfg) (R/19) must have written
# outputs/tables/eye/stress_ucell_per_cell.csv and ..._per_pseudobulk.csv.
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# Published artifacts. The harness refuses to write any file with these
# basenames, whatever directory it is pointed at.
.PUBLISHED_PROTECTED <- c(
  "DGE_pseudobulk_Autoimmune_vs_Viral.csv",
  "pca_subject_scores.csv", "pca_gene_loadings.csv",
  "pca_variance_explained.csv", "pca_pc1_significance.csv",
  "pc1_loadings_by_program.csv",
  "pc1_bridge_scores.csv", "pc1_bridge_correlation.csv",
  "pc1_bridge_per_substate.csv")

.assert_not_published_path <- function(path) {
  if (basename(path) %in% .PUBLISHED_PROTECTED)
    stop("stress_sensitivity refuses to write a published artifact: ", path,
         call. = FALSE)
  invisible(path)
}

.sens_paths <- function(cfg) {
  p <- cfg$paths_stress_sensitivity %||% list(
    tables = "outputs/tables/stress_sensitivity",
    viz    = "outputs/viz/stress_sensitivity")
  p$snapshot <- file.path(p$tables, "published_snapshot")
  p
}

.sens_write <- function(x, dir, stem) {
  ensure_dir(dir)
  path <- file.path(dir, paste0(stem, ".csv"))
  .assert_not_published_path(path)
  utils::write.csv(x, path, row.names = FALSE)
  log_message("Wrote ", path)
  invisible(path)
}

# Freeze the baseline ONCE. Re-snapshotting mid-analysis would silently swap the
# comparison target if a pipeline step ran in between, so the guard is a hard
# skip rather than an overwrite.
.sens_snapshot <- function(cfg) {
  sp <- .sens_paths(cfg)$snapshot
  if (dir.exists(sp) && length(list.files(sp, recursive = TRUE))) {
    log_message("  snapshot already present at ", sp, "; leaving untouched.")
    return(invisible(sp))
  }
  ensure_dir(sp)
  srcs <- list(
    eye     = get_target_paths(cfg, "eye")$results_tables,
    myeloid = get_target_paths(cfg, "myeloid")$results_tables,
    tcell   = get_target_paths(cfg, "tcell")$results_tables,
    bcell   = get_target_paths(cfg, "bcell")$results_tables,
    cross   = (cfg$paths_cross_compartment %||%
                 list(tables = "outputs/tables/cross_compartment"))$tables)
  n <- 0L
  for (k in names(srcs)) {
    d <- file.path(sp, k); ensure_dir(d)
    for (f in .PUBLISHED_PROTECTED) {
      s <- file.path(srcs[[k]], f)
      if (file.exists(s)) { file.copy(s, file.path(d, f), overwrite = FALSE); n <- n + 1L }
    }
  }
  log_message("  snapshotted ", n, " published CSVs to ", sp)
  invisible(sp)
}

.snap_read <- function(cfg, which, file) {
  f <- file.path(.sens_paths(cfg)$snapshot, which, file)
  if (!file.exists(f)) return(NULL)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Concordance primitives
# ---------------------------------------------------------------------------

# AMENDMENT 2026-08-01, recorded rather than applied silently.
#
# The gating metrics are therefore now DIRECTIONAL: they fire only on
# degradation. Two changes, both narrowing what counts as MATERIAL:
#   * cohens_d_ratio (adjusted / published) replaces abs relative change, and
#     a ratio above 1 can never be MATERIAL.
#   * signed_delta_r replaces abs_delta_r on the bridge, so a coupling that
#     strengthens under adjustment cannot trip the gate.
# PC1 rotation (pc1_cor_to_published) is demoted to non-gating: an axis that
# rotates toward better disease separation is a finding, not a failure. It is
# still computed and reported.
#
.verdict <- function(value, minor, material, higher_is_better = TRUE) {
  if (!is.finite(value)) return("indeterminate")
  if (higher_is_better) {
    if (value >= minor) "MINOR" else if (value < material) "MATERIAL" else "EQUIVOCAL"
  } else {
    if (value <= minor) "MINOR" else if (value > material) "MATERIAL" else "EQUIVOCAL"
  }
}

.vrow <- function(target, arm, covariate, metric, value, minor, material,
                  higher_is_better = TRUE, gating = TRUE) {
  data.frame(target = target, arm = arm, covariate = covariate, metric = metric,
             value = value, threshold_minor = minor,
             threshold_material = material,
             verdict = .verdict(value, minor, material, higher_is_better),
             gating = gating, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Target: eye pseudobulk DGE (Figure 2)
# ---------------------------------------------------------------------------

# The audit that needs no re-run and makes no identifiability claim. Answers
# directly: are the intraocular DEGs actually dissociation genes? This is the
# sentence for the response letter.
.sens_gene_membership <- function(cfg, gene_sets) {
  pub <- .snap_read(cfg, "eye", "DGE_pseudobulk_Autoimmune_vs_Viral.csv")
  if (is.null(pub)) {
    log_message("  membership audit: published eye DGE not in the snapshot.")
    return(NULL)
  }
  padj_thr <- as.numeric(cfg$dge$padj_threshold %||% 0.05)
  lfc_thr  <- as.numeric(cfg$dge$lfc_threshold  %||% 0.5)
  rows <- list()
  for (cl in unique(pub$cluster)) {
    d <- pub[pub$cluster == cl, , drop = FALSE]
    sig <- !is.na(d$padj) & d$padj < padj_thr & abs(d$log2FoldChange) > lfc_thr
    if (sum(sig) < 5L) next
    for (nm in names(gene_sets)) {
      mem <- toupper(d$gene) %in% toupper(gene_sets[[nm]])
      a <- sum(sig & mem); b <- sum(sig & !mem)
      cc <- sum(!sig & mem); dd <- sum(!sig & !mem)
      ft <- tryCatch(stats::fisher.test(matrix(c(a, b, cc, dd), 2),
                                        alternative = "greater"),
                     error = function(e) NULL)
      rows[[length(rows) + 1]] <- data.frame(
        cluster = as.character(cl), gene_set = nm,
        n_tested = nrow(d), n_significant = sum(sig),
        n_sig_in_set = a, n_set_tested = sum(mem),
        frac_sig_in_set = a / max(1L, sum(sig)),
        frac_background_in_set = sum(mem) / max(1L, nrow(d)),
        odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
        fisher_p   = if (is.null(ft)) NA_real_ else ft$p.value,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- dplyr::bind_rows(rows)
  out$fisher_fdr <- stats::p.adjust(out$fisher_p, method = "BH")
  out
}

.sens_eye_dge <- function(cfg, pb_cov, covariate, sens_dir, gene_sets) {
  pub <- .snap_read(cfg, "eye", "DGE_pseudobulk_Autoimmune_vs_Viral.csv")
  if (is.null(pub)) return(NULL)
  paths_eye <- get_target_paths(cfg, "eye")
  op <- file.path(paths_eye$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(op)) {
    log_message("  eye_dge: eye object missing; skipping the adjusted re-run.")
    return(NULL)
  }
  cov_df <- pb_cov |>
    dplyr::filter(.data$target == "eye") |>
    dplyr::transmute(cluster = .data$substate, sample = .data$sample,
                     !!covariate := .data[[covariate]])
  log_message("  eye_dge: loading eye object for the stress-adjusted re-run")
  obj <- readRDS(op)
  cluster_col <- if ("knn.leiden.cluster" %in% colnames(obj[[]]))
                   "knn.leiden.cluster" else "seurat_clusters"
  adj <- tryCatch(
    run_pseudobulk_deseq2(obj, group_col = "Phenotype_2",
                          group1 = "NIU", group2 = "Viral",
                          cluster_col = cluster_col, cfg = cfg,
                          target = "eye", force_simple = TRUE,
                          covariate_df = cov_df, covariate_cols = covariate),
    error = function(e) { log_message("  eye_dge adjusted run failed: ",
                                      conditionMessage(e)); NULL })
  rm(obj); gc(verbose = FALSE)
  if (is.null(adj) || !nrow(adj)) return(NULL)
  .sens_write(adj, sens_dir,
              paste0("DGE_pseudobulk_Autoimmune_vs_Viral_stressadj_", covariate))

  crit <- cfg$stress_sensitivity$criteria$eye_dge %||% list()
  padj_thr <- as.numeric(cfg$dge$padj_threshold %||% 0.05)
  lfc_thr  <- as.numeric(cfg$dge$lfc_threshold  %||% 0.5)
  named    <- toupper(as.character(crit$named_genes %||% character(0)))
  vr <- list(); detail <- list()
  for (cl in intersect(unique(pub$cluster), unique(adj$cluster))) {
    p <- pub[pub$cluster == cl, ]; q <- adj[adj$cluster == cl, ]
    m <- dplyr::inner_join(
      dplyr::select(p, gene = "gene", lfc_pub = "log2FoldChange", padj_pub = "padj"),
      dplyr::select(q, gene = "gene", lfc_adj = "log2FoldChange", padj_adj = "padj"),
      by = "gene")
    if (nrow(m) < 50L) next
    sp <- suppressWarnings(stats::cor(m$lfc_pub, m$lfc_adj, method = "spearman"))
    sig_p <- !is.na(m$padj_pub) & m$padj_pub < padj_thr & abs(m$lfc_pub) > lfc_thr
    sig_a <- !is.na(m$padj_adj) & m$padj_adj < padj_thr & abs(m$lfc_adj) > lfc_thr
    jac <- if (sum(sig_p | sig_a)) sum(sig_p & sig_a) / sum(sig_p | sig_a) else NA_real_
    flip <- if (sum(sig_p)) mean(sign(m$lfc_pub[sig_p]) != sign(m$lfc_adj[sig_p])) else NA_real_
    top <- m[sig_p, ][order(m$padj_pub[sig_p]), ]
    top <- utils::head(top, 50L)
    ret <- if (nrow(top)) mean(!is.na(top$padj_adj) & top$padj_adj < padj_thr) else NA_real_
    nm_bad <- if (length(named)) {
      k <- toupper(m$gene) %in% named & sig_p
      if (any(k)) sum(sign(m$lfc_pub[k]) != sign(m$lfc_adj[k]) |
                        is.na(m$padj_adj[k]) | m$padj_adj[k] >= padj_thr) else 0L
    } else 0L
    detail[[length(detail) + 1]] <- data.frame(
      cluster = as.character(cl), n_common = nrow(m),
      n_sig_published = sum(sig_p), n_sig_adjusted = sum(sig_a),
      lfc_spearman = sp, deg_jaccard = jac, sign_flip_rate = flip,
      top50_retention = ret, n_named_genes_lost = nm_bad,
      stringsAsFactors = FALSE)
    if (!identical(as.character(cl), "global")) next
    vr[[length(vr) + 1]] <- .vrow("eye_dge", "stress_adjusted", covariate,
      "lfc_spearman_global", sp,
      crit$lfc_spearman_minor %||% 0.90, crit$lfc_spearman_material %||% 0.80)
    vr[[length(vr) + 1]] <- .vrow("eye_dge", "stress_adjusted", covariate,
      "deg_jaccard_global", jac,
      crit$deg_jaccard_minor %||% 0.70, crit$deg_jaccard_material %||% 0.50)
    vr[[length(vr) + 1]] <- .vrow("eye_dge", "stress_adjusted", covariate,
      "sign_flip_rate_global", flip,
      (crit$sign_flip_rate_material %||% 0.02) / 2,
      crit$sign_flip_rate_material %||% 0.02, higher_is_better = FALSE)
    vr[[length(vr) + 1]] <- .vrow("eye_dge", "stress_adjusted", covariate,
      "top50_retention_global", ret,
      crit$top50_retention_minor %||% 0.80, crit$top50_retention_material %||% 0.60)
    if (length(named))
      vr[[length(vr) + 1]] <- .vrow("eye_dge", "stress_adjusted", covariate,
        "named_genes_lost_global", nm_bad, 0, 0.5, higher_is_better = FALSE)
  }
  if (length(detail))
    .sens_write(dplyr::bind_rows(detail), sens_dir,
                paste0("eye_dge_concordance_", covariate))
  if (length(vr)) dplyr::bind_rows(vr) else NULL
}

# ---------------------------------------------------------------------------
# Targets: myeloid / T cell per-substate PC1 (Figures 3E, 4D)
# ---------------------------------------------------------------------------

.sens_pca <- function(cfg, target, covariate, arms, sens_dir, gene_sets) {
  pub_sc  <- .snap_read(cfg, target, "pca_subject_scores.csv")
  pub_sig <- .snap_read(cfg, target, "pca_pc1_significance.csv")
  pub_var <- .snap_read(cfg, target, "pca_variance_explained.csv")
  if (is.null(pub_sc) || is.null(pub_sig)) {
    log_message("  ", target, "_pca: published snapshot missing; skipping.")
    return(NULL)
  }
  p <- get_target_paths(cfg, target)
  op <- file.path(p$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(op)) {
    log_message("  ", target, "_pca: object missing; skipping.")
    return(NULL)
  }
  log_message("  ", target, "_pca: loading object and stamping stress scores")
  obj <- .stress_scores_for_object(readRDS(op), cfg)
  cpcfg <- cfg$compartment_pca %||% list()
  vif_abort <- as.numeric(cfg$stress_sensitivity$collinearity$vif_abort %||% 10)
  crit <- cfg$stress_sensitivity$criteria$pca %||% list()
  gating_arms <- as.character(cfg$stress_sensitivity$gating_arms %||%
                                "rbe_group_protected")
  # read.csv gives logical or the strings "TRUE"/"True"; %in% handles both and
  # never returns NA, unlike a bare == comparison.
  sep_pub <- as.character(
    pub_sig$substate[pub_sig$separating %in% c(TRUE, "TRUE", "True")])

  # Build the pseudobulks ONCE, carrying the covariate, and reuse across arms.
  pbs <- build_per_substate_pseudobulks(
    obj, min_cells_per_pb = .pca_min_cells(cfg, target),
    covariate_cols = covariate)
  vr <- list(); detail <- list()
  for (arm in setdiff(arms, "published")) {
    log_message("    arm: ", arm)
    res <- tryCatch(compute_per_substate_pca(
      obj, min_cells_per_pb = .pca_min_cells(cfg, target),
      min_gene_count = as.integer(cpcfg$min_gene_count %||% 10L),
      hvg_n          = as.integer(cpcfg$hvg_n %||% 2000L),
      n_pcs          = as.integer(cpcfg$n_pcs %||% 5L),
      vst_blind      = isTRUE(cpcfg$vst_blind),
      pc1_split_fdr  = as.numeric(cpcfg$pc1_split_fdr %||% 0.05),
      covariate_col  = covariate, covariate_mode = arm,
      exclude_genes  = unique(unlist(gene_sets)),
      vif_abort      = vif_abort, pbs = pbs),
      error = function(e) { log_message("    arm ", arm, " failed: ",
                                        conditionMessage(e)); NULL })
    if (is.null(res)) next
    sfx <- paste0("_", arm, "_", covariate)
    .sens_write(res$scores,       file.path(sens_dir, target),
                paste0("pca_subject_scores", sfx))
    .sens_write(res$significance, file.path(sens_dir, target),
                paste0("pca_pc1_significance", sfx))
    .sens_write(res$variance,     file.path(sens_dir, target),
                paste0("pca_variance_explained", sfx))

    cors <- c(); dds <- c()
    for (ss in unique(res$scores$substate)) {
      a <- res$scores[res$scores$substate == ss, ]
      b <- pub_sc[pub_sc$substate == ss, ]
      j <- dplyr::inner_join(
        dplyr::select(a, sample = "sample", adj = "PC1_oriented",
                      grp = "Phenotype_2"),
        dplyr::select(b, sample = "sample", pub = "PC1_oriented"), by = "sample")
      if (nrow(j) < 4L) next
      r <- suppressWarnings(stats::cor(j$adj, j$pub))
      # Sign-align to the published axis rather than trusting the centroid flip:
      # once a substate loses separation under adjustment the flip is arbitrary.
      cohd <- function(v, g) {
        x <- v[g == "NIU"]; y <- v[g == "Viral"]
        if (length(x) < 2 || length(y) < 2) return(NA_real_)
        sp <- sqrt(((length(x) - 1) * stats::var(x) +
                    (length(y) - 1) * stats::var(y)) / (length(x) + length(y) - 2))
        if (!is.finite(sp) || sp == 0) return(NA_real_)
        (mean(y) - mean(x)) / sp
      }
      d_adj <- cohd(j$adj, j$grp)
      d_pub <- cohd(j$pub, j$grp)
      va <- res$variance$var_explained[res$variance$substate == ss &
                                         res$variance$PC == "PC1"]
      vp <- if (!is.null(pub_var))
              pub_var$var_explained[pub_var$substate == ss & pub_var$PC == "PC1"]
            else numeric(0)
      sa <- res$significance[res$significance$substate == ss, ]
      sb <- pub_sig[pub_sig$substate == ss, ]
      was_sep <- nrow(sb) > 0 && isTRUE(sb$separating[1])
      cors <- c(cors, stats::setNames(abs(r), ss))
      # Directional: ratio of adjusted to published effect size. >= 1 means the
      # adjusted axis separates at least as well, and must never trip the gate.
      d_ratio <- if (is.finite(d_pub) && d_pub != 0 && is.finite(d_adj))
                   abs(d_adj) / abs(d_pub) else NA_real_
      if (is.finite(d_ratio)) dds <- c(dds, stats::setNames(d_ratio, ss))
      detail[[length(detail) + 1]] <- data.frame(
        target = target, arm = arm, covariate = covariate,
        substate = as.character(ss), n_shared = nrow(j),
        pc1_cor_to_published = r, pc1_sign_flipped = isTRUE(r < 0),
        cohens_d_published = d_pub, cohens_d_adjusted = d_adj,
        cohens_d_ratio = d_ratio,
        cohens_d_rel_change = if (is.finite(d_pub) && d_pub != 0)
                                abs((d_adj - d_pub) / d_pub) else NA_real_,
        var_explained_published = if (length(vp)) vp[1] else NA_real_,
        var_explained_adjusted  = if (length(va)) va[1] else NA_real_,
        separating_published = was_sep,
        separating_adjusted  = nrow(sa) > 0 && isTRUE(sa$separating[1]),
        arm_status = if (nrow(sa)) sa$arm_status[1] else NA_character_,
        arm_vif    = if (nrow(sa)) sa$arm_vif[1] else NA_real_,
        stringsAsFactors = FALSE)
    }
    gate <- arm %in% gating_arms
    sep_adj <- as.character(
      res$significance$substate[res$significance$separating %in% TRUE])
    jac <- if (length(union(sep_pub, sep_adj)))
             length(intersect(sep_pub, sep_adj)) / length(union(sep_pub, sep_adj))
           else NA_real_
    # Gate on the substates that separated in the PUBLISHED data. A substate that
    # was never separating cannot "lose" separation, and including it would
    # dilute the metric with noise from underpowered facets.
    cors_sep <- cors[names(cors) %in% sep_pub]
    dds_sep  <- dds[names(dds) %in% sep_pub]
    # --- NON-GATING, informational -----------------------------------------
    # PC1 rotation. Demoted from gating by the 2026-08-01 amendment: an axis
    # that rotates toward BETTER disease separation is a result, not a failure,
    # and this metric cannot tell the two apart. Still reported.
    vr[[length(vr) + 1]] <- .vrow(paste0(target, "_pca"), arm, covariate,
      "min_abs_pc1_cor_separating",
      if (length(cors_sep)) min(cors_sep, na.rm = TRUE) else NA_real_,
      crit$pc1_cor_minor %||% 0.85, crit$pc1_cor_material %||% 0.60,
      gating = FALSE)
    vr[[length(vr) + 1]] <- .vrow(paste0(target, "_pca"), arm, covariate,
      "median_abs_pc1_cor", if (length(cors)) stats::median(cors) else NA_real_,
      crit$pc1_cor_minor %||% 0.85, crit$pc1_cor_median_material %||% 0.70,
      gating = FALSE)
    # Jaccard penalises GAINING a separating substate as much as losing one, so
    # it is informational too; n_separating_lost below is the directional form.
    vr[[length(vr) + 1]] <- .vrow(paste0(target, "_pca"), arm, covariate,
      "separating_jaccard", jac,
      crit$separating_jaccard_minor %||% 0.80,
      crit$separating_jaccard_material %||% 0.60, gating = FALSE)
    vr[[length(vr) + 1]] <- .vrow(paste0(target, "_pca"), arm, covariate,
      "n_separating_gained",
      length(setdiff(sep_adj, sep_pub)), 0, Inf,
      higher_is_better = FALSE, gating = FALSE)

    # --- GATING, directional ------------------------------------------------
    # Did a published-separating substate LOSE separation?
    vr[[length(vr) + 1]] <- .vrow(paste0(target, "_pca"), arm, covariate,
      "n_separating_lost",
      length(setdiff(sep_pub, sep_adj)), 0, 0.5,
      higher_is_better = FALSE, gating = gate)
    # Did the effect size SHRINK? Ratio of adjusted to published |d|, minimum
    # across published-separating substates. >= 1 means it grew.
    vr[[length(vr) + 1]] <- .vrow(paste0(target, "_pca"), arm, covariate,
      "min_cohens_d_ratio",
      if (length(dds_sep)) min(dds_sep, na.rm = TRUE) else NA_real_,
      crit$cohens_d_ratio_minor %||% 0.75,
      crit$cohens_d_ratio_material %||% 0.50, gating = gate)
  }
  rm(obj, pbs); gc(verbose = FALSE)
  if (length(detail))
    .sens_write(dplyr::bind_rows(detail), sens_dir,
                paste0(target, "_pca_concordance_", covariate))
  if (length(vr)) dplyr::bind_rows(vr) else NULL
}

# ---------------------------------------------------------------------------
# Target: the myeloid -> T cell bridge (Figures 4E / 4F)
# ---------------------------------------------------------------------------

# Two conceptually different adjustments, both run and both labelled:
#   partial_cov  keep the published PC1 scores, add subject-level stress as an
#                extra continuous nuisance term. "Is the residual coupling just
#                between-subject stress variation?"
#   pca_<arm>    re-run the bridge on the adjusted PC1 scores with the published
#                partial model. "If the substate axes are de-stressed, does the
#                coupling survive?"
.sens_bridge <- function(cfg, subj_cov, covariate, arms, sens_dir) {
  crit <- cfg$stress_sensitivity$criteria$bridge %||% list()
  gating_arms <- c("partial_cov",
                   as.character(cfg$stress_sensitivity$gating_arms %||%
                                  "rbe_group_protected"))
  bdir <- file.path(sens_dir, "bridge")
  n_boot <- as.integer(cfg$stress_sensitivity$n_bootstrap %||% 2000L)
  n_perm <- as.integer(cfg$stress_sensitivity$n_permutation %||% 1000L)

  pub_cor <- .snap_read(cfg, "cross", "pc1_bridge_correlation.csv")
  pub_ps  <- .snap_read(cfg, "cross", "pc1_bridge_per_substate.csv")
  r_pub <- as.numeric(crit$published_r_unweighted %||% 0.700)
  ci_lo <- as.numeric(crit$published_ci_lo_unweighted %||% 0.439)
  ci_hi <- as.numeric(crit$published_ci_hi_unweighted %||% 0.841)

  vr <- list()
  score_of <- function(arm, target) {
    if (identical(arm, "partial_cov"))
      return(file.path(.sens_paths(cfg)$snapshot, target, "pca_subject_scores.csv"))
    f <- file.path(sens_dir, target,
                   paste0("pca_subject_scores_", arm, "_", covariate, ".csv"))
    if (file.exists(f)) f else NULL
  }
  for (arm in c("partial_cov", setdiff(arms, "published"))) {
    mf <- score_of(arm, "myeloid"); tf <- score_of(arm, "tcell")
    if (is.null(mf) || is.null(tf) || !file.exists(mf) || !file.exists(tf)) {
      log_message("    bridge arm ", arm, ": PC1 scores unavailable; skipping.")
      next
    }
    log_message("    bridge arm: ", arm)
    ok <- tryCatch(run_cross_compartment_pc1_bridge(
      cfg, n_boot = n_boot, n_perm = n_perm,
      covariate_df   = if (identical(arm, "partial_cov")) subj_cov else NULL,
      covariate_cols = if (identical(arm, "partial_cov")) covariate else NULL,
      myel_scores_csv = mf, tcel_scores_csv = tf,
      out_dir = bdir, out_suffix = paste0("_stressadj_", arm, "_", covariate)),
      error = function(e) { log_message("    bridge arm ", arm, " failed: ",
                                        conditionMessage(e)); FALSE })
    if (!isTRUE(ok)) next

    cf <- file.path(bdir, paste0("pc1_bridge_correlation_stressadj_", arm, "_",
                                 covariate, ".csv"))
    if (!file.exists(cf)) next
    adj <- utils::read.csv(cf, stringsAsFactors = FALSE)
    row <- adj[grepl("^partial_controlling_Phenotype_2", adj$stratum) &
                 adj$weighting == "unweighted", , drop = FALSE]
    if (!nrow(row)) next
    gate <- arm %in% gating_arms
    vr[[length(vr) + 1]] <- .vrow("bridge", arm, covariate,
      "partial_r_unweighted", row$pearson_r[1], ci_lo, ci_lo, gating = gate)
    vr[[length(vr) + 1]] <- .vrow("bridge", arm, covariate,
      # Signed, not absolute: a coupling that STRENGTHENS under adjustment must
      # not trip the gate. Reported as the drop from published, so negative
      # values mean the adjusted r is larger. (2026-08-01 amendment.)
      "delta_r_drop_vs_published", r_pub - row$pearson_r[1],
      crit$delta_r_minor %||% 0.10, crit$delta_r_material %||% 0.20,
      higher_is_better = FALSE, gating = gate)
    vr[[length(vr) + 1]] <- .vrow("bridge", arm, covariate,
      "permutation_p", row$permutation_p[1],
      (crit$perm_p_material %||% 0.05) / 2, crit$perm_p_material %||% 0.05,
      higher_is_better = FALSE, gating = gate)
    roww <- adj[grepl("^partial_controlling_Phenotype_2", adj$stratum) &
                  adj$weighting == "weighted", , drop = FALSE]
    if (nrow(roww))
      vr[[length(vr) + 1]] <- .vrow("bridge", arm, covariate,
        "partial_r_weighted", roww$pearson_r[1],
        crit$published_ci_lo_weighted %||% 0.346,
        crit$published_ci_lo_weighted %||% 0.346, gating = gate)

    # F4F: does the per-substate-pair coupling map hold its shape?
    pf <- file.path(bdir, paste0("pc1_bridge_per_substate_stressadj_", arm, "_",
                                 covariate, ".csv"))
    if (!is.null(pub_ps) && file.exists(pf)) {
      ps <- utils::read.csv(pf, stringsAsFactors = FALSE)
      j <- dplyr::inner_join(
        dplyr::select(pub_ps, "myeloid_substate", "tcell_substate",
                      r_pub = "partial_pearson_r"),
        dplyr::select(ps, "myeloid_substate", "tcell_substate",
                      r_adj = "partial_pearson_r"),
        by = c("myeloid_substate", "tcell_substate"))
      if (nrow(j) >= 5L) {
        rho <- suppressWarnings(stats::cor(j$r_pub, j$r_adj, method = "spearman"))
        top <- j[order(-abs(j$r_pub)), ][seq_len(max(1L, floor(nrow(j) / 10))), ]
        vr[[length(vr) + 1]] <- .vrow("bridge", arm, covariate,
          "per_substate_r_spearman", rho,
          crit$per_substate_spearman_minor %||% 0.85,
          crit$per_substate_spearman_material %||% 0.70, gating = gate)
        vr[[length(vr) + 1]] <- .vrow("bridge", arm, covariate,
          "per_substate_top_decile_sign_flips",
          sum(sign(top$r_pub) != sign(top$r_adj), na.rm = TRUE),
          0, 0.5, higher_is_better = FALSE, gating = gate)
        .sens_write(j, bdir, paste0("per_substate_concordance_", arm, "_", covariate))
      }
    }
  }
  if (length(vr)) dplyr::bind_rows(vr) else NULL
}

# ---------------------------------------------------------------------------
# Target: blood replication — the falsification arm
# ---------------------------------------------------------------------------

# Blood is never passed through the 35 um strainer, in either arm. If the eye
# NIU-vs-Viral signature is recovered directionally in blood, the filtration
# explanation is contradicted along an axis where filtration does not vary. This
# is the only target here that can falsify rather than bound.
#
.sens_blood_replication <- function(cfg, sens_dir) {
  op <- file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(op)) {
    log_message("  blood_replication: full atlas object missing; skipping.")
    return(NULL)
  }
  pub <- .snap_read(cfg, "eye", "DGE_pseudobulk_Autoimmune_vs_Viral.csv")
  log_message("  blood_replication: loading the full atlas object (~3.5 GB)")
  obj <- readRDS(op)
  meta <- obj[[]]
  ctb <- resolve_celltype_broad(meta)
  if (is.null(ctb)) { rm(obj); gc(verbose = FALSE); return(NULL) }

  run_one <- function(cells, gcol, g1, g2, stem) {
    if (length(cells) < 100L) return(NULL)
    o <- subset(obj, cells = cells)
    o$.pb_cluster <- as.character(o[[]][[ctb]])
    r <- tryCatch(run_pseudobulk_deseq2(o, group_col = gcol, group1 = g1,
                                        group2 = g2, cluster_col = ".pb_cluster",
                                        cfg = cfg, target = "all",
                                        force_simple = TRUE),
                  error = function(e) { log_message("    ", stem, " failed: ",
                                                    conditionMessage(e)); NULL })
    rm(o); gc(verbose = FALSE)
    if (is.null(r) || !nrow(r)) return(NULL)
    r$stratum <- stem
    .sens_write(r, sens_dir, stem)
    r
  }

  is_blood <- meta$Tissue_1 %in% "Blood"
  bl <- run_one(rownames(meta)[is_blood & meta$Phenotype_2 %in% c("NIU", "Viral")],
                "Phenotype_2", "NIU", "Viral",
                "DGE_pseudobulk_Autoimmune_vs_Viral_blood_only")
  # Healthy blood, US vs Japan: etiology held at zero, so this is a pure
  # site/cohort estimator. n = 3 vs 4 — underpowered for discovery, and the
  # caveat belongs in the same sentence as any number taken from it.
  hc <- run_one(rownames(meta)[is_blood &
                                 meta$Phenotype_2 %in% c("Healthy", "HC")],
                "Cohort", "US", "Japan",
                "DGE_pseudobulk_Cohort_healthy_blood_only")
  rm(obj); gc(verbose = FALSE)

  concord <- function(x, label) {
    if (is.null(x) || is.null(pub)) return(NULL)
    a <- pub[pub$cluster == "global", ]; b <- x[x$cluster == "global", ]
    if (!nrow(a) || !nrow(b)) return(NULL)
    m <- dplyr::inner_join(
      dplyr::select(a, gene = "gene", lfc_eye = "log2FoldChange", padj_eye = "padj"),
      dplyr::select(b, gene = "gene", lfc_alt = "log2FoldChange", padj_alt = "padj"),
      by = "gene")
    if (nrow(m) < 50L) return(NULL)
    sig <- !is.na(m$padj_eye) & m$padj_eye < 0.05
    agree <- if (sum(sig)) mean(sign(m$lfc_eye[sig]) == sign(m$lfc_alt[sig])) else NA_real_
    st <- if (sum(sig) >= 5)
      stats::binom.test(sum(sign(m$lfc_eye[sig]) == sign(m$lfc_alt[sig])),
                        sum(sig), 0.5, alternative = "greater")$p.value else NA_real_
    data.frame(comparison = label, n_common = nrow(m),
               n_eye_significant = sum(sig),
               lfc_spearman = suppressWarnings(
                 stats::cor(m$lfc_eye, m$lfc_alt, method = "spearman")),
               sign_agreement_rate = agree, sign_test_p = st,
               stringsAsFactors = FALSE)
  }
  out <- dplyr::bind_rows(
    concord(bl, "eye_NIUvsViral_vs_blood_NIUvsViral"),
    concord(hc, "eye_NIUvsViral_vs_healthyblood_USvsJapan"))
  if (!is.null(out) && nrow(out))
    .sens_write(out, sens_dir, "blood_replication_concordance")
  out
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

run_stress_sensitivity <- function(cfg) {
  scfg <- cfg$stress_sensitivity %||% list()
  sp <- .sens_paths(cfg)
  sens_dir <- sp$tables
  ensure_dir(sens_dir)
  log_message("=== stress_sensitivity (Reviewer 1 Major 1) ===")
  .sens_snapshot(cfg)

  paths_eye <- get_target_paths(cfg, "eye")
  f_pc <- file.path(paths_eye$results_tables, "stress_ucell_per_cell.csv")
  f_pb <- file.path(paths_eye$results_tables, "stress_ucell_per_pseudobulk.csv")
  f_sj <- file.path(paths_eye$results_tables, "stress_ucell_per_subject.csv")
  f_co <- file.path(paths_eye$results_tables, "stress_ucell_collinearity.csv")
  if (!file.exists(f_pc)) {
    log_message("stress_sensitivity: ", f_pc,
                " missing. Run steps$stress_qc first.")
    return(invisible(FALSE))
  }
  gene_sets <- .stress_gene_sets(cfg)
  covariate <- scfg$primary_covariate %||% "stress_dissoc_pruned_ucell"
  arms      <- as.character(scfg$pca_arms %||%
                              c("rbe_group_protected", "rbe_naive", "hvg_drop"))
  targets   <- as.character(scfg$targets %||%
                              c("eye_dge", "myeloid_pca", "tcell_pca", "bridge"))

  # The membership audit runs first and unconditionally: it is the one number
  # that is completely unaffected by the identifiability problem.
  mem <- .sens_gene_membership(cfg, gene_sets)
  if (!is.null(mem)) {
    .sens_write(mem, sens_dir, "dge_stress_gene_membership")
    g <- mem[mem$cluster == "global" & mem$gene_set == "stress_dissoc", ]
    if (nrow(g))
      log_message(sprintf(
        "  MEMBERSHIP: %d of %d significant intraocular DEGs (%.1f%%) are van den Brink dissociation genes, vs %.1f%% of the tested background (OR = %.2f, Fisher p = %.3g).",
        g$n_sig_in_set[1], g$n_significant[1], 100 * g$frac_sig_in_set[1],
        100 * g$frac_background_in_set[1], g$odds_ratio[1], g$fisher_p[1]))
  }

  # Identifiability gate. If the covariate is near-aliased with the group term
  # everywhere, saying so is a STRONGER position than quoting an unstable
  # coefficient, so the adjusted arms are skipped rather than fudged.
  collin <- if (file.exists(f_co)) utils::read.csv(f_co, stringsAsFactors = FALSE) else NULL
  if (!is.null(collin)) {
    cp <- collin[collin$score == covariate & collin$level == "subject", ]
    if (nrow(cp) && is.finite(cp$vif[1]) &&
        cp$vif[1] > as.numeric(scfg$collinearity$vif_abort %||% 10)) {
      log_message(sprintf(
        "  NON-IDENTIFIABLE: subject-level VIF for %s is %.1f, above the abort threshold. Skipping all adjusted arms; report the membership audit and the blood replication instead.",
        covariate, cp$vif[1]))
      targets <- intersect(targets, "blood_replication")
    }
  }

  pb_cov <- if (file.exists(f_pb))
    utils::read.csv(f_pb, stringsAsFactors = FALSE) else NULL
  subj_cov <- NULL
  if (file.exists(f_sj)) {
    s <- utils::read.csv(f_sj, stringsAsFactors = FALSE)
    s <- s[s$level == "Subject", , drop = FALSE]
    subj_cov <- dplyr::transmute(s, subject = .data$Subject,
                                 !!covariate := .data[[covariate]])
  }

  vr <- list()
  if ("eye_dge" %in% targets && !is.null(pb_cov))
    vr[["eye_dge"]] <- .sens_eye_dge(cfg, pb_cov, covariate, sens_dir, gene_sets)
  for (tg in c("myeloid", "tcell"))
    if (paste0(tg, "_pca") %in% targets)
      vr[[tg]] <- .sens_pca(cfg, tg, covariate, arms, sens_dir, gene_sets)
  if ("bridge" %in% targets && !is.null(subj_cov))
    vr[["bridge"]] <- .sens_bridge(cfg, subj_cov, covariate, arms, sens_dir)
  if ("blood_replication" %in% targets)
    .sens_blood_replication(cfg, sens_dir)

  verdict <- dplyr::bind_rows(Filter(Negate(is.null), vr))
  if (!nrow(verdict)) {
    log_message("stress_sensitivity: no verdict rows produced.")
    return(invisible(FALSE))
  }
  .sens_write(verdict, sens_dir, "stress_sensitivity_verdict")

  # Global gate. Re-run the pipeline iff any GATING arm returns MATERIAL. F4E is
  # the mechanistic hinge of Figure 4, so it escalates on EQUIVOCAL too; the
  # script deliberately does not decide that case on its own.
  g <- verdict[verdict$gating %in% TRUE, , drop = FALSE]
  n_mat <- sum(g$verdict == "MATERIAL")
  n_equ <- sum(g$verdict == "EQUIVOCAL")
  bridge_equ <- any(g$target == "bridge" & g$verdict == "EQUIVOCAL")
  overall <- if (n_mat > 0L) "MATERIAL -> RE-RUN THE PIPELINE"
             else if (bridge_equ || n_equ >= 2L) "EQUIVOCAL -> HUMAN DECISION"
             else "MINOR -> supplemental figure only, published numbers stand"
  log_message("  ---------------------------------------------------------")
  log_message("  GLOBAL GATE: ", overall)
  log_message(sprintf("  gating rows: %d MATERIAL, %d EQUIVOCAL, %d MINOR",
                      n_mat, n_equ, sum(g$verdict == "MINOR")))
  if (n_mat > 0L)
    for (i in which(g$verdict == "MATERIAL"))
      log_message(sprintf("    MATERIAL: %s / %s / %s = %.4g",
                          g$target[i], g$arm[i], g$metric[i], g$value[i]))
  log_message("  ---------------------------------------------------------")
  invisible(list(verdict = verdict, overall = overall, membership = mem))
}
