# R/19_stress_qc.R
#
# Per-cell acute-stress / dissociation scoring for the eye sub-atlas, plus the
# identifiability diagnostic that governs whether any stress-adjusted contrast
# in R/46b is estimable at all.
#
# Naming: occupies slot 19, the name the revision plan asked for. Structurally
# it follows
# R/15_lens_qc.R (UCell scoring, config-driven panels, read-only w.r.t. the
# Seurat objects, re-runnable at no cost).
#
# Entry points:
#   * run_stress_qc(cfg) — score every eye cell, write the per-cell / per-
#     subject / per-pseudobulk score tables, run the subject-level NIU vs Viral
#     test, emit the collinearity diagnostic, render the QC figures.
#   * .stress_scores_for_object(obj, cfg) — stamp the eye-derived per-cell
#     scores onto a compartment object by barcode. Used by R/46b so that F2, F3
#     and F4 are all adjusted on literally the same covariate.
#
# Outputs (all under outputs/):
#   tables/eye/stress_ucell_per_cell.csv
#   tables/eye/stress_ucell_gene_coverage.csv
#   tables/eye/stress_ucell_per_subject.csv
#   tables/eye/stress_ucell_per_pseudobulk.csv
#   tables/eye/stress_ucell_group_test.csv
#   tables/eye/stress_ucell_by_celltype.csv
#   tables/eye/stress_ucell_by_etiology.csv
#   tables/eye/stress_ucell_collinearity.csv
#   viz/eye/01_qc/stress_ucell_{distribution,subjects,etiology,depth,celltype,identifiability}.{pdf,png}
#
# Metadata assumptions (set in R/02 + R/03): Subject, Subject_Timepoint,
# Phenotype_2 (NIU | Viral), Tissue_1, Tissue_2, Cohort, Etiology,
# knn.leiden.cluster, celltype / celltype_broad (possibly _full suffixed on the
# eye object — always resolve via resolve_celltype_broad()). Missing columns are
# tolerated with NA fallbacks. HLA_B27 lives in inputs/data/metadata.csv and is
# NOT propagated into the objects; join on Subject if you need it.
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------

# van den Brink et al., Nat Methods 2017, "Single-cell sequencing reveals
# dissociation-induced gene expression in tissue subpopulations", Table S1.
# The published set is mouse; this is the human ortholog set. Most entries are a
# plain case change. The ones that are not:
#   Ccrn4l  -> NOCT       (gene renamed)
#   Cyr61   -> CCN1       (renamed; CYR61 also listed for older references)
#   Fam132b -> ERFE       (renamed)
#   H3f3b   -> H3-3B      (GRCh38-2020-A) or H3F3B (older refs) — both listed,
#                          intersect() picks whichever the reference carries
#   Mt1/Mt2 -> no 1:1 human ortholog; expanded to the human MT1 family + MT2A
# .compute_stress_ucell() intersects against rownames(obj) and writes the
# realized coverage to stress_ucell_gene_coverage.csv, so any deviation from the
# published panel is auditable rather than silent.
.STRESS_DISSOC_VDB <- c(
  "ACTG1","ANKRD1","ARID5A","ATF3","ATF4","BAG3","BHLHE40","BRD2","BTG1",
  "BTG2","CCNL1","NOCT","CEBPB","CEBPD","CEBPG","CSRNP1","CXCL1","CCN1",
  "CYR61","DCN","DDX3X","DDX5","DES","DNAJA1","DNAJB1","DNAJB4","DUSP1",
  "DUSP8","EGR1","EGR2","EIF1","EIF5","ERF","ERRFI1","ERFE","FOS","FOSB",
  "FOSL2","GADD45B","GADD45G","GCC1","GEM","H3-3B","H3F3B","HIPK3",
  "HSP90AA1","HSP90AB1","HSPA1A","HSPA1B","HSPA5","HSPA8","HSPB1","HSPE1",
  "HSPH1","ID3","IDI1","IER2","IER3","IER5","IL6","IRF1","IRF8","ITPKC",
  "JUN","JUNB","JUND","KLF2","KLF4","KLF6","KLF9","LITAF","LMNA","MAFF",
  "MAFK","MCL1","MIDN","MIR22HG","MT1X","MT1E","MT1F","MT1G","MT2A","MYADM",
  "MYC","MYD88","NCKAP5L","NCOA7","NFKBIA","NFKBIZ","NOP58","NPPC","NR4A1",
  "ODC1","OSGIN1","OXNAD1","PCF11","PDE4B","PER1","PHLDA1","PNP","PNRC1",
  "PPP1CC","PPP1R15A","PXDC1","RAP1B","RASSF1","RHOB","RHOH","RIPK1","SAT1",
  "SBNO2","SDC4","SERPINE1","SKIL","SLC10A6","SLC38A2","SLC41A1","SOCS3",
  "SQSTM1","SRF","SRSF5","SRSF7","STAT3","TAGLN2","TIPARP","TNFAIP3",
  "TNFAIP6","TPM3","TPPP3","TRA2A","TRA2B","TRIB1","TUBB4B","TUBB6","UBC",
  "USP2","WAC","ZC3H12A","ZFAND5","ZFP36","ZFP36L1","ZFP36L2","ZYX")

# Subtracted from the above to form `dissociation_pruned`.
#
# This is the single most consequential analytic choice in the module, so state
# it plainly: the van den Brink panel is not a clean artifact list. Roughly
# twenty of its members are canonical NF-kB and inflammatory-response effectors
# that are expected to differ between granulomatous autoimmune uveitis and acute
# necrotizing viral retinitis for real biological reasons. Using the unpruned
# set as an ADJUSTMENT covariate does not remove an artifact, it removes the
# finding.
#
# Therefore: the FULL set is what we report descriptively (it is the set the
# reviewer named, and the descriptive comparison should use it unmodified). The
# PRUNED set is the primary adjustment covariate. The full-set adjustment is
# also run, and reported as a deliberately over-corrected worst-case bound.
.STRESS_IMMUNE_PRUNE <- c(
  "IL6","CXCL1","IRF1","IRF8","SOCS3","NFKBIA","NFKBIZ","TNFAIP3","TNFAIP6",
  "STAT3","MYD88","RIPK1","ZC3H12A","SERPINE1","CEBPB","CEBPD","MCL1",
  "SDC4","LITAF","SBNO2","TRIB1","PDE4B")

# Immediate-early genes, as specified in the revision plan.
#
# Caveat that belongs in the manuscript, not just here: FOS, FOSB, JUN, JUNB and
# EGR1 are the immediate downstream of TCR and Fc-receptor ligation in vivo, not
# only of ex vivo handling. The IEG-adjusted arm is therefore over-conservative
# for the T cell compartment specifically — it strips activation signal along
# with any handling signal. That cuts in the paper's favour (if F4D survives IEG
# adjustment it has survived a deliberately hostile correction), so report it,
# but never make it primary.
.STRESS_IEG_CORE <- c("FOS","FOSB","JUN","JUNB","EGR1",
                      "HSPA1A","HSPA1B","HSPB1","DNAJB1")

# Reviewer 1's mechanism (1), scored rather than argued.
#
# The reviewer's first proposed mechanism is biological, not technical: the
# vitreous sits behind the blood-retinal barrier, exchanges material slowly, and
# is claimed to hold more tissue-resident memory cells and retina-derived glia
# than the rapidly turned-over aqueous compartment. That is a testable claim,
# and unlike the filtration argument it is not something an adjustment can
# remove -- if it is true, it is part of the biology the paper is describing.
#
# Direction is deliberately NOT pre-committed. The manuscript currently reports
# the ITGA1+ effector memory cluster as showing only mild, condition-specific
# shifts, so there is no published TRM enrichment to defend and the contrast is
# reported whichever way it falls.
#
# Naming caution for the manuscript: "TRM" is currently overloaded in the text,
# meaning tissue-resident macrophage in the myeloid section and tissue-resident
# memory T cell here. Both usages need disambiguating.
.TRM_CORE <- c("ITGA1", "ITGAE", "CXCR6", "ZNF683", "CD69")

# Retinal / glial contamination panel. The reviewer explicitly raises glial
# migration from the adjacent retina into the vitreous. GFAP, SLC1A3 and AQP4
# are astrocyte / Muller glia; RLBP1 is Muller glia and RPE; RPE65 and PMEL are
# retinal pigment epithelium. If these are present at non-trivial frequency in
# the viral (vitreous) arm, that is a real compositional difference between the
# compartments and belongs in the limitation paragraph as biology rather than
# being silently absorbed into a "stress" score.
.GLIAL_CORE <- c("GFAP", "RLBP1", "AQP4", "SLC1A3", "PMEL", "RPE65")

.stress_gene_sets <- function(cfg) {
  g   <- cfg$stress_qc$gene_sets %||% list()
  vdb <- unique(toupper(as.character(g$dissociation_vdb %||% .STRESS_DISSOC_VDB)))
  prn <- unique(toupper(as.character(g$dissociation_immune_prune %||%
                                       .STRESS_IMMUNE_PRUNE)))
  ieg <- unique(toupper(as.character(g$ieg_core %||% .STRESS_IEG_CORE)))
  trm <- unique(toupper(as.character(g$trm_core %||% .TRM_CORE)))
  gli <- unique(toupper(as.character(g$glial_core %||% .GLIAL_CORE)))
  list(stress_dissoc        = vdb,
       stress_dissoc_pruned = setdiff(vdb, prn),
       stress_ieg           = ieg,
       trm                  = trm,
       glial                = gli)
}

# Only the dissociation / IEG panels are stress covariates. TRM and glial are
# scored on the same pass (the UCell rank matrix is computed once) but describe
# compartment biology, so they are reported and never used as adjustments.
.STRESS_SCORE_COLS <- c("stress_dissoc_ucell", "stress_dissoc_pruned_ucell",
                        "stress_ieg_ucell", "stress_combined_z",
                        "trm_ucell", "glial_ucell")
.STRESS_COVARIATE_COLS <- c("stress_dissoc_ucell", "stress_dissoc_pruned_ucell",
                            "stress_ieg_ucell", "stress_combined_z")

# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

# All sets go into ONE AddModuleScore_UCell call so the per-cell rank matrix is
# computed once. maxRank is passed explicitly (R/15 relies on the 1500 default)
# and recorded in the coverage table: vitreous cells may be systematically
# shallower, and a score saturating against maxRank would be a *different*
# artifact from the one the reviewer named. The depth panel exists to check it.
.compute_stress_ucell <- function(obj, gene_sets, seed = 42L,
                                  maxrank = 1500L, min_genes_present = 3L,
                                  min_frac_present = 0.5) {
  if (!requireNamespace("UCell", quietly = TRUE))
    stop("R/19_stress_qc: UCell required. BiocManager::install('UCell').")

  present <- lapply(gene_sets, function(g) intersect(g, rownames(obj)))
  coverage <- do.call(rbind, lapply(names(gene_sets), function(nm) {
    data.frame(set = nm,
               n_requested = length(gene_sets[[nm]]),
               n_present   = length(present[[nm]]),
               frac_present = round(length(present[[nm]]) /
                                      max(1L, length(gene_sets[[nm]])), 4),
               maxRank = as.integer(maxrank),
               missing = paste(setdiff(gene_sets[[nm]], present[[nm]]),
                               collapse = ";"),
               stringsAsFactors = FALSE)
  }))
  for (nm in names(present))
    log_message(sprintf("  stress set '%s': %d / %d genes present", nm,
                        length(present[[nm]]), length(gene_sets[[nm]])))

  # Coverage gate. Two conditions, deliberately: an absolute floor so a panel
  # reduced to a handful of genes is not scored, and a FRACTION floor so a panel
  # that mostly survived is kept. A pure absolute floor would silently drop
  # ieg_core, which is 9 genes by design and fully present — exactly the panel
  # the reviewer asked for.
  n_req <- vapply(gene_sets[names(present)], length, integer(1))
  n_pre <- vapply(present, length, integer(1))
  keep <- names(present)[n_pre >= min_genes_present &
                           (n_pre / pmax(1L, n_req)) >= min_frac_present]
  if (length(keep) == 0L)
    stop("R/19_stress_qc: no gene set met the coverage gate (>= ",
         min_genes_present, " genes and >= ", 100 * min_frac_present,
         "% of the requested panel).")
  dropped <- setdiff(names(present), keep)
  if (length(dropped))
    log_message("  WARN: dropping under-covered set(s): ",
                paste(dropped, collapse = ", "))

  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  set_global_seed(seed)
  scored <- UCell::AddModuleScore_UCell(obj, features = present[keep],
                                        name = "_ucell", assay = "RNA",
                                        maxRank = as.integer(maxrank))
  sm <- scored[[]][, paste0(keep, "_ucell"), drop = FALSE]
  scores <- data.frame(cell_id = colnames(scored), sm,
                       stringsAsFactors = FALSE, check.names = FALSE)

  # Combined score: mean of the z-scored pruned-dissociation and IEG scores,
  # z-scored across all eye cells. Deliberately NOT a UCell score and named
  # differently so nobody treats it as one.
  zz <- function(v) { s <- stats::sd(v, na.rm = TRUE)
                      if (!is.finite(s) || s == 0) rep(0, length(v))
                      else (v - mean(v, na.rm = TRUE)) / s }
  a <- if ("stress_dissoc_pruned_ucell" %in% names(scores))
         zz(scores$stress_dissoc_pruned_ucell) else NULL
  b <- if ("stress_ieg_ucell" %in% names(scores))
         zz(scores$stress_ieg_ucell) else NULL
  scores$stress_combined_z <- if (!is.null(a) && !is.null(b)) (a + b) / 2
                              else if (!is.null(a)) a else b

  list(scores = scores, coverage = coverage)
}

# Stamp the eye-derived per-cell scores onto any compartment object by barcode.
#
# Scoring once on the eye object and propagating is deliberate. UCell scores are
# per-cell ranks over the object's gene universe, and the compartment objects
# were re-integrated with their own feature selection. Re-scoring per compartment
# would give subtly different values for the same cell, and the claim that F2, F3
# and F4 were adjusted on "the same covariate" would be false. The hard match
# assert exists because a barcode rename during re-integration would otherwise
# turn into a silent all-NA join.
.stress_scores_for_object <- function(obj, cfg, score_cols = NULL,
                                      min_match = 0.99) {
  paths_eye <- get_target_paths(cfg, "eye")
  f <- file.path(paths_eye$results_tables, "stress_ucell_per_cell.csv")
  if (!file.exists(f))
    stop("R/19_stress_qc: ", f, " missing. Run stress_qc first.")
  pc <- utils::read.csv(f, stringsAsFactors = FALSE)
  score_cols <- intersect(score_cols %||% .STRESS_SCORE_COLS, colnames(pc))
  idx <- match(colnames(obj), pc$cell_id)
  frac <- mean(!is.na(idx))
  if (frac < min_match)
    stop(sprintf(paste("R/19_stress_qc: only %.1f%% of cells matched the eye",
                       "score table (need >= %.0f%%). Barcodes likely diverged",
                       "between the eye object and this compartment."),
                 100 * frac, 100 * min_match))
  for (cn in score_cols) obj[[cn]] <- pc[[cn]][idx]
  obj
}

# ---------------------------------------------------------------------------
# Subject-level summaries and tests
# ---------------------------------------------------------------------------

# Two levels: Subject_Timepoint (the pseudobulk unit) and Subject (the
# independence unit, and the one the tests use). Timepoints collapse into
# subjects by a cell-count-weighted mean so a subject with two visits does not
# get double weight.
.stress_subject_summary <- function(per_cell, score_cols) {
  wmean <- function(v, w) {
    ok <- is.finite(v) & is.finite(w) & w > 0
    if (!any(ok)) return(NA_real_)
    stats::weighted.mean(v[ok], w[ok])
  }
  tp <- per_cell |>
    dplyr::group_by(.data$Subject_Timepoint, .data$Subject, .data$Phenotype_2,
                    .data$Tissue_2, .data$Cohort, .data$Etiology) |>
    dplyr::summarise(n_cells = dplyr::n(),
                     nFeature_RNA_mean = mean(.data$nFeature_RNA, na.rm = TRUE),
                     dplyr::across(dplyr::all_of(score_cols),
                                   \(x) mean(x, na.rm = TRUE)),
                     .groups = "drop")
  sj <- tp |>
    dplyr::group_by(.data$Subject, .data$Phenotype_2, .data$Tissue_2,
                    .data$Cohort, .data$Etiology) |>
    dplyr::summarise(n_cells = sum(.data$n_cells),
                     nFeature_RNA_mean = wmean(.data$nFeature_RNA_mean,
                                               .data$n_cells),
                     dplyr::across(dplyr::all_of(score_cols),
                                   \(x) wmean(x, .data$n_cells)),
                     .groups = "drop")
  dplyr::bind_rows(
    dplyr::mutate(tp, level = "Subject_Timepoint",
                  id = .data$Subject_Timepoint) |>
      dplyr::select(-"Subject_Timepoint"),
    dplyr::mutate(sj, level = "Subject", id = .data$Subject))
}

# Primary test. Deliberately not a t-test: UCell scores are bounded on [0, 1]
# and right-skewed, and n = 21 vs 13 is small enough that the distributional
# assumption would be doing real work. AUC = P(Viral > NIU) is the headline
# because a reviewer can interpret it without trusting any assumption.
.stress_group_test <- function(subject_df, score_cols,
                               group_col = "Phenotype_2",
                               groups = c("NIU", "Viral"),
                               n_boot = 5000L, seed = 42L) {
  d <- subject_df[subject_df$level == "Subject" &
                    subject_df[[group_col]] %in% groups, , drop = FALSE]
  auc_of <- function(a, b) {   # P(b > a), ties at 0.5
    if (!length(a) || !length(b)) return(NA_real_)
    mean(outer(b, a, ">") + 0.5 * outer(b, a, "=="))
  }
  hedges_g <- function(a, b) {
    na <- length(a); nb <- length(b)
    if (na < 2 || nb < 2) return(NA_real_)
    sp <- sqrt(((na - 1) * stats::var(a) + (nb - 1) * stats::var(b)) /
                 (na + nb - 2))
    if (!is.finite(sp) || sp == 0) return(NA_real_)
    d <- (mean(b) - mean(a)) / sp
    d * (1 - 3 / (4 * (na + nb) - 9))          # small-sample correction
  }
  set_global_seed(seed)
  rows <- lapply(score_cols, function(sc) {
    a <- d[[sc]][d[[group_col]] == groups[1]]; a <- a[is.finite(a)]
    b <- d[[sc]][d[[group_col]] == groups[2]]; b <- b[is.finite(b)]
    auc <- auc_of(a, b); g <- hedges_g(a, b)
    bs <- if (length(a) >= 3 && length(b) >= 3 && n_boot >= 100L)
      vapply(seq_len(n_boot), function(i) {
        aa <- sample(a, length(a), TRUE); bb <- sample(b, length(b), TRUE)
        c(auc_of(aa, bb), hedges_g(aa, bb))
      }, numeric(2)) else NULL
    qs <- function(v) if (is.null(bs)) c(NA_real_, NA_real_) else
      unname(stats::quantile(v[is.finite(v)], c(0.025, 0.975), na.rm = TRUE))
    wq <- if (length(a) >= 3 && length(b) >= 3)
      suppressWarnings(stats::wilcox.test(b, a, exact = TRUE)$p.value) else NA_real_
    data.frame(score = sc, group_ref = groups[1], group_alt = groups[2],
               n_ref = length(a), n_alt = length(b),
               mean_ref = mean(a), mean_alt = mean(b),
               median_ref = stats::median(a), median_alt = stats::median(b),
               auc_alt_gt_ref = auc,
               auc_ci_lo = qs(if (is.null(bs)) NA else bs[1, ])[1],
               auc_ci_hi = qs(if (is.null(bs)) NA else bs[1, ])[2],
               hedges_g = g,
               hedges_g_ci_lo = qs(if (is.null(bs)) NA else bs[2, ])[1],
               hedges_g_ci_hi = qs(if (is.null(bs)) NA else bs[2, ])[2],
               wilcox_p = wq,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  # Pre-specified premise verdict. This is a check on the REVIEWER'S factual
  # claim, not on our result, so it is stated before any adjustment is run.
  out$premise <- dplyr::case_when(
    is.na(out$auc_alt_gt_ref) ~ "indeterminate",
    out$auc_alt_gt_ref >= 0.35 & out$auc_alt_gt_ref <= 0.65 &
      abs(out$hedges_g) < 0.5 ~ "refuted",
    out$auc_alt_gt_ref > 0.85 | abs(out$hedges_g) > 1.5 ~ "confirmed_large",
    TRUE ~ "confirmed_moderate")
  out
}

# Secondary test. The random intercept on Subject is what stops ~10^5 cells from
# manufacturing p ~ 1e-300 by pseudoreplication. celltype_broad is not optional:
# vitreous is myeloid-shifted, so a marginal per-cell comparison confounds
# composition with per-cell stress. nFeature_RNA guards the UCell maxRank/depth
# interaction.
.stress_cell_lmm <- function(per_cell, score_cols, group_col = "Phenotype_2") {
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    log_message("  lmerTest not installed; skipping the mixed-model secondary.")
    return(NULL)
  }
  d <- per_cell[per_cell[[group_col]] %in% c("NIU", "Viral"), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$.grp <- factor(as.character(d[[group_col]]), levels = c("NIU", "Viral"))
  d$.ct  <- factor(ifelse(is.na(d$celltype_broad), "NA", d$celltype_broad))
  d$.nf  <- as.numeric(scale(d$nFeature_RNA))
  do.call(rbind, lapply(score_cols, function(sc) {
    d$.y <- d[[sc]]
    fit <- tryCatch(
      lmerTest::lmer(.y ~ .grp + .ct + .nf + (1 | Subject), data = d),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    co <- summary(fit)$coefficients
    rn <- grep("^\\.grpViral$", rownames(co), value = TRUE)
    icc <- tryCatch({
      n0 <- lme4::lmer(.y ~ (1 | Subject), data = d)
      vc <- as.data.frame(lme4::VarCorr(n0))
      vc$vcov[vc$grp == "Subject"] / sum(vc$vcov)
    }, error = function(e) NA_real_)
    data.frame(score = sc,
               term = "Phenotype_2 Viral vs NIU",
               estimate = if (length(rn)) co[rn, "Estimate"] else NA_real_,
               std_error = if (length(rn)) co[rn, "Std. Error"] else NA_real_,
               df = if (length(rn)) co[rn, "df"] else NA_real_,
               p_value = if (length(rn)) co[rn, "Pr(>|t|)"] else NA_real_,
               icc_subject = icc,
               singular = lme4::isSingular(fit),
               n_cells = nrow(d), n_subjects = length(unique(d$Subject)),
               stringsAsFactors = FALSE)
  }))
}

# The gate. Everything R/46b does rests on the stress score carrying usable
# within-group variance. Reported at both levels where the covariate is actually
# used: per subject, and per pseudobulk column within each compartment substate.
.stress_collinearity <- function(subject_df, pb_df, score_cols,
                                 group_col = "Phenotype_2",
                                 vif_warn = 5, vif_abort = 10) {
  one <- function(v, g, level, target, substate) {
    ok <- is.finite(v) & !is.na(g)
    v <- v[ok]; g <- as.character(g[ok])
    if (length(v) < 4L || length(unique(g)) < 2L) return(NULL)
    r2 <- suppressWarnings(summary(stats::lm(v ~ factor(g)))$r.squared)
    vif <- if (!is.finite(r2)) NA_real_ else if (r2 >= 1 - 1e-12) Inf else 1 / (1 - r2)
    mu <- tapply(v, g, mean)
    sw <- sqrt(mean(tapply(v, g, stats::var)[is.finite(tapply(v, g, stats::var))],
                    na.rm = TRUE))
    gap <- if (length(mu) >= 2) abs(diff(range(mu))) else NA_real_
    data.frame(level = level, target = target, substate = substate,
               n = length(v), r2_on_group = r2, vif = vif,
               within_group_sd = sw, between_group_gap = gap,
               sd_over_gap = if (is.finite(gap) && gap > 0) sw / gap else NA_real_,
               flag = if (!is.finite(vif)) "abort"
                      else if (vif > vif_abort) "abort"
                      else if (vif > vif_warn) "warn" else "ok",
               stringsAsFactors = FALSE)
  }
  # one() returns NULL for strata too small or single-group to characterise.
  # cbind(score = sc, NULL) would collapse to a 1-column frame and blow up the
  # rbind, so the NULL check has to happen BEFORE the cbind, not after.
  push <- function(rows, sc, res) {
    if (is.null(res)) return(rows)
    rows[[length(rows) + 1]] <- cbind(score = sc, res, stringsAsFactors = FALSE)
    rows
  }
  rows <- list()
  sj <- subject_df[subject_df$level == "Subject", , drop = FALSE]
  for (sc in intersect(score_cols, colnames(sj)))
    rows <- push(rows, sc,
                 one(sj[[sc]], sj[[group_col]], "subject", "eye", NA_character_))
  if (!is.null(pb_df) && nrow(pb_df)) {
    for (sc in intersect(score_cols, colnames(pb_df))) {
      for (tg in unique(pb_df$target)) {
        for (ss in unique(pb_df$substate[pb_df$target == tg])) {
          k <- pb_df$target == tg & pb_df$substate == ss
          rows <- push(rows, sc,
                       one(pb_df[[sc]][k], pb_df$Phenotype_2[k],
                           "pseudobulk", tg, as.character(ss)))
        }
      }
    }
  }
  if (!length(rows)) return(NULL)
  out <- dplyr::bind_rows(rows)
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Compartment biology: tissue residency and retinal glial carryover
# ---------------------------------------------------------------------------

# These answer Reviewer 1's mechanism (1), which is biological rather than
# technical and therefore is NOT something to adjust away. Two differences from
# the stress scoring above:
#
#   * The TRM contrast is restricted to the T cell compartment. Scoring ITGA1 /
#     ITGAE / CXCR6 / ZNF683 / CD69 across all eye cells would dilute it with
#     myeloid and B cells that do not express the panel, and the subject-level
#     mean would then be driven by composition rather than residency.
#   * Glial carryover is reported as a FRACTION of cells above a detection
#     threshold, not as a mean score. "What proportion of the vitreous cells are
#     retina-derived?" is the question the reviewer is actually asking, and a
#     mean UCell score over a population that is >99% leukocyte does not answer
#     it. The threshold is a high quantile of the panel score across all eye
#     cells, so it is defined relative to this dataset rather than imported.
.stress_compartment_biology <- function(per_cell, cfg, q = 0.99) {
  out <- list()
  ctb <- per_cell$celltype_broad

  if ("trm_ucell" %in% colnames(per_cell)) {
    tcells <- per_cell[!is.na(ctb) & ctb %in% c("T cell", "T/NK"), , drop = FALSE]
    if (nrow(tcells) > 100L) {
      out$trm <- tcells |>
        dplyr::group_by(.data$Subject, .data$Phenotype_2, .data$Tissue_2,
                        .data$Etiology) |>
        dplyr::summarise(n_tcells = dplyr::n(),
                         trm_mean = mean(.data$trm_ucell, na.rm = TRUE),
                         trm_median = stats::median(.data$trm_ucell, na.rm = TRUE),
                         .groups = "drop")
    }
  }

  if ("glial_ucell" %in% colnames(per_cell)) {
    thr <- stats::quantile(per_cell$glial_ucell, q, na.rm = TRUE)
    per_cell$.glial_pos <- per_cell$glial_ucell > thr
    out$glial <- per_cell |>
      dplyr::group_by(.data$Subject, .data$Phenotype_2, .data$Tissue_2,
                      .data$Etiology) |>
      dplyr::summarise(n_cells = dplyr::n(),
                       n_glial_pos = sum(.data$.glial_pos, na.rm = TRUE),
                       frac_glial_pos = mean(.data$.glial_pos, na.rm = TRUE),
                       glial_mean = mean(.data$glial_ucell, na.rm = TRUE),
                       .groups = "drop")
    out$glial$threshold <- unname(thr)
    out$glial$threshold_quantile <- q
  }
  out
}

# ---------------------------------------------------------------------------
# Per-pseudobulk covariate table (the load-bearing artifact for R/46b)
# ---------------------------------------------------------------------------

# One row per (target, substate, sample) with the mean score over exactly the
# cells that form that pseudobulk column. Computed straight off metadata rather
# than via build_per_substate_pseudobulks(): the means are identical and this
# skips a Seurat::AggregateExpression pass per compartment.
#
# Note the unit. A per-SAMPLE lookup would be wrong — the pseudobulk unit is
# (substate x sample), and the same subject can sit at very different stress
# levels in the myeloid and T cell compartments.
.stress_per_pseudobulk <- function(cfg, eye_meta, score_cols) {
  rows <- list()
  # cluster_col is explicit rather than resolved by name: the eye frame passed in
  # here is the per-cell table, whose cluster column is called `cluster`, while
  # the compartment objects carry `knn.leiden.cluster`. Resolving by name only
  # would silently return NULL for the eye target and leave the F2 DGE with no
  # covariate at all.
  add <- function(meta, target, floor, cluster_col = NULL, sample_col = NULL) {
    cluster_col <- cluster_col %||%
      (if ("knn.leiden.cluster" %in% colnames(meta)) "knn.leiden.cluster"
       else "seurat_clusters")
    # Same precedence as build_per_substate_pseudobulks() and
    # run_pseudobulk_deseq2(). On this cohort NO object carries
    # Subject_Timepoint, so both resolve to orig.ident.
    #
    # The eye frame passed in here is the per-cell TABLE, which does carry a
    # Subject_Timepoint column — but only as a synthesized fallback filled with
    # Subject. Resolving against it would key the covariate on subject IDs while
    # the DGE keys on orig.ident, and every join would miss. Callers therefore
    # pass sample_col explicitly for the eye frame.
    sample_col <- sample_col %||%
      (if ("Subject_Timepoint" %in% colnames(meta)) "Subject_Timepoint"
       else if ("orig.ident" %in% colnames(meta)) "orig.ident"
       else "Subject")
    need <- c(cluster_col, sample_col, "Phenotype_2", "Subject")
    if (!all(need %in% colnames(meta))) {
      log_message("  per-pseudobulk: ", target, " missing ",
                  paste(setdiff(need, colnames(meta)), collapse = "/"),
                  "; skipping.")
      return(NULL)
    }
    meta$.sub <- as.character(meta[[cluster_col]])
    meta$.smp <- as.character(meta[[sample_col]])
    base <- meta |>
      dplyr::filter(.data$Phenotype_2 %in% c("NIU", "Viral"),
                    !is.na(.data$.sub), !is.na(.data$.smp))
    if (!nrow(base)) return(NULL)
    per_sub <- base |>
      dplyr::group_by(substate = .data$.sub, sample = .data$.smp,
                      subject = .data$Subject, Phenotype_2 = .data$Phenotype_2) |>
      dplyr::summarise(n_cells = dplyr::n(),
                       nFeature_RNA_mean = mean(.data$nFeature_RNA, na.rm = TRUE),
                       dplyr::across(dplyr::all_of(score_cols),
                                     \(x) mean(x, na.rm = TRUE)),
                       .groups = "drop") |>
      dplyr::filter(.data$n_cells >= floor)
    # A "global" row per sample so the eye DGE (which loops c("global",
    # clusters)) can key on the same (cluster, sample) pair.
    glb <- base |>
      dplyr::group_by(sample = .data$.smp, subject = .data$Subject,
                      Phenotype_2 = .data$Phenotype_2) |>
      dplyr::summarise(n_cells = dplyr::n(),
                       nFeature_RNA_mean = mean(.data$nFeature_RNA, na.rm = TRUE),
                       dplyr::across(dplyr::all_of(score_cols),
                                     \(x) mean(x, na.rm = TRUE)),
                       .groups = "drop") |>
      dplyr::mutate(substate = "global")
    dplyr::bind_rows(per_sub, glb) |>
      dplyr::mutate(target = target, min_cells_per_pb = floor)
  }

  rows[["eye"]] <- add(eye_meta, "eye", 1L, cluster_col = "cluster",
                       sample_col = "orig.ident")
  for (tg in c("myeloid", "tcell", "bcell")) {
    p <- get_target_paths(cfg, tg)
    f <- file.path(p$results_objects, "IntegratedSeuratObject.rds")
    if (!file.exists(f)) {
      log_message("  per-pseudobulk: ", tg, " object missing; skipping.")
      next
    }
    log_message("  per-pseudobulk: loading ", tg)
    o <- readRDS(f)
    m <- o[[]]
    idx <- match(rownames(m), eye_meta$cell_id)
    if (mean(!is.na(idx)) < 0.99) {
      log_message("  WARN: only ", round(100 * mean(!is.na(idx)), 1),
                  "% of ", tg, " barcodes found in the eye score table; ",
                  "skipping this compartment rather than emitting NA means.")
      rm(o, m); gc(verbose = FALSE); next
    }
    for (sc in score_cols) m[[sc]] <- eye_meta[[sc]][idx]
    if (!"nFeature_RNA" %in% colnames(m)) m$nFeature_RNA <- NA_real_
    rows[[tg]] <- add(m, tg, .pca_min_cells(cfg, tg))
    rm(o, m); gc(verbose = FALSE)
  }
  out <- dplyr::bind_rows(Filter(Negate(is.null), rows))
  if (!nrow(out)) return(out)
  dplyr::select(out, "target", "substate", "sample", "subject", "Phenotype_2",
                "n_cells", "min_cells_per_pb", "nFeature_RNA_mean",
                dplyr::all_of(score_cols))
}

# ---------------------------------------------------------------------------
# Figures — one file per panel (no patchwork composites)
# ---------------------------------------------------------------------------

.stress_plots <- function(per_cell, subject_df, tests, collin, paths, cfg) {
  primary <- cfg$stress_qc$primary_score %||% "stress_dissoc_pruned_ucell"
  if (!primary %in% colnames(per_cell)) return(invisible(FALSE))
  vd <- viz_subdir(paths, "qc"); ensure_dir(vd)
  pal <- tryCatch(ETIOLOGY_GROUP_COLORS, error = function(e) NULL) %||%
           c(NIU = "#E21F26", Viral = "#397FB9")
  sj <- subject_df[subject_df$level == "Subject", , drop = FALSE]
  tt <- tests[tests$score == primary, , drop = FALSE]
  sub_lab <- if (nrow(tt))
    sprintf("AUC P(Viral>NIU) = %.2f [%.2f, %.2f]   Hedges g = %.2f   Wilcoxon p = %.3g   n = %d vs %d",
            tt$auc_alt_gt_ref[1], tt$auc_ci_lo[1], tt$auc_ci_hi[1],
            tt$hedges_g[1], tt$wilcox_p[1], tt$n_ref[1], tt$n_alt[1]) else ""

  # A. per-cell distribution by sampling site
  save_pdf_png(
    ggplot(per_cell, aes(.data[[primary]], fill = .data$Phenotype_2)) +
      geom_density(alpha = 0.55, colour = NA) +
      facet_wrap(~ Tissue_2, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = pal, name = NULL) +
      labs(title = "Per-cell dissociation score by sampling site",
           subtitle = paste0("Score: ", primary,
                             ". Aqueous is NIU-only and vitreous Viral-only ",
                             "by design, so the facets are not independent."),
           x = primary, y = "density") +
      theme_bw(base_size = 10),
    file.path(vd, "stress_ucell_distribution"), w = 7.5, h = 6)

  # B. subject level — the test that is actually reported
  save_pdf_png(
    ggplot(sj, aes(.data$Phenotype_2, .data[[primary]],
                   colour = .data$Phenotype_2)) +
      geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40",
                   fill = NA) +
      geom_jitter(width = 0.12, height = 0, size = 2.6, alpha = 0.9) +
      scale_colour_manual(values = pal, guide = "none") +
      labs(title = "Subject-level dissociation score, NIU vs Viral",
           subtitle = sub_lab, x = NULL, y = primary) +
      theme_bw(base_size = 10),
    file.path(vd, "stress_ucell_subjects"), w = 5.5, h = 5)

  # C. within-group negative control. If NIU sub-diagnoses spread as widely as
  # the NIU/Viral gap, the score is tracking biology, not just handling.
  save_pdf_png(
    ggplot(sj, aes(stats::reorder(.data$Etiology, .data[[primary]],
                                  FUN = stats::median),
                   .data[[primary]], colour = .data$Phenotype_2)) +
      geom_point(size = 2.6, alpha = 0.9) +
      scale_colour_manual(values = pal, name = NULL) +
      coord_flip() +
      labs(title = "Dissociation score by clinical etiology (subject level)",
           subtitle = "Within-NIU spread comparable to the NIU/Viral gap argues the score is not purely technical",
           x = NULL, y = primary) +
      theme_bw(base_size = 10),
    file.path(vd, "stress_ucell_etiology"), w = 7, h = 5.5)

  # D. depth confound. UCell scores rank genes up to maxRank, so a systematic
  # nFeature_RNA difference between sites could masquerade as stress.
  save_pdf_png(
    ggplot(per_cell, aes(.data$nFeature_RNA, .data[[primary]])) +
      geom_hex(bins = 60) +
      scale_fill_viridis_c(trans = "log10", name = "cells") +
      facet_wrap(~ Tissue_2, ncol = 2) +
      labs(title = "Dissociation score versus sequencing depth",
           subtitle = "Checks whether the site difference is partly a depth artifact",
           x = "nFeature_RNA", y = primary) +
      theme_bw(base_size = 10),
    file.path(vd, "stress_ucell_depth"), w = 9, h = 4.5)

  # E. is the stress global or cell-type-specific?
  ct <- per_cell |>
    dplyr::filter(!is.na(.data$celltype_broad)) |>
    dplyr::group_by(.data$celltype_broad, .data$Phenotype_2) |>
    dplyr::summarise(mean_score = mean(.data[[primary]], na.rm = TRUE),
                     n_cells = dplyr::n(), .groups = "drop")
  save_pdf_png(
    ggplot(ct, aes(.data$Phenotype_2, .data$celltype_broad,
                   fill = .data$mean_score)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.3f", .data$mean_score)), size = 3) +
      scale_fill_viridis_c(name = primary) +
      labs(title = "Mean dissociation score by broad cell type",
           x = NULL, y = NULL) +
      theme_bw(base_size = 10),
    file.path(vd, "stress_ucell_celltype"), w = 6.5, h = 5)

  # F. the identifiability panel. This is the one that decides whether any
  # adjusted number in R/46b is worth reporting.
  if (!is.null(collin) && nrow(collin)) {
    cc <- collin[collin$score == primary, , drop = FALSE]
    cc$lab <- ifelse(cc$level == "subject", "eye (subject)",
                     paste0(cc$target, " ", cc$substate, " (pb)"))
    save_pdf_png(
      ggplot(cc, aes(stats::reorder(.data$lab, .data$vif), .data$vif,
                     fill = .data$flag)) +
        geom_col() +
        geom_hline(yintercept = c(5, 10), linetype = "dashed",
                   colour = "grey40") +
        scale_fill_manual(values = c(ok = "#397FB9", warn = "#F0A202",
                                     abort = "#E21F26"), name = NULL) +
        coord_flip() +
        labs(title = "Is the stress covariate identifiable?",
             subtitle = paste("VIF of the stress score against Phenotype_2.",
                              "Above 10 the covariate is near-aliased with the",
                              "group term and R/46b refuses to adjust on it."),
             x = NULL, y = "variance inflation factor") +
        theme_bw(base_size = 10),
      file.path(vd, "stress_ucell_identifiability"),
      w = 8, h = max(4, 0.22 * nrow(cc) + 2))
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


# reuse_scores = TRUE skips the UCell pass and reads stress_ucell_per_cell.csv
# back instead. Scoring the full eye object takes ~9 minutes and is entirely
# deterministic given the seed and the gene sets, so iterating on the downstream
# summaries should not re-pay it. Any change to the gene sets, the seed or
# maxRank invalidates the cache and it must be re-run with reuse_scores = FALSE.
run_stress_qc <- function(cfg, reuse_scores = FALSE) {
  scfg <- cfg$stress_qc %||% list()
  paths_eye <- get_target_paths(cfg, "eye")
  eye_path <- file.path(paths_eye$results_objects, "IntegratedSeuratObject.rds")
  cache <- file.path(paths_eye$results_tables, "stress_ucell_per_cell.csv")
  if (isTRUE(reuse_scores) && file.exists(cache)) {
    log_message("=== stress_qc (reusing cached per-cell scores) ===")
    return(.run_stress_qc_downstream(
      cfg, utils::read.csv(cache, stringsAsFactors = FALSE), NULL))
  }
  if (!file.exists(eye_path)) {
    log_message("stress_qc: eye object not found at ", eye_path, ". Skipping.")
    return(invisible(FALSE))
  }
  log_message("=== stress_qc (dissociation / acute-stress signature) ===")
  eye <- readRDS(eye_path)

  gene_sets <- .stress_gene_sets(cfg)
  sc <- .compute_stress_ucell(eye, gene_sets,
                              seed = as.integer(scfg$seed %||% cfg$seed %||% 42L),
                              maxrank = as.integer(scfg$ucell_maxrank %||% 1500L),
                              min_genes_present = as.integer(
                                scfg$min_genes_present %||% 3L),
                              min_frac_present = as.numeric(
                                scfg$min_frac_present %||% 0.5))
  score_cols <- intersect(.STRESS_SCORE_COLS, colnames(sc$scores))

  meta <- eye[[]]
  ct_col  <- resolve_celltype(meta)
  ctb_col <- resolve_celltype_broad(meta)
  pick <- function(col, default = NA_character_) {
    if (!is.null(col) && col %in% colnames(meta)) as.character(meta[[col]])
    else rep(default, nrow(meta))
  }
  per_cell <- data.frame(
    cell_id           = colnames(eye),
    Subject           = pick("Subject"),
    Subject_Timepoint = if ("Subject_Timepoint" %in% colnames(meta))
                          as.character(meta$Subject_Timepoint) else pick("Subject"),
    # orig.ident matters: neither the eye object nor the compartment objects
    # carry Subject_Timepoint, so build_per_substate_pseudobulks() and
    # run_pseudobulk_deseq2() both fall through to orig.ident as the pseudobulk
    # sample key. The covariate table has to key on the same thing or the join
    # silently matches nothing.
    orig.ident        = pick("orig.ident"),
    Phenotype_2       = pick("Phenotype_2"),
    Tissue_1          = pick("Tissue_1"),
    Tissue_2          = pick("Tissue_2"),
    Cohort            = pick("Cohort"),
    Etiology          = pick("Etiology"),
    cluster           = pick(if ("knn.leiden.cluster" %in% colnames(meta))
                               "knn.leiden.cluster" else "seurat_clusters"),
    celltype          = pick(ct_col),
    celltype_broad    = pick(ctb_col),
    nCount_RNA        = if ("nCount_RNA" %in% colnames(meta))
                          as.numeric(meta$nCount_RNA) else NA_real_,
    nFeature_RNA      = if ("nFeature_RNA" %in% colnames(meta))
                          as.numeric(meta$nFeature_RNA) else NA_real_,
    stringsAsFactors  = FALSE)
  stopifnot(identical(per_cell$cell_id, sc$scores$cell_id))
  per_cell <- cbind(per_cell, sc$scores[, score_cols, drop = FALSE])
  rm(eye); gc(verbose = FALSE)

  ensure_dir(paths_eye$results_tables)
  utils::write.csv(per_cell,
                   file.path(paths_eye$results_tables, "stress_ucell_per_cell.csv"),
                   row.names = FALSE)
  log_message("Wrote ", file.path(paths_eye$results_tables,
                                  "stress_ucell_per_cell.csv"))
  .run_stress_qc_downstream(cfg, per_cell, sc$coverage)
}

# Everything after the UCell pass. Split out so run_stress_qc(reuse_scores=TRUE)
# can re-derive the summaries from the cached per-cell table.
.run_stress_qc_downstream <- function(cfg, per_cell, coverage) {
  scfg <- cfg$stress_qc %||% list()
  paths_eye <- get_target_paths(cfg, "eye")
  score_cols <- intersect(.STRESS_SCORE_COLS, colnames(per_cell))
  # Barcodes are written as "<orig.ident>.<10x barcode>" by R/02_ingest_data.R,
  # and that prefix is exactly the sample key the pseudobulk builders use
  # (verified against pca_subject_scores.csv$sample). Deriving it keeps an older
  # cached per-cell table usable and gives a fallback if the column is absent.
  if (!"orig.ident" %in% colnames(per_cell) || all(is.na(per_cell$orig.ident))) {
    per_cell$orig.ident <- sub("[.].*$", "", per_cell$cell_id)
    log_message("  orig.ident derived from the cell-barcode prefix (",
                length(unique(per_cell$orig.ident)), " samples).")
  }
  ensure_dir(paths_eye$results_tables)
  wr <- function(x, stem) {
    p <- file.path(paths_eye$results_tables, paste0(stem, ".csv"))
    utils::write.csv(x, p, row.names = FALSE); log_message("Wrote ", p)
  }
  if (!is.null(coverage)) wr(coverage, "stress_ucell_gene_coverage")

  subject_df <- .stress_subject_summary(per_cell, score_cols)
  wr(subject_df, "stress_ucell_per_subject")

  tests <- .stress_group_test(
    subject_df, score_cols,
    n_boot = as.integer(scfg$n_boot %||% 5000L),
    seed   = as.integer(scfg$seed %||% cfg$seed %||% 42L))
  lmm <- .stress_cell_lmm(per_cell, score_cols)
  if (!is.null(lmm)) {
    tests <- dplyr::left_join(
      tests,
      dplyr::select(lmm, "score", lmm_estimate = "estimate",
                    lmm_p = "p_value", lmm_icc_subject = "icc_subject",
                    lmm_singular = "singular"),
      by = "score")
  }
  wr(tests, "stress_ucell_group_test")

  by_ct <- per_cell |>
    dplyr::group_by(.data$celltype_broad, .data$cluster, .data$Phenotype_2) |>
    dplyr::summarise(n_cells = dplyr::n(),
                     dplyr::across(dplyr::all_of(score_cols),
                                   list(mean = \(x) mean(x, na.rm = TRUE),
                                        median = \(x) stats::median(x, na.rm = TRUE)),
                                   .names = "{.col}_{.fn}"),
                     .groups = "drop")
  wr(by_ct, "stress_ucell_by_celltype")

  by_eti <- subject_df |>
    dplyr::filter(.data$level == "Subject") |>
    dplyr::group_by(.data$Phenotype_2, .data$Etiology) |>
    dplyr::summarise(n_subjects = dplyr::n(),
                     dplyr::across(dplyr::all_of(score_cols),
                                   list(mean = \(x) mean(x, na.rm = TRUE),
                                        sd = \(x) stats::sd(x, na.rm = TRUE)),
                                   .names = "{.col}_{.fn}"),
                     .groups = "drop")
  wr(by_eti, "stress_ucell_by_etiology")

  # Reviewer 1 mechanism (1): compartment biology, reported not adjusted.
  bio <- .stress_compartment_biology(per_cell, cfg)
  if (!is.null(bio$trm)) {
    wr(bio$trm, "stress_trm_per_subject")
    tt <- .stress_group_test(
      dplyr::mutate(bio$trm, level = "Subject", id = .data$Subject,
                    trm_ucell = .data$trm_mean),
      "trm_ucell", n_boot = as.integer(scfg$n_boot %||% 5000L))
    wr(tt, "stress_trm_group_test")
    log_message(sprintf(
      "  TISSUE RESIDENCY (T cells only): AUC P(Viral>NIU) = %.3f [%.3f, %.3f], Hedges g = %.2f, Wilcoxon p = %.3g",
      tt$auc_alt_gt_ref[1], tt$auc_ci_lo[1], tt$auc_ci_hi[1],
      tt$hedges_g[1], tt$wilcox_p[1]))
  }
  if (!is.null(bio$glial)) {
    wr(bio$glial, "stress_glial_per_subject")
    g <- bio$glial
    log_message(sprintf(
      "  GLIAL CARRYOVER (fraction above the %.0fth percentile of the panel): NIU %.3f%% vs Viral %.3f%% of cells",
      100 * g$threshold_quantile[1],
      100 * mean(g$frac_glial_pos[g$Phenotype_2 == "NIU"], na.rm = TRUE),
      100 * mean(g$frac_glial_pos[g$Phenotype_2 == "Viral"], na.rm = TRUE)))
  }

  pb <- .stress_per_pseudobulk(cfg, per_cell, score_cols)
  if (nrow(pb)) wr(pb, "stress_ucell_per_pseudobulk")

  vcfg <- cfg$stress_sensitivity$collinearity %||% list()
  collin <- .stress_collinearity(
    subject_df, pb, score_cols,
    vif_warn  = as.numeric(vcfg$vif_warn %||% 5),
    vif_abort = as.numeric(vcfg$vif_abort %||% 10))
  if (!is.null(collin)) wr(collin, "stress_ucell_collinearity")

  tryCatch(.stress_plots(per_cell, subject_df, tests, collin, paths_eye, cfg),
           error = function(e)
             log_message("  stress_qc figures failed: ", conditionMessage(e)))

  # Console summary. These two numbers decide whether R/46b is worth running.
  primary <- scfg$primary_score %||% "stress_dissoc_pruned_ucell"
  tp <- tests[tests$score == primary, , drop = FALSE]
  if (nrow(tp))
    log_message(sprintf(
      "  PREMISE [%s]: %s. AUC P(Viral>NIU) = %.3f [%.3f, %.3f], Hedges g = %.2f, Wilcoxon p = %.3g",
      primary, toupper(tp$premise[1]), tp$auc_alt_gt_ref[1], tp$auc_ci_lo[1],
      tp$auc_ci_hi[1], tp$hedges_g[1], tp$wilcox_p[1]))
  if (!is.null(collin)) {
    cp <- collin[collin$score == primary, , drop = FALSE]
    n_ab <- sum(cp$flag == "abort", na.rm = TRUE)
    log_message(sprintf(
      "  IDENTIFIABILITY [%s]: subject VIF = %.2f; %d / %d pseudobulk strata above the abort threshold.",
      primary,
      suppressWarnings(cp$vif[cp$level == "subject"][1]),
      n_ab, nrow(cp)))
    if (n_ab > 0L)
      log_message("  WARN: strata flagged 'abort' will NOT be stress-adjusted ",
                  "by R/46b. Report them as non-identifiable rather than ",
                  "quoting an unstable coefficient.")
  }

  invisible(list(per_cell = per_cell, subject = subject_df,
                 tests = tests, collinearity = collin, per_pseudobulk = pb))
}
