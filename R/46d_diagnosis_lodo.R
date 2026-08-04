# R/46d_diagnosis_lodo.R
#
# Leave-one-diagnosis-out sensitivity on the NIU arm. 
#
# What this does: drop one NIU sub-diagnosis at a time, re-run the eye
# pseudobulk DGE and the myeloid / T cell per-substate PC1, and report how far
# each result moves. The claim to defend is that no single sub-diagnosis is
# load-bearing, and in particular that the result does not depend on the two
# BSCR subjects.
#
# Subject counts in the eye arm, which set what is testable:
#   Idiopathic 11, HLA_B27 4, BSCR 2, VKH 2, JIA 1, SLE 1   (viral: VZV_ARN 9,
#   CMV_CRN 1, HSV1 1, HSV2 1, HTLV1 1)
# Dropping Idiopathic removes over half the NIU arm, so that fold is reported as
# a bound rather than a sensitivity test. Dropping JIA or SLE removes a single
# subject and is nearly a leave-one-subject-out fold, already covered by R/46c.
# The informative folds are BSCR, VKH and HLA_B27.
#
# Nothing here overwrites a published table. Everything routes through
# .sens_write() and lands under outputs/tables/stress_sensitivity/lodo/.
#
# Entry point: run_diagnosis_lodo(cfg), gated on cfg$steps$diagnosis_lodo.
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

.lodo_dir <- function(cfg) {
  file.path((cfg$paths_stress_sensitivity %||%
               list(tables = "outputs/tables/stress_sensitivity"))$tables, "lodo")
}

# Which NIU sub-diagnoses are worth a fold, and how to read each one.
.lodo_folds <- function(cfg, meta) {
  niu <- as.character(cfg$etiology_groups$niu %||%
                        c("Idiopathic", "HLA_B27", "VKH", "BSCR", "JIA", "SLE"))
  m <- meta[meta$Phenotype_2 %in% "NIU" & meta$Etiology %in% niu, , drop = FALSE]
  if (!nrow(m)) return(NULL)
  tab <- tapply(m$Subject, m$Etiology, function(x) length(unique(x)))
  n_total <- length(unique(m$Subject))
  data.frame(
    etiology = names(tab),
    n_subjects = as.integer(tab),
    frac_of_niu = as.numeric(tab) / n_total,
    # A fold that removes most of the arm is not a sensitivity test, and a fold
    # that removes one subject is a jackknife fold we already have. Label both
    # so the summary is not read as though every row carries equal weight.
    fold_class = ifelse(as.numeric(tab) / n_total > 0.40, "majority_bound",
                 ifelse(as.integer(tab) < 2L, "single_subject", "informative")),
    stringsAsFactors = FALSE)[order(-as.integer(tab)), ]
}

# --- PC1 folds --------------------------------------------------------------

.lodo_pca <- function(cfg, target, folds, out_dir, floor = NULL) {
  p <- get_target_paths(cfg, target)
  op <- file.path(p$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(op)) {
    log_message("  lodo[", target, "]: object missing; skipping.")
    return(NULL)
  }
  log_message("  lodo[", target, "]: loading object")
  obj <- readRDS(op)
  meta <- obj[[]]
  cpcfg <- cfg$compartment_pca %||% list()
  floor <- as.integer(floor %||% .pca_min_cells(cfg, target))
  args <- list(min_cells_per_pb = floor,
               min_gene_count = as.integer(cpcfg$min_gene_count %||% 10L),
               hvg_n          = as.integer(cpcfg$hvg_n %||% 2000L),
               n_pcs          = as.integer(cpcfg$n_pcs %||% 5L),
               vst_blind      = isTRUE(cpcfg$vst_blind),
               pc1_split_fdr  = as.numeric(cpcfg$pc1_split_fdr %||% 0.05))

  full <- do.call(compute_per_substate_pca, c(list(obj = obj), args))
  if (is.null(full)) { rm(obj); gc(verbose = FALSE); return(NULL) }

  rows <- list()
  for (i in seq_len(nrow(folds))) {
    eti <- folds$etiology[i]
    keep <- rownames(meta)[!(meta$Etiology %in% eti)]
    if (length(keep) < 200L) next
    log_message("    drop ", eti, " (", folds$n_subjects[i], " subjects)")
    sub <- subset(obj, cells = keep)
    res <- tryCatch(do.call(compute_per_substate_pca, c(list(obj = sub), args)),
                    error = function(e) NULL)
    rm(sub); gc(verbose = FALSE)
    if (is.null(res)) next
    for (ss in unique(res$scores$substate)) {
      a <- res$scores[res$scores$substate == ss, ]
      b <- full$scores[full$scores$substate == ss, ]
      j <- dplyr::inner_join(
        dplyr::select(a, sample = "sample", fold = "PC1_oriented",
                      grp = "Phenotype_2"),
        dplyr::select(b, sample = "sample", full = "PC1_oriented"), by = "sample")
      if (nrow(j) < 4L) next
      sa <- res$significance[res$significance$substate == ss, ]
      sf <- full$significance[full$significance$substate == ss, ]
      rows[[length(rows) + 1]] <- data.frame(
        dropped_etiology = eti, fold_class = folds$fold_class[i],
        n_subjects_dropped = folds$n_subjects[i],
        target = target, min_cells_per_pb = floor,
        substate = as.character(ss),
        n_NIU = if (nrow(sa)) sa$n_NIU[1] else NA_integer_,
        n_Viral = if (nrow(sa)) sa$n_Viral[1] else NA_integer_,
        q_value = if (nrow(sa)) sa$q_value[1] else NA_real_,
        separating = if (nrow(sa)) isTRUE(sa$separating[1]) else NA,
        separating_full = if (nrow(sf)) isTRUE(sf$separating[1]) else NA,
        cohens_d = .loo_cohens_d(j$fold, j$grp),
        cohens_d_full = .loo_cohens_d(j$full, j$grp),
        pc1_cor_to_full = suppressWarnings(stats::cor(j$fold, j$full)),
        stringsAsFactors = FALSE)
    }
  }
  rm(obj); gc(verbose = FALSE)
  if (!length(rows)) return(NULL)
  out <- dplyr::bind_rows(rows)
  out$pc1_sign_flipped <- out$pc1_cor_to_full < 0
  out$cohens_d_ratio <- abs(out$cohens_d) / abs(out$cohens_d_full)
  out
}

# --- Eye DGE folds ----------------------------------------------------------

.lodo_eye_dge <- function(cfg, folds, out_dir) {
  pub <- .snap_read(cfg, "eye", "DGE_pseudobulk_Autoimmune_vs_Viral.csv")
  if (is.null(pub)) {
    log_message("  lodo[eye_dge]: published snapshot missing; skipping.")
    return(NULL)
  }
  paths_eye <- get_target_paths(cfg, "eye")
  op <- file.path(paths_eye$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(op)) return(NULL)
  log_message("  lodo[eye_dge]: loading eye object")
  obj <- readRDS(op)
  meta <- obj[[]]
  cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                   "knn.leiden.cluster" else "seurat_clusters"
  padj_thr <- as.numeric(cfg$dge$padj_threshold %||% 0.05)
  lfc_thr  <- as.numeric(cfg$dge$lfc_threshold  %||% 0.5)

  rows <- list()
  for (i in seq_len(nrow(folds))) {
    eti <- folds$etiology[i]
    keep <- rownames(meta)[!(meta$Etiology %in% eti)]
    if (length(keep) < 200L) next
    log_message("    drop ", eti)
    sub <- subset(obj, cells = keep)
    d <- tryCatch(run_pseudobulk_deseq2(sub, group_col = "Phenotype_2",
                                        group1 = "NIU", group2 = "Viral",
                                        cluster_col = cluster_col, cfg = cfg,
                                        target = "eye", force_simple = TRUE),
                  error = function(e) NULL)
    rm(sub); gc(verbose = FALSE)
    if (is.null(d) || !nrow(d)) next
    d$dropped_etiology <- eti
    .sens_write(d, out_dir,
                paste0("DGE_pseudobulk_Autoimmune_vs_Viral_drop_", eti))
    a <- pub[pub$cluster == "global", ]; b <- d[d$cluster == "global", ]
    if (!nrow(a) || !nrow(b)) next
    m <- dplyr::inner_join(
      dplyr::select(a, gene = "gene", lfc_pub = "log2FoldChange", padj_pub = "padj"),
      dplyr::select(b, gene = "gene", lfc_fold = "log2FoldChange", padj_fold = "padj"),
      by = "gene")
    if (nrow(m) < 50L) next
    sp <- !is.na(m$padj_pub) & m$padj_pub < padj_thr & abs(m$lfc_pub) > lfc_thr
    sf <- !is.na(m$padj_fold) & m$padj_fold < padj_thr & abs(m$lfc_fold) > lfc_thr
    top <- m[sp, ][order(m$padj_pub[sp]), ]; top <- utils::head(top, 50L)
    rows[[length(rows) + 1]] <- data.frame(
      dropped_etiology = eti, fold_class = folds$fold_class[i],
      n_subjects_dropped = folds$n_subjects[i],
      n_common = nrow(m), n_sig_published = sum(sp), n_sig_fold = sum(sf),
      lfc_spearman = suppressWarnings(
        stats::cor(m$lfc_pub, m$lfc_fold, method = "spearman")),
      deg_jaccard = if (sum(sp | sf)) sum(sp & sf) / sum(sp | sf) else NA_real_,
      sign_flip_rate = if (sum(sp))
        mean(sign(m$lfc_pub[sp]) != sign(m$lfc_fold[sp])) else NA_real_,
      top50_retention = if (nrow(top))
        mean(!is.na(top$padj_fold) & top$padj_fold < padj_thr) else NA_real_,
      stringsAsFactors = FALSE)
  }
  rm(obj); gc(verbose = FALSE)
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}


viz_diagnosis_lodo <- function(cfg, dge = NULL, pca = NULL) {
  out_dir <- .lodo_dir(cfg)
  vd <- (cfg$paths_stress_sensitivity %||%
           list(viz = "outputs/viz/stress_sensitivity"))$viz
  ensure_dir(vd)
  rd <- function(f) { p <- file.path(out_dir, f)
                      if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) }
  dge <- dge %||% rd("lodo_eye_dge_concordance.csv")
  pca <- pca %||% rd("lodo_pca_pc1.csv")
  cls_pal <- c(informative = "#1B7837", majority_bound = "#E08214",
               single_subject = "grey60")

  # --- A. eye DGE agreement ------------------------------------------------
  if (!is.null(dge) && nrow(dge)) {
    d <- dge
    d$lab <- paste0(d$dropped_etiology, "\n(-", d$n_subjects_dropped, " pts)")
    long <- rbind(
      data.frame(lab = d$lab, fold_class = d$fold_class,
                 metric = "log2FC agreement (Spearman rho)",
                 value = d$lfc_spearman, stringsAsFactors = FALSE),
      data.frame(lab = d$lab, fold_class = d$fold_class,
                 metric = "DEG set overlap (Jaccard)",
                 value = d$deg_jaccard, stringsAsFactors = FALSE),
      data.frame(lab = d$lab, fold_class = d$fold_class,
                 metric = "Top 50 genes still significant",
                 value = d$top50_retention, stringsAsFactors = FALSE))
    long$lab <- factor(long$lab, levels = d$lab[order(-d$n_subjects_dropped)])
    pa <- ggplot2::ggplot(long, ggplot2::aes(.data$lab, .data$value)) +
      ggplot2::geom_col(ggplot2::aes(fill = .data$fold_class), width = 0.65) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$value)),
                         vjust = -0.4, size = 2.6) +
      ggplot2::facet_wrap(~ .data$metric) +
      ggplot2::scale_fill_manual(values = cls_pal, name = "fold type") +
      ggplot2::coord_cartesian(ylim = c(0, 1.12)) +
      ggplot2::labs(
        title = "Eye pseudobulk DGE, leave-one-diagnosis-out",
        subtitle = paste("Agreement with the full NIU versus viral analysis",
                         "after removing each NIU sub-diagnosis.",
                         "\nOnly the green folds are true sensitivity tests.",
                         "Dropping idiopathic removes over half the NIU arm and",
                         "is a power bound."),
        x = NULL, y = NULL) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = 7.5),
                     axis.text.x = ggplot2::element_text(size = 7),
                     legend.position = "bottom")
    save_pdf_png(pa, file.path(vd, "lodo_eye_dge_concordance"), w = 10, h = 4.8)
  }

  # --- B. per-substate PC1 separation grid ---------------------------------
  if (!is.null(pca) && nrow(pca)) {
    p <- pca[pca$separating_full %in% TRUE, , drop = FALSE]
    if (nrow(p)) {
      p$panel <- if ("min_cells_per_pb" %in% names(p))
        paste0(p$target, " (floor ", p$min_cells_per_pb, ")") else p$target
      p$state <- ifelse(is.na(p$separating), "not testable",
                 ifelse(p$separating %in% TRUE, "still separates",
                        "loses separation"))
      p$lab <- paste0(p$dropped_etiology, "\n(-", p$n_subjects_dropped, ")")
      p$lab <- factor(p$lab, levels = unique(
        p$lab[order(-p$n_subjects_dropped)]))
      pb <- ggplot2::ggplot(p, ggplot2::aes(.data$lab,
                                            factor(.data$substate))) +
        ggplot2::geom_tile(ggplot2::aes(fill = .data$state), colour = "white",
                           linewidth = 0.6) +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(.data$q_value), "",
                             ifelse(.data$q_value < 0.001, "<.001",
                                    sprintf("%.3f", .data$q_value)))),
                           size = 2.5) +
        ggplot2::facet_wrap(~ .data$panel, scales = "free_y") +
        ggplot2::scale_fill_manual(
          values = c(`still separates` = "#CDE8CD",
                     `loses separation` = "#E8453B",
                     `not testable` = "grey85"), name = NULL) +
        ggplot2::labs(
          title = "Per-substate PC1 separation after dropping each NIU sub-diagnosis",
          subtitle = paste("Restricted to substates that separate in the full",
                           "data. Cell text is the BH q value."),
          x = "dropped NIU sub-diagnosis (subjects removed)", y = "substate") +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                       plot.subtitle = ggplot2::element_text(size = 7.5),
                       axis.text.x = ggplot2::element_text(size = 7),
                       panel.grid = ggplot2::element_blank(),
                       legend.position = "bottom")
      save_pdf_png(pb, file.path(vd, "lodo_pc1_separation_grid"), w = 10, h = 5.5)
    }
  }
  invisible(NULL)
}

run_diagnosis_lodo <- function(cfg, targets = NULL) {
  targets <- as.character(targets %||% cfg$diagnosis_lodo$targets %||%
                            c("eye_dge", "myeloid", "tcell"))
  out_dir <- .lodo_dir(cfg)
  ensure_dir(out_dir)

  paths_eye <- get_target_paths(cfg, "eye")
  f_pc <- file.path(paths_eye$results_tables, "stress_ucell_per_cell.csv")
  meta <- if (file.exists(f_pc))
    utils::read.csv(f_pc, stringsAsFactors = FALSE)[, c("Subject", "Phenotype_2",
                                                        "Etiology")] else NULL
  if (is.null(meta)) {
    log_message("  lodo: need stress_ucell_per_cell.csv for the subject table.")
    return(invisible(FALSE))
  }
  folds <- .lodo_folds(cfg, meta)
  if (is.null(folds) || !nrow(folds)) return(invisible(FALSE))
  .sens_write(folds, out_dir, "lodo_fold_definition")
  for (i in seq_len(nrow(folds)))
    log_message(sprintf("  fold %-11s n = %2d subjects (%.0f%% of NIU) [%s]",
                        folds$etiology[i], folds$n_subjects[i],
                        100 * folds$frac_of_niu[i], folds$fold_class[i]))

  if ("eye_dge" %in% targets) {
    dge <- .lodo_eye_dge(cfg, folds, out_dir)
    if (!is.null(dge)) {
      .sens_write(dge, out_dir, "lodo_eye_dge_concordance")
      for (i in seq_len(nrow(dge)))
        log_message(sprintf(
          "  EYE DGE drop %-11s [%s]: lfc rho = %.3f, DEG Jaccard = %.3f, sign flips = %.1f%%, top50 kept = %.0f%%",
          dge$dropped_etiology[i], dge$fold_class[i], dge$lfc_spearman[i],
          dge$deg_jaccard[i], 100 * dge$sign_flip_rate[i],
          100 * dge$top50_retention[i]))
    }
  }

  jobs <- list()
  for (tg in intersect(targets, c("myeloid", "tcell"))) {
    fl <- .pca_min_cells(cfg, tg)
    jobs[[length(jobs) + 1]] <- list(target = tg, floor = fl)
    if (identical(tg, "myeloid")) {
      f3e <- as.integer(
        (cfg$compartment_pca$min_cells_per_subject_substate_myeloid_f3e) %||% fl)
      if (!identical(f3e, fl))
        jobs[[length(jobs) + 1]] <- list(target = tg, floor = f3e)
    }
  }
  pca <- dplyr::bind_rows(Filter(Negate(is.null), lapply(
    jobs, function(j) .lodo_pca(cfg, j$target, folds, out_dir, floor = j$floor))))
  if (nrow(pca)) {
    .sens_write(pca, out_dir, "lodo_pca_pc1")
    summ <- pca |>
      dplyr::filter(.data$separating_full %in% TRUE) |>
      dplyr::group_by(.data$target, .data$min_cells_per_pb, .data$substate) |>
      dplyr::summarise(
        n_folds = dplyr::n(),
        frac_folds_separating = mean(.data$separating %in% TRUE),
        min_abs_pc1_cor = min(abs(.data$pc1_cor_to_full), na.rm = TRUE),
        n_sign_flips = sum(.data$pc1_sign_flipped %in% TRUE),
        min_cohens_d_ratio = min(.data$cohens_d_ratio, na.rm = TRUE),
        worst_fold = .data$dropped_etiology[which.min(.data$cohens_d_ratio)],
        .groups = "drop")
    .sens_write(summ, out_dir, "lodo_pca_summary")

    # The specific claim the reviewer asked us to defend.
    bscr <- pca[pca$dropped_etiology == "BSCR" & pca$separating_full %in% TRUE, ]
    if (nrow(bscr))
      log_message(sprintf(
        "  BSCR CHECK: dropping both BSCR subjects leaves %d of %d separating substates separating; min |PC1 r| to full = %.3f; min effect-size ratio = %.2f",
        sum(bscr$separating %in% TRUE), nrow(bscr),
        min(abs(bscr$pc1_cor_to_full), na.rm = TRUE),
        min(bscr$cohens_d_ratio, na.rm = TRUE)))
    lost <- summ[summ$frac_folds_separating < 1, , drop = FALSE]
    if (nrow(lost))
      for (i in seq_len(nrow(lost)))
        log_message(sprintf(
          "  SENSITIVE: %s substate %s separates in only %.0f%% of folds, worst = drop %s",
          lost$target[i], lost$substate[i],
          100 * lost$frac_folds_separating[i], lost$worst_fold[i]))
    else log_message("  No separating substate depends on any single NIU sub-diagnosis.")
  }
  tryCatch(viz_diagnosis_lodo(cfg),
           error = function(e)
             log_message("  lodo figures failed: ", conditionMessage(e)))
  invisible(list(folds = folds, pca = pca))
}
