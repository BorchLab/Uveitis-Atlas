# R/79_tcr_sharpened_bridge.R
# Figure 5 Panel K: sharpened myeloid <-> T cell PC1 bridge restricted to
# disease-associated TCR clone groups.
#
# Motivation: Figure 4 panel E established that subject myeloid PC1 residual
# correlates with subject T cell PC1 residual (r = 0.70, p = 0.001, n = 20).
# That correlation uses ALL T cells per subject. The Figure 5 question is:
# does the bridge sharpen when we restrict the T cell pseudobulk to cells
# carrying disease-associated TCR motifs (GLIPH-NIU-enriched, expanded,
# HLA-B27 pathogenic)? If yes, the motif-defined clones are the ones being
# instructed by the myeloid axis the LIANA analysis surfaced. If no, the
# bridge is a substate-level phenomenon and the TCR is downstream.
#
# Implementation:
#   1. Read the eye T cell compartment object + per-cell clone-group flags
#      (computed identically to R/70_tcr_genex_signatures.R::.tgs_label_clones).
#   2. For each clone group, build per-(subject, substate) pseudobulks
#      restricted to cells in the group. Drop subject_substate cells <
#      cfg$tcr_advanced$fig5$panel_k_min_cells_per_subject. Run the same
#      DESeq2 VST + PCA pipeline that R/45_compartment_pca uses so the PC1
#      orientation matches the all-cells baseline.
#   3. Read the myeloid all-cells PC1 from outputs/tables/eye/myeloid/
#      pca_subject_scores.csv (subject-level aggregate per .aggregate_subject_pc1
#      in R/46) and join with each clone-group T cell PC1.
#   4. Run the same partial-correlation + bootstrap battery as R/46 so the
#      panel reports comparable r / p / CI to Figure 4 panel E.
#   5. Write per-(clone_group, weighting) subject scatter and correlation
#      table; render Panel K as a multi-facet scatter.
#
# Outputs:
#   outputs/tables/cross_compartment/pc1_bridge_tcr_sharpened_scores.csv
#   outputs/tables/cross_compartment/pc1_bridge_tcr_sharpened_correlation.csv
#   outputs/viz/eye/tcell/tcr_advanced/sharpened_bridge_scatter.{pdf,png}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Compute the same clone-group flags R/70_tcr_genex_signatures.R uses, plus
# HLA-B27 pathogenic. Returns the cell-level metadata frame with appended
# columns: grp_expanded, grp_public, grp_gliph_niu, grp_gliph_viral, grp_hla_b27.
.tsb_label_clones <- function(obj, gliph_niu_cdr3s, gliph_viral_cdr3s) {
  meta <- obj@meta.data
  cs <- as.character(meta$cloneSize)
  meta$grp_expanded <- grepl("^Large|^Hyperexpanded", cs)
  meta$grp_expanded[is.na(meta$grp_expanded)] <- FALSE

  ct <- as.character(meta$CTstrict)
  sb <- as.character(meta$Subject)
  ct_subj <- tapply(sb, ct, function(x) length(unique(x[!is.na(x)])))
  meta$grp_public <- !is.na(ct) & ct_subj[ct] >= 2L
  meta$grp_public[is.na(meta$grp_public)] <- FALSE

  trb <- stringr::str_split(as.character(meta$CTaa), "_", simplify = TRUE)
  trb_aa <- if (ncol(trb) >= 2) trb[, 2] else rep(NA_character_, nrow(meta))
  meta$grp_gliph_niu   <- !is.na(trb_aa) & trb_aa %in% gliph_niu_cdr3s
  meta$grp_gliph_viral <- !is.na(trb_aa) & trb_aa %in% gliph_viral_cdr3s

  meta$grp_hla_b27 <- flag_hla_b27_pathogenic_tcr(meta$CTgene, meta$CTaa)
  meta$grp_hla_b27[is.na(meta$grp_hla_b27)] <- FALSE

  meta
}

# Reuse R/46's helpers. Re-implementing the per-(subject, substate) PC1
# pipeline from R/45_compartment_pca would re-derive the same VST+PCA, but
# the cleaner path is to leverage the cell-level filter. We pseudobulk only
# the cells in the clone group, fit DESeq2 VST + PCA per substate, orient
# PC1 by Phenotype_2, then aggregate to subject-level the same way R/46 does.
.tsb_pseudobulk_pca_per_substate <- function(obj, cells_in_group, cfg,
                                              min_cells_per_subject = 30L) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    log_message("  DESeq2 not installed; sharpened bridge skipped.")
    return(NULL)
  }
  meta <- obj@meta.data
  cells_in_group <- intersect(cells_in_group, colnames(obj))
  if (length(cells_in_group) == 0L) return(NULL)
  obj_sub <- obj[, cells_in_group]
  meta_sub <- obj_sub@meta.data

  # substate_key fallback identical to R/70_tcr_genex_signatures.R
  if (!"substate_key" %in% colnames(meta_sub)) {
    meta_sub$substate_key <- get_substate_key_vector(obj_sub, "tcell")
  }
  if (!"Subject_Timepoint" %in% colnames(meta_sub)) {
    meta_sub$Subject_Timepoint <- paste0(meta_sub$Subject, "_",
                                         meta_sub$orig.ident)
  }

  # Count cells per (Subject_Timepoint, substate); drop below floor.
  sample_col <- "Subject_Timepoint"
  group_tbl <- meta_sub |>
    dplyr::group_by(.data$substate_key, .data[[sample_col]],
                    .data$Subject, .data$Phenotype_2) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::filter(.data$n_cells >= min_cells_per_subject)
  if (nrow(group_tbl) == 0L) {
    log_message("  No (subject, substate) bins met n >= ",
                min_cells_per_subject, " floor in clone group.")
    return(NULL)
  }

  # Per-substate pseudobulk + PCA. PC1 oriented so NIU subjects score on the
  # positive pole (matches R/45 convention).
  rows <- list()
  counts <- as.matrix(SeuratObject::GetAssayData(obj_sub, layer = "counts"))
  for (sub_id in unique(group_tbl$substate_key)) {
    keep <- group_tbl[group_tbl$substate_key == sub_id, ]
    if (nrow(keep) < 4L) next
    sample_keep <- keep[[sample_col]]
    cells_use <- meta_sub[[sample_col]] %in% sample_keep &
                 meta_sub$substate_key == sub_id
    if (sum(cells_use) < min_cells_per_subject) next
    cnt <- counts[, cells_use, drop = FALSE]
    s_id <- meta_sub[[sample_col]][cells_use]
    pb <- vapply(unique(s_id),
                 function(s) rowSums(cnt[, s_id == s, drop = FALSE]),
                 numeric(nrow(cnt)))
    colnames(pb) <- unique(s_id)
    keep_gene <- rowSums(pb >= cfg$compartment_pca$min_gene_count %||% 10) >= 2L
    pb <- pb[keep_gene, , drop = FALSE]
    if (nrow(pb) < 50L || ncol(pb) < 4L) next

    col_meta <- keep[match(colnames(pb), keep[[sample_col]]), ]
    col_meta$Phenotype_2 <- factor(col_meta$Phenotype_2,
                                   levels = c("Viral", "NIU"))
    dds <- tryCatch(
      DESeq2::DESeqDataSetFromMatrix(countData = pb,
                                     colData = col_meta,
                                     design = ~ Phenotype_2),
      error = function(e) NULL)
    if (is.null(dds)) next
    vsd <- tryCatch(DESeq2::vst(dds, blind = FALSE),
                    error = function(e) NULL)
    if (is.null(vsd)) next
    mat <- SummarizedExperiment::assay(vsd)
    hvg_n <- cfg$compartment_pca$hvg_n %||% 2000L
    if (nrow(mat) > hvg_n) {
      vv <- matrixStats::rowVars(mat)
      mat <- mat[order(-vv)[1:hvg_n], , drop = FALSE]
    }
    pr <- prcomp(t(mat), center = TRUE, scale. = FALSE)
    pc1 <- pr$x[, 1]
    # Orient: positive => NIU pole
    niu_mean   <- mean(pc1[col_meta$Phenotype_2 == "NIU"],   na.rm = TRUE)
    viral_mean <- mean(pc1[col_meta$Phenotype_2 == "Viral"], na.rm = TRUE)
    sign_flip  <- if (niu_mean < viral_mean) -1 else 1
    pc1 <- pc1 * sign_flip

    rows[[sub_id]] <- data.frame(
      substate    = sub_id,
      sample      = colnames(pb),
      subject     = col_meta$Subject,
      Phenotype_2 = as.character(col_meta$Phenotype_2),
      n_cells     = col_meta$n_cells,
      PC1_oriented = unname(pc1),
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0L) return(NULL)
  dplyr::bind_rows(rows)
}

# Subject-level aggregate matches R/46 .aggregate_subject_pc1 (weighted).
.tsb_aggregate_subject <- function(scores) {
  st1 <- scores |>
    dplyr::group_by(.data$subject, .data$substate, .data$Phenotype_2) |>
    dplyr::summarise(
      PC1 = stats::weighted.mean(.data$PC1_oriented, w = .data$n_cells,
                                 na.rm = TRUE),
      n_cells = sum(.data$n_cells, na.rm = TRUE),
      .groups = "drop")
  st1 |>
    dplyr::group_by(.data$subject, .data$Phenotype_2) |>
    dplyr::summarise(
      PC1     = stats::weighted.mean(.data$PC1, w = .data$n_cells,
                                     na.rm = TRUE),
      n_cells = sum(.data$n_cells, na.rm = TRUE),
      n_substates = dplyr::n(),
      .groups = "drop")
}

# Lightweight partial-correlation battery copied from R/46 (avoids importing
# the private helper); keeps numerics comparable to the Figure 4 panel E
# baseline.
.tsb_cor_battery <- function(x, y, n_boot = 500L, seed = 42L) {
  out <- tibble::tibble(
    n = length(x),
    pearson_r = NA_real_, pearson_p = NA_real_,
    pearson_ci_lo = NA_real_, pearson_ci_hi = NA_real_)
  if (length(x) < 3L) return(out)
  pt <- tryCatch(stats::cor.test(x, y, method = "pearson"),
                 error = function(e) NULL)
  if (!is.null(pt)) {
    out$pearson_r <- unname(pt$estimate)
    out$pearson_p <- pt$p.value
  }
  set.seed(seed)
  if (n_boot >= 50L) {
    n <- length(x)
    rs <- vapply(seq_len(n_boot), function(i) {
      idx <- sample.int(n, n, replace = TRUE)
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
  out
}

# Partial r controlling for Phenotype_2 (same idea as R/46's helper).
.tsb_partial_cor <- function(x, y, z) {
  if (length(x) < 4L || length(unique(z)) < 2L)
    return(list(r = NA_real_, p = NA_real_))
  rx <- tryCatch(stats::resid(stats::lm(x ~ factor(z))),
                 error = function(e) NULL)
  ry <- tryCatch(stats::resid(stats::lm(y ~ factor(z))),
                 error = function(e) NULL)
  if (is.null(rx) || is.null(ry))
    return(list(r = NA_real_, p = NA_real_))
  ct <- tryCatch(stats::cor.test(rx, ry, method = "pearson"),
                 error = function(e) NULL)
  if (is.null(ct)) return(list(r = NA_real_, p = NA_real_))
  list(r = unname(ct$estimate), p = ct$p.value)
}

run_tcr_sharpened_bridge <- function(cfg) {
  if (!isTRUE(cfg$steps$tcr_sharpened_bridge)) {
    log_message("tcr_sharpened_bridge disabled. Skipping.")
    return(invisible(TRUE))
  }
  log_message("Starting TCR-sharpened cross-compartment bridge (Panel K)...")

  paths_tcell <- get_target_paths(cfg, "tcell")
  paths_myel  <- get_target_paths(cfg, "myeloid")
  cc_paths    <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/")
  ensure_dir(cc_paths$tables)
  base_viz <- file.path(paths_tcell$viz_dir, "tcr_advanced")
  ensure_dir(base_viz)

  tcell_rds <- file.path(paths_tcell$results_objects,
                         "IntegratedSeuratObject.rds")
  myel_csv  <- file.path(paths_myel$results_tables, "pca_subject_scores.csv")
  if (!file.exists(tcell_rds) || !file.exists(myel_csv)) {
    log_message("  Required inputs missing (tcell object or myeloid PCA ",
                "scores); skipping.")
    return(invisible(FALSE))
  }

  obj <- readRDS(tcell_rds)
  myel <- utils::read.csv(myel_csv, stringsAsFactors = FALSE)

  # Load GLIPH enrichment to define clone groups (NIU- and Viral-enriched).
  gliph_rds <- file.path(cfg$paths$results_objects, "ImmGLIPHResults.rds")
  gliph_niu_cdr3s   <- character(0)
  gliph_viral_cdr3s <- character(0)
  fdr_cut <- cfg$tcr_advanced$gliph_enrichment_fdr %||% 0.1
  if (file.exists(gliph_rds)) {
    g <- readRDS(gliph_rds)
    if (!is.null(g$enrich) && !is.null(g$clusters)) {
      en <- g$enrich
      cls <- g$clusters
      niu_ids   <- en$cluster_id[!is.na(en$FDR) & en$FDR < fdr_cut &
                                  grepl("NIU",   en$direction)]
      viral_ids <- en$cluster_id[!is.na(en$FDR) & en$FDR < fdr_cut &
                                  grepl("Viral", en$direction)]
      gliph_niu_cdr3s   <- unique(cls$CDR3b[cls$cluster_id %in% niu_ids])
      gliph_viral_cdr3s <- unique(cls$CDR3b[cls$cluster_id %in% viral_ids])
    }
  }
  log_message("  GLIPH-enriched CDR3b counts: NIU=",   length(gliph_niu_cdr3s),
              ", Viral=", length(gliph_viral_cdr3s))

  meta <- .tsb_label_clones(obj, gliph_niu_cdr3s, gliph_viral_cdr3s)
  obj@meta.data <- meta

  # Subject-level all-cells (baseline) myeloid PC1, weighted aggregate.
  m_subj_all <- .tsb_aggregate_subject(
    myel |> dplyr::rename(subject = "subject"))
  if (!nrow(m_subj_all)) {
    log_message("  Myeloid PCA score CSV has no usable rows; skipping.")
    return(invisible(FALSE))
  }

  # All-cells T baseline (matches R/46 output for sanity-checking).
  tcel_csv <- file.path(paths_tcell$results_tables, "pca_subject_scores.csv")
  t_subj_all <- if (file.exists(tcel_csv)) {
    tcel_all <- utils::read.csv(tcel_csv, stringsAsFactors = FALSE)
    .tsb_aggregate_subject(tcel_all)
  } else NULL

  clone_groups <- cfg$tcr_advanced$fig5$panel_k_clone_groups %||%
                   c("gliph_niu", "gliph_viral", "expanded", "hla_b27")
  min_cells <- cfg$tcr_advanced$fig5$panel_k_min_cells_per_subject %||%
                 cfg$compartment_pca$min_cells_per_subject_substate_tcell %||% 30L

  scatter_rows <- list()
  cor_rows     <- list()

  add_scatter <- function(label, x, y, z, subj, eti) {
    if (length(x) < 3L) return(NULL)
    cb   <- .tsb_cor_battery(x, y, n_boot = 500L)
    part <- .tsb_partial_cor(x, y, z)
    cor_rows[[label]] <<- tibble::tibble(
      clone_group = label,
      n = length(x),
      pearson_r = cb$pearson_r,
      pearson_p = cb$pearson_p,
      pearson_ci_lo = cb$pearson_ci_lo,
      pearson_ci_hi = cb$pearson_ci_hi,
      partial_r = part$r,
      partial_p = part$p)
    scatter_rows[[label]] <<- tibble::tibble(
      clone_group = label,
      subject = subj, Phenotype_2 = z,
      myeloid_pc1 = x, tcell_pc1 = y)
  }

  if (!is.null(t_subj_all)) {
    jn0 <- dplyr::inner_join(
      m_subj_all |> dplyr::select("subject", "Phenotype_2",
                                  myeloid_pc1 = "PC1"),
      t_subj_all |> dplyr::select("subject", "Phenotype_2",
                                   tcell_pc1 = "PC1"),
      by = c("subject", "Phenotype_2"))
    add_scatter("all_cells_baseline",
                jn0$myeloid_pc1, jn0$tcell_pc1,
                jn0$Phenotype_2, jn0$subject, jn0$Phenotype_2)
  }

  for (cg in clone_groups) {
    flag_col <- paste0("grp_", cg)
    if (!flag_col %in% colnames(meta)) {
      log_message("  Skipping clone_group '", cg, "': flag column missing.")
      next
    }
    cells <- rownames(meta)[meta[[flag_col]]]
    log_message("  Clone group '", cg, "': ", length(cells), " cells.")
    if (length(cells) < min_cells * 3L) {
      log_message("    Too few cells (< 3 x min_cells_per_subject); skipping.")
      next
    }
    scores <- .tsb_pseudobulk_pca_per_substate(
      obj, cells_in_group = cells, cfg = cfg,
      min_cells_per_subject = min_cells)
    if (is.null(scores) || nrow(scores) == 0L) {
      log_message("    Pseudobulk PCA returned empty for '", cg, "'.")
      next
    }
    t_subj_cg <- .tsb_aggregate_subject(scores)
    if (!nrow(t_subj_cg)) next
    jn <- dplyr::inner_join(
      m_subj_all |> dplyr::select("subject", "Phenotype_2",
                                  myeloid_pc1 = "PC1"),
      t_subj_cg  |> dplyr::select("subject", "Phenotype_2",
                                   tcell_pc1 = "PC1"),
      by = c("subject", "Phenotype_2"))
    log_message("    Paired subjects: ", nrow(jn))
    add_scatter(cg, jn$myeloid_pc1, jn$tcell_pc1,
                jn$Phenotype_2, jn$subject, jn$Phenotype_2)
  }

  scatter <- if (length(scatter_rows)) dplyr::bind_rows(scatter_rows) else
               tibble::tibble()
  cor_tbl <- if (length(cor_rows))     dplyr::bind_rows(cor_rows)     else
               tibble::tibble()

  utils::write.csv(scatter,
                   file.path(cc_paths$tables,
                             "pc1_bridge_tcr_sharpened_scores.csv"),
                   row.names = FALSE)
  utils::write.csv(cor_tbl,
                   file.path(cc_paths$tables,
                             "pc1_bridge_tcr_sharpened_correlation.csv"),
                   row.names = FALSE)
  log_message("  Saved sharpened-bridge tables (",
              nrow(scatter), " scatter rows, ",
              nrow(cor_tbl), " correlation rows).")

  # Render Panel K
  if (nrow(scatter) > 0L && requireNamespace("ggplot2", quietly = TRUE)) {
    cor_lbl <- cor_tbl |>
      dplyr::mutate(label = sprintf(
        "n=%d  r=%.2f  p=%.2g\\npartial r=%.2f  p=%.2g",
        .data$n, .data$pearson_r, .data$pearson_p,
        .data$partial_r, .data$partial_p))
    scatter$clone_group <- factor(scatter$clone_group,
      levels = c("all_cells_baseline", clone_groups))
    cor_lbl$clone_group <- factor(cor_lbl$clone_group,
      levels = c("all_cells_baseline", clone_groups))
    pal <- if (exists("ETIOLOGY_GROUP_COLORS", inherits = TRUE))
             ETIOLOGY_GROUP_COLORS else
             c(NIU = "#E21F26", Viral = "#397FB9", Healthy = "#A8DADC")
    p <- ggplot2::ggplot(scatter,
                         ggplot2::aes(.data$myeloid_pc1, .data$tcell_pc1)) +
      ggplot2::geom_point(ggplot2::aes(color = .data$Phenotype_2),
                          size = 2.5, alpha = 0.85) +
      ggplot2::geom_smooth(method = "lm", se = TRUE,
                           color = "black", linewidth = 0.7) +
      ggplot2::geom_text(data = cor_lbl,
        ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
        hjust = -0.05, vjust = 1.1, size = 3, inherit.aes = FALSE) +
      ggplot2::facet_wrap(~ .data$clone_group, scales = "free", ncol = 3) +
      ggplot2::scale_color_manual(values = pal, na.value = "grey60",
                                  name = "Phenotype") +
      ggplot2::labs(
        title    = paste0("Figure 5K: myeloid <-> T cell PC1 bridge ",
                          "restricted to TCR-defined clone groups"),
        subtitle = paste0("Per-subject pseudobulk PCA within the eye T cell ",
                          "compartment restricted to cells in each clone ",
                          "group. Subject floor = ", min_cells, " cells in ",
                          "the filtered set per (subject, substate)."),
        x = "Myeloid PC1 (subject-level aggregate, all cells)",
        y = "T cell PC1 (clone-group restricted)") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = 9),
                     strip.text    = ggplot2::element_text(face = "bold"))
    save_pdf_png(p, file.path(base_viz, "sharpened_bridge_scatter"),
                 w = 11, h = 7.5)
  }
  log_message("TCR-sharpened cross-compartment bridge complete.")
  invisible(TRUE)
}
