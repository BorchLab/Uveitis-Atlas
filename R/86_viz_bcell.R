# R/86_viz_bcell.R
# B / plasma compartment figure block. Uses shared helpers in
# R/81_viz_compartment_helpers.R for cluster panels, GSEA heatmap, pathway
# bars, etiology heatmaps, public-clone count, CDR3 length, IGHV usage,
# per-substate volcano grid, and the NIU sub-contrast heatmap.
#
# B cell-specific panels:
#   - Eye-blood BCR overlap boxplot
#   - Isotype usage stacked bar by substate (NIU vs Viral)
#   - Heavy-chain SHM frequency NIU vs Viral within each substate (violin)
#   - RNA UMAP: SHM hex split by Phenotype_2 (NIU vs Viral) +
#     matching NIU-vs-Viral SHM boxplot on the same bcell cell set
#
# All filenames use snake_case descriptions. The full MiloR viz block from
# 82_viz_dispatch.R is wired in at the end.

suppressPackageStartupMessages({
  library(Seurat)
  library(scplotter)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Filter Healthy etiology consistently across Fig 6 panels. Pass cfg so we
# pull cfg$etiology_groups$healthy (defaults to "Healthy"). Restricts both
# the etiology column and the phenotype column where present.
.drop_healthy <- function(df, cfg, etiology_col = "Etiology",
                          phenotype_col = "Phenotype_2") {
  healthy_set <- as.character(cfg$etiology_groups$healthy %||% "Healthy")
  if (!is.null(df) && etiology_col %in% colnames(df))
    df <- df[!(df[[etiology_col]] %in% healthy_set), , drop = FALSE]
  if (!is.null(df) && phenotype_col %in% colnames(df))
    df <- df[df[[phenotype_col]] %in% c("NIU","Viral"), , drop = FALSE]
  df
}

.shm_total_col <- function(meta) {
  # Sum heavy-chain CDR + FWR replacement and silent SHM frequency.
  cols <- intersect(c("mu_freq_cdr_r_heavy", "mu_freq_cdr_s_heavy",
                      "mu_freq_fwr_r_heavy", "mu_freq_fwr_s_heavy"),
                    colnames(meta))
  if (length(cols) == 0) return(NULL)
  rowSums(as.matrix(meta[, cols, drop = FALSE]), na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Isotype stacked bar by substate, NIU vs Viral panels
# ---------------------------------------------------------------------------
.bcell_isotype_by_substate <- function(obj, paths, cfg) {
  meta <- obj@meta.data
  isotype_col <- intersect(c("c_call", "isotype", "Isotype"),
                           colnames(meta))[1]
  if (is.na(isotype_col) || is.null(isotype_col)) return(invisible())
  meta$substate <- substate_labels(cfg, "bcell", meta$knn.leiden.cluster)
  meta <- .drop_healthy(meta, cfg)
  iso <- meta |>
    dplyr::filter(!is.na(.data[[isotype_col]])) |>
    dplyr::count(substate, Phenotype_2, isotype = .data[[isotype_col]])
  if (nrow(iso) == 0) return(invisible())
  p <- ggplot(iso, aes(x = substate, y = n, fill = isotype)) +
       geom_col(position = "fill") +
       facet_wrap(~ Phenotype_2) +
       viridis::scale_fill_viridis(option = "viridis", discrete = TRUE,
                                   name = "Isotype") +
       theme_classic() +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
       labs(title = "B/plasma isotype distribution per substate, NIU vs Viral",
            x = "Substate", y = "Proportion of cells")
  save_pdf_png(p, file.path(viz_subdir(paths, "repertoire"),
                            "bcell_isotype_by_substate_niu_vs_viral"),
               w = 12, h = 6)
}

# Boxplot quantification of heavy-chain SHM, NIU vs Viral, at the SUBJECT
# level. Each point is one subject's mean SHM across their mutated B/plasma
# cells (SHM > 0, i.e. excluding naive IGHM/D cells whose unmutated state
# would otherwise pin the cell-level Q1 at zero and mask the per-subject
# shift). Unpaired Wilcoxon between groups since subjects are NIU XOR Viral.
.bcell_shm_boxplot_niu_vs_viral <- function(obj, paths, cfg,
                                            min_cells = 5L) {
  shm_total <- .shm_total_col(obj@meta.data)
  if (is.null(shm_total)) {
    log_message("  SHM boxplot: no mu_freq columns; skipping.")
    return(invisible())
  }
  meta <- obj@meta.data
  meta$SHM_total <- shm_total
  subj_col <- intersect(c("Subject", "subject", "orig.ident"),
                        colnames(meta))[1]
  if (is.na(subj_col)) {
    log_message("  SHM boxplot: no Subject column; skipping.")
    return(invisible())
  }
  cells <- .drop_healthy(meta, cfg) |>
    dplyr::filter(!is.na(SHM_total), SHM_total > 0,
                  Phenotype_2 %in% c("NIU", "Viral"))
  if (nrow(cells) < 30L) {
    log_message("  SHM boxplot: <30 mutated cells; skipping.")
    return(invisible())
  }
  subj <- cells |>
    dplyr::group_by(.data[[subj_col]], Phenotype_2) |>
    dplyr::summarise(mean_shm = mean(SHM_total, na.rm = TRUE),
                     n_cells  = dplyr::n(),
                     .groups  = "drop") |>
    dplyr::filter(n_cells >= min_cells)
  if (nrow(subj) < 4L ||
      dplyr::n_distinct(subj$Phenotype_2) < 2L) {
    log_message("  SHM boxplot: too few subjects after >=", min_cells,
                "-cell filter; skipping.")
    return(invisible())
  }
  subj$Phenotype_2 <- factor(subj$Phenotype_2, levels = c("NIU", "Viral"))
  wp <- tryCatch(stats::wilcox.test(mean_shm ~ Phenotype_2,
                                    data = subj)$p.value,
                 error = function(e) NA_real_)
  n_by <- subj |>
    dplyr::count(Phenotype_2, name = "n") |>
    dplyr::arrange(Phenotype_2)
  sub_lab <- sprintf("Subject Wilcoxon p = %.3g  (n NIU=%d, n Viral=%d)",
                     wp,
                     n_by$n[n_by$Phenotype_2 == "NIU"]   %||% 0L,
                     n_by$n[n_by$Phenotype_2 == "Viral"] %||% 0L)
  p <- ggplot(subj, aes(x = Phenotype_2, y = mean_shm,
                        fill = Phenotype_2)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
    geom_jitter(width = 0.12, height = 0, size = 2.4, alpha = 0.85,
                color = "black") +
    scale_fill_manual(values = ETIOLOGY_GROUP_COLORS) +
    labs(title = "Heavy-chain SHM: NIU vs Viral (per subject)",
         subtitle = sub_lab,
         x = "Phenotype",
         y = "Mean heavy-chain SHM frequency (mutated cells)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  save_pdf_png(p,
    file.path(viz_subdir(paths, "repertoire"), "bcell_shm_boxplot_niu_vs_viral"),
    w = 5, h = 6)
}

# ---------------------------------------------------------------------------
# Cross-compartment B cell panels: PC1 partial coupling (B vs myeloid /
# B vs T cell) and LIANA LR heatmap per pair. Consume the CSVs written by
# R/57_bcell_cross_compartment_bridge.R + R/49_liana_bcell.R so the panels
# are recomputable from disk without re-reading any Seurat object.
#
# File stems (all under paths$viz_dir/11_pca_coupling = outputs/viz/eye/bcell/11_pca_coupling/):
#   bcell_fig6_pc1_partial_<partner>            primary scatter — partial r
#                                               residual-residual (matches
#                                               F4_F_partial style)
#   bcell_fig6_pc1_within_etiology_<partner>    within-arm faceted scatter,
#                                               honest within-patient view
#   bcell_fig6_pc1_per_substate_<partner>       per-(B substate x partner
#                                               substate) partial r heatmap
# ---------------------------------------------------------------------------

.bcell_cc_paths <- function(cfg) {
  cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/",
    viz    = "outputs/viz/cross_compartment/")
}

# Subject-level PC1 partial coupling — partial residual scatter. Mirrors the
# F4_F_partial panel produced for myeloid x T cell. partner in {"myeloid","tcell"}.
.bcell_fig6_pc1_partial <- function(paths, cfg, partner,
                                    weighting_choice = "weighted") {
  cc <- .bcell_cc_paths(cfg)
  scores_csv <- file.path(cc$tables,
                          paste0("pc1_bridge_bcell_", partner, "_scores.csv"))
  corr_csv   <- file.path(cc$tables,
                          paste0("pc1_bridge_bcell_", partner, "_correlation.csv"))
  if (!file.exists(scores_csv) || !file.exists(corr_csv)) {
    log_message("  bcell pc1_partial[", partner,
                "]: scores or correlation CSV missing; ",
                "run run_bcell_cross_compartment_bridge first.")
    return(invisible())
  }
  scores <- utils::read.csv(scores_csv, stringsAsFactors = FALSE)
  corr   <- utils::read.csv(corr_csv,   stringsAsFactors = FALSE)
  df <- scores[scores$weighting == weighting_choice, , drop = FALSE]
  if (nrow(df) < 3L) {
    log_message("  bcell pc1_partial[", partner,
                "]: <3 paired subjects; skipping.")
    return(invisible())
  }
  pal <- if (exists("ETIOLOGY_GROUP_COLORS", inherits = TRUE))
            ETIOLOGY_GROUP_COLORS else c(NIU = "#E21F26", Viral = "#397FB9")

  partial <- corr[corr$weighting == weighting_choice &
                  corr$stratum   == "partial_controlling_Phenotype_2",
                  , drop = FALSE]
  partial <- if (nrow(partial) > 0L) partial[1L, , drop = FALSE] else partial

  rx <- tryCatch(stats::resid(stats::lm(df$bcell_pc1   ~ df$Phenotype_2)),
                 error = function(e) NULL)
  ry <- tryCatch(stats::resid(stats::lm(df$partner_pc1 ~ df$Phenotype_2)),
                 error = function(e) NULL)
  if (is.null(rx) || is.null(ry)) {
    log_message("  bcell pc1_partial[", partner,
                "]: residual computation failed; skipping.")
    return(invisible())
  }
  df$resid_bcell   <- rx
  df$resid_partner <- ry

  partial_subtitle <- if (nrow(partial) > 0L && !is.na(partial$pearson_r))
    sprintf("Residuals after lm(PC1 ~ Phenotype_2).  Partial r=%.2f  95%% CI [%.2f, %.2f]  perm p=%.3g",
            partial$pearson_r, partial$pearson_ci_lo,
            partial$pearson_ci_hi, partial$permutation_p)
  else "Partial correlation (stats unavailable)"

  partner_title <- if (partner == "myeloid") "myeloid" else "T cell"
  p <- ggplot(df, aes(.data$resid_bcell, .data$resid_partner,
                      color = .data$Phenotype_2)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.3, color = "grey70") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6,
                formula = y ~ x,
                aes(group = 1), color = "black", alpha = 0.15) +
    geom_point(size = 4.5, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = .data$subject), size = 3.4,
                              box.padding = 0.35, max.overlaps = Inf,
                              show.legend = FALSE) +
    scale_color_manual(values = pal, name = NULL) +
    labs(title = paste0("Partial coupling: residual B cell vs ",
                        partner_title, " PC1"),
         subtitle = partial_subtitle,
         x = "B cell PC1 residual (after Phenotype_2)",
         y = paste0(partner_title, " PC1 residual (after Phenotype_2)")) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, lineheight = 1.2),
          aspect.ratio  = 1)
  save_pdf_png(p, file.path(viz_subdir(paths, "pca_coupling"),
                            paste0("bcell_fig6_pc1_partial_", partner)),
               w = 7.5, h = 7.5)
}

# Within-etiology faceted scatter — the honest within-patient view kept as
# a supplement to the partial residual panel. Mirrors F4_F_within_etiology.
.bcell_fig6_pc1_within_etiology <- function(paths, cfg, partner,
                                            weighting_choice = "weighted") {
  cc <- .bcell_cc_paths(cfg)
  scores_csv <- file.path(cc$tables,
                          paste0("pc1_bridge_bcell_", partner, "_scores.csv"))
  corr_csv   <- file.path(cc$tables,
                          paste0("pc1_bridge_bcell_", partner, "_correlation.csv"))
  if (!file.exists(scores_csv) || !file.exists(corr_csv)) return(invisible())
  scores <- utils::read.csv(scores_csv, stringsAsFactors = FALSE)
  corr   <- utils::read.csv(corr_csv,   stringsAsFactors = FALSE)
  df <- scores[scores$weighting == weighting_choice, , drop = FALSE]
  if (nrow(df) < 3L) return(invisible())
  pal <- if (exists("ETIOLOGY_GROUP_COLORS", inherits = TRUE))
            ETIOLOGY_GROUP_COLORS else c(NIU = "#E21F26", Viral = "#397FB9")
  by_eti <- corr[corr$weighting == weighting_choice &
                 grepl("^etiology:", corr$stratum), , drop = FALSE]
  if (nrow(by_eti) > 0L)
    by_eti$etiology <- sub("^etiology:", "", by_eti$stratum)
  eti_subtitle <- if (nrow(by_eti) > 0L)
    paste(vapply(seq_len(nrow(by_eti)), function(i) {
      sprintf("%s  n=%d  r=%.2f  95%% CI [%.2f, %.2f]  perm p=%.3g",
              by_eti$etiology[i], by_eti$n[i], by_eti$pearson_r[i],
              by_eti$pearson_ci_lo[i], by_eti$pearson_ci_hi[i],
              by_eti$permutation_p[i])
    }, character(1)), collapse = "\n")
  else ""
  partner_title <- if (partner == "myeloid") "myeloid" else "T cell"
  p <- ggplot(df, aes(.data$bcell_pc1, .data$partner_pc1,
                      color = .data$Phenotype_2)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.3, color = "grey70") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6,
                formula = y ~ x, alpha = 0.15) +
    geom_point(size = 2.6, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = .data$subject), size = 2.6,
                              box.padding = 0.25, max.overlaps = Inf,
                              show.legend = FALSE) +
    scale_color_manual(values = pal, name = NULL, guide = "none") +
    facet_wrap(~ .data$Phenotype_2, scales = "free", ncol = 2) +
    labs(title = paste0("Within-etiology coupling of B cell + ",
                        partner_title, " PC1"),
         subtitle = eti_subtitle,
         x = "B cell PC1 (NIU <-> Viral)",
         y = paste0(partner_title, " PC1 (NIU <-> Viral)")) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, lineheight = 1.2),
          strip.text    = element_text(face = "bold"),
          aspect.ratio  = 1)
  save_pdf_png(p, file.path(viz_subdir(paths, "pca_coupling"),
                            paste0("bcell_fig6_pc1_within_etiology_",
                                   partner)),
               w = 11, h = 6.5)
}

# Per-(B substate x partner substate) partial r heatmap. Mirrors
# F4_F_per_substate_heatmap. partner in {"myeloid","tcell"}.
.bcell_fig6_per_substate_heatmap <- function(paths, cfg, partner) {
  cc <- .bcell_cc_paths(cfg)
  per_sub_csv <- file.path(cc$tables,
                           paste0("pc1_bridge_bcell_", partner,
                                  "_per_substate.csv"))
  if (!file.exists(per_sub_csv)) {
    log_message("  bcell per_substate[", partner, "]: ", per_sub_csv,
                " missing.")
    return(invisible())
  }
  ps <- utils::read.csv(per_sub_csv, stringsAsFactors = FALSE)
  if (nrow(ps) == 0L) return(invisible())

  ps$bcell_display   <- vapply(as.character(ps$bcell_substate),
                               function(id) get_substate_display(cfg, "bcell", id),
                               character(1))
  ps$partner_display <- vapply(as.character(ps$partner_substate),
                               function(id) get_substate_display(cfg, partner, id),
                               character(1))
  ps$bcell_display   <- factor(ps$bcell_display,
                               levels = sort(unique(ps$bcell_display)))
  ps$partner_display <- factor(ps$partner_display,
                               levels = sort(unique(ps$partner_display)))
  ps$sig <- !is.na(ps$partial_permutation_p) & ps$partial_permutation_p < 0.05
  ps$label <- ifelse(
    is.na(ps$partial_pearson_r), "",
    sprintf("%s%.2f\n(n=%d)",
            ifelse(ps$sig, "* ", ""),
            ps$partial_pearson_r, ps$n))
  lim <- max(abs(ps$partial_pearson_r), na.rm = TRUE)
  lim <- if (is.finite(lim)) lim else 1
  ps_sig <- dplyr::filter(ps, .data$sig)
  partner_title <- if (partner == "myeloid") "myeloid" else "T cell"

  p <- ggplot(ps, aes(.data$partner_display, .data$bcell_display,
                      fill = .data$partial_pearson_r)) +
    geom_tile(color = "grey95", linewidth = 0.4) +
    geom_tile(data = ps_sig, fill = NA, color = "black", linewidth = 1.0) +
    geom_text(aes(label = .data$label,
                  fontface = ifelse(.data$sig, "bold", "plain")),
              size = 2.8) +
    viridis::scale_fill_viridis(option = "viridis",
                                limits = c(-lim, lim),
                                name = "Partial Pearson r\n(within etiology)") +
    coord_fixed(ratio = 1) +
    labs(title = paste0("Per-substate B cell vs ", partner_title,
                        " PC1 coupling"),
         subtitle = paste0("Partial r controlling Phenotype_2. * marks ",
                           "permutation p < 0.05. ",
                           sum(ps$sig, na.rm = TRUE), " of ",
                           nrow(ps), " cells significant."),
         x = paste0(partner_title, " substate"),
         y = "B cell substate") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9),
          panel.grid    = element_blank(),
          axis.text.x   = element_text(angle = 30, hjust = 1))
  n_x <- length(levels(ps$partner_display))
  n_y <- length(levels(ps$bcell_display))
  cell_in <- 0.85
  save_pdf_png(p,
               file.path(viz_subdir(paths, "pca_coupling"),
                         paste0("bcell_fig6_pc1_per_substate_", partner)),
               w = n_x * cell_in + 4.5,
               h = n_y * cell_in + 3.5)
}

# Per-substate pseudobulk PCA scatter (PC1 x PC2) with NIU + Viral 95%
# ellipses. Mirrors R/88_viz_tcell.R .tcell_F4_pc_scatter (F4_D panel) for
# the B cell compartment. Substates that didn't reach FDR<0.05 on the Welch
# t are tagged "(underpowered)" in the strip label so the panel stays
# honest about which substates actually carry the signal.
.bcell_fig6_pc_scatter <- function(paths, cfg) {
  scores_csv <- file.path(paths$results_tables, "pca_subject_scores.csv")
  sig_csv    <- file.path(paths$results_tables, "pca_pc1_significance.csv")
  if (!file.exists(scores_csv) || !file.exists(sig_csv)) {
    log_message("  bcell pc_scatter: pca_subject_scores or _pc1_significance ",
                "CSV missing; skipping.")
    return(invisible())
  }
  scores <- utils::read.csv(scores_csv, stringsAsFactors = FALSE)
  sig    <- utils::read.csv(sig_csv,    stringsAsFactors = FALSE)
  if (nrow(scores) == 0L) return(invisible())

  sig$separating <- as.logical(sig$separating)
  scores <- dplyr::left_join(scores,
                             sig |> dplyr::select("substate", "separating",
                                                  "q_value"),
                             by = "substate")
  scores$substate_label <- vapply(as.character(scores$substate),
                                  function(id) get_substate_display(cfg, "bcell", id),
                                  character(1))
  scores$substate_label <- ifelse(!is.na(scores$separating) & scores$separating,
                                  scores$substate_label,
                                  paste0(scores$substate_label, "\n(underpowered)"))
  scores$substate_label <- factor(scores$substate_label,
                                  levels = unique(scores$substate_label))

  pal <- if (exists("ETIOLOGY_GROUP_COLORS", inherits = TRUE))
            ETIOLOGY_GROUP_COLORS else c(NIU = "#E21F26", Viral = "#397FB9")

  p <- ggplot(scores, aes(.data$PC1_oriented, .data$PC2,
                          color = .data$Phenotype_2)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.25, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.25, color = "grey70") +
    geom_point(size = 2.2, alpha = 0.9) +
    stat_ellipse(level = 0.95, linewidth = 0.5,
                 aes(group = .data$Phenotype_2)) +
    scale_color_manual(values = pal, name = NULL) +
    facet_wrap(~ .data$substate_label, scales = "free", ncol = 2) +
    labs(title = "B cell per-substate pseudobulk PCA",
         subtitle = "PC1 oriented so Viral centroid is positive",
         x = "PC1", y = "PC2") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          aspect.ratio  = 1,
          panel.grid.minor = element_blank())

  n_facets <- length(unique(scores$substate))
  save_pdf_png(p, file.path(viz_subdir(paths, "pca_coupling"), "bcell_fig6_pc12_disease"),
               w = 8.5,
               h = 3.4 * ceiling(n_facets / 2) + 2)
}

# Single-shot dispatcher — both partners x all four panels. Idempotent and
# tryCatch-wrapped so a failure in one panel doesn't bring the rest down.
.bcell_fig6_cross_compartment_block <- function(paths, cfg) {
  tryCatch(.bcell_fig6_pc_scatter(paths, cfg),
           error = function(e)
             log_message("  bcell pc_scatter failed: ",
                         conditionMessage(e)))
  for (partner in c("myeloid", "tcell")) {
    tryCatch(.bcell_fig6_pc1_partial(paths, cfg, partner),
             error = function(e)
               log_message("  bcell pc1_partial[", partner, "] failed: ",
                           conditionMessage(e)))
    tryCatch(.bcell_fig6_pc1_within_etiology(paths, cfg, partner),
             error = function(e)
               log_message("  bcell pc1_within[", partner, "] failed: ",
                           conditionMessage(e)))
    tryCatch(.bcell_fig6_per_substate_heatmap(paths, cfg, partner),
             error = function(e)
               log_message("  bcell pc1_per_substate[", partner, "] failed: ",
                           conditionMessage(e)))
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
run_visualizations_bcell <- function(cfg) {
  paths <- get_target_paths(cfg, "bcell")
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Bcell IntegratedSeuratObject.rds not found. Skipping bcell viz.")
    return(invisible(TRUE))
  }
  log_message("=== B/plasma compartment visualizations ===")
  ensure_dir(paths$viz_dir)
  obj <- readRDS(obj_path)

  # --- Shared cross-cutting
  viz_compartment_umap(obj, "bcell", paths, cfg)
  viz_compartment_umap_etiology(obj, "bcell", paths, cfg)
  viz_compartment_dotplot(obj, "bcell", paths, cfg)
  viz_compartment_milo_da_box("bcell", paths, cfg)
  viz_compartment_composition(obj, "bcell", paths, cfg)
  viz_compartment_volcano("bcell", paths, cfg)
  viz_compartment_volcano_per_substate("bcell", paths, cfg)
  viz_compartment_niu_subcontrast_heatmap("bcell", paths, cfg)
  viz_compartment_functional_gestalt_full_heatmap(obj, "bcell", paths, cfg)

  # --- B-cell-specific panels
  .bcell_isotype_by_substate(obj, paths, cfg)
  .bcell_shm_boxplot_niu_vs_viral(obj, paths, cfg)

  # --- Figure 6 panels (local antigen-driven B cell response)
  ensure_dir("outputs/tables/eye/bcell")
  # Cross-compartment PC1 coupling + LIANA LR heatmaps (B <-> myeloid,
  # B <-> T cell). Consumes outputs/tables/cross_compartment/pc1_bridge_bcell_*
  # written by R/57_bcell_cross_compartment_bridge.R.
  tryCatch(.bcell_fig6_cross_compartment_block(paths, cfg),
           error = function(e) log_message("  fig6 cross-compartment failed: ",
                                           conditionMessage(e)))

  # --- Full MiloR block
  if (exists("viz_milo")) {
    tryCatch(viz_milo(obj, cfg, paths, target = "bcell"),
             error = function(e)
               log_message("  viz_milo failed for bcell: ",
                           conditionMessage(e)))
  }

  log_message("=== B/plasma visualizations complete ===")
  invisible(TRUE)
}

# Standalone Figure 6 dispatcher. Lets the pipeline render just the Figure 6
# panels (cover, substate-by-tissue, gini, top clones, encoder benchmark,
# LIANA dotplot) without re-running the heavy cross-cutting + MiloR blocks
# from run_visualizations_bcell(). Used by run_pipeline.R Phase 6.
run_visualizations_bcell_fig6 <- function(cfg) {
  paths <- get_target_paths(cfg, "bcell")
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Bcell object missing at ", obj_path, "; skipping Fig6 viz.")
    return(invisible(FALSE))
  }
  ensure_dir(paths$viz_dir); ensure_dir("outputs/tables/eye/bcell")
  log_message("=== Figure 6 (B cell) — biology panels ===")
  obj <- readRDS(obj_path)

  # Panel C: isotype distribution per substate, NIU vs Viral
  tryCatch(.bcell_isotype_by_substate(obj, paths, cfg),
           error = function(e) log_message("  Panel C: ",
                                            conditionMessage(e)))
  # Panels H-K: cross-compartment PC1 coupling + LIANA LR heatmaps for
  # B<->myeloid and B<->T cell. Four sub-panels per pair:
  #   pc1_partial, pc1_within_etiology, pc1_per_substate, liana_heatmap.
  tryCatch(.bcell_fig6_cross_compartment_block(paths, cfg),
           error = function(e) log_message("  Panels H-K: ",
                                            conditionMessage(e)))

  log_message("=== Figure 6 (B cell) complete ===")
  invisible(TRUE)
}

# run_visualizations_bcell renders the B/plasma compartment figure block:
# cross-cutting cluster panels (substate UMAP labeled+stripped, UMAP by
# etiology, marker dot plot, NIU-vs-Viral composition, compartment-global
# and per-substate volcanoes, GSEA pathway-by-substate heatmap, NIU-vs-Viral
# top pathway bar, pathway-by-etiology heatmap, NIU sub-contrast heatmap)
# plus adaptive-receptor panels (CDR3-H length NIU vs Viral, IGHV usage
# heatmap NIU vs Viral, public-clone count at multiple sharing thresholds),
# and B-cell-specific blocks: eye-blood BCR overlap boxplot, SHM violin per
# substate, SHM NIU-vs-Viral split by substate, isotype distribution per
# substate split by NIU/Viral. The full MiloR viz block runs last so the
# bcell milo nhood-graph, beeswarm, boxplot, summary, composition, volcano,
# and nhood-size hist all land under outputs/viz/eye/bcell/07_milo/. All
# filenames are snake_case descriptions; no F4 panel-letter scheme is used.
