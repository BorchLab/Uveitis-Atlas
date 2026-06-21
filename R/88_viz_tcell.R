# R/88_viz_tcell.R
# T cell compartment figure block. Uses shared helpers in
# R/81_viz_compartment_helpers.R for cluster panels, GSEA heatmap, pathway
# bars, etiology heatmaps, public-clone count, CDR3 length, TRBV usage,
# per-substate volcano grid, and the NIU sub-contrast heatmap.
#
# T cell-specific panels:
#   - Eye-blood TCR overlap boxplot
#   - HLA-B27-stratified TCR features: V-gene usage and CDR3 length within
#     HLA-B27+ vs other NIU vs Viral (B27 has n=10 in NIU, large enough)
#   - Alluvial of top expanded eye clones eye -> blood substate
#   - UCell UMAP overlay for top GSEA pathways
#
# All filenames use snake_case descriptions. The full MiloR viz block from
# 82_viz_dispatch.R runs at the end.

suppressPackageStartupMessages({
  library(Seurat)
  library(scplotter)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  tryCatch(install.packages("ggalluvial"), error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# Eye-blood TCR overlap boxplot (Wilcoxon NIU vs Viral on per-subject
# fraction of eye clones detected in blood)
# ---------------------------------------------------------------------------
.tcell_eye_blood_overlap <- function(cfg) {
  paths <- get_target_paths(cfg, "tcell")
  ov_path <- "outputs/tables/repertoire/TCR_eye_blood_overlap.csv"
  if (!file.exists(ov_path)) return(invisible())
  ov <- read.csv(ov_path, stringsAsFactors = FALSE)
  group_col <- if ("phenotype" %in% colnames(ov)) "phenotype" else "etiology"
  ov <- ov |> dplyr::filter(!is.na(frac_eye_shared_with_blood),
                            .data[[group_col]] %in% c("NIU", "Viral"))
  if (nrow(ov) == 0) return(invisible())
  wp <- tryCatch(wilcox.test(as.formula(paste("frac_eye_shared_with_blood ~",
                                              group_col)),
                              data = ov)$p.value,
                 error = function(e) NA_real_)
  p <- ggplot(ov, aes(x = .data[[group_col]],
                      y = frac_eye_shared_with_blood,
                      fill = .data[[group_col]])) +
       geom_boxplot(outlier.shape = NA, alpha = 0.5) +
       geom_jitter(width = 0.1, alpha = 0.7) +
       scale_fill_manual(values = ETIOLOGY_GROUP_COLORS) +
       labs(title = "Fraction of eye TCR clones detected in blood",
            subtitle = sprintf("Wilcoxon p = %.3g", wp),
            x = "Phenotype", y = "Fraction shared (eye -> blood)") +
       theme_classic() +
       theme(legend.position = "none",
             plot.title = element_text(face = "bold"))
  save_pdf_png(p, file.path(viz_subdir(paths, "repertoire"), "tcell_eye_blood_tcr_overlap"),
               w = 6, h = 6)
}

# ---------------------------------------------------------------------------
# F4 panels D / E: T cell pseudobulk PCA (mirrors F3 panel F + E framing).
# Both helpers read the on-disk CSVs that R/45_compartment_pca writes; they
# are no-ops when those CSVs are missing so the rest of the T cell viz block
# still runs.
# ---------------------------------------------------------------------------

# Panel D — Per-substate PC1 vs PC2 scatter, NIU (red) + Viral (blue) 95%
# ellipses. Substates that didn't reach FDR<0.05 on the Welch t are greyed
# out with an "(underpowered)" tag in the strip label so the figure stays
# honest about which substates carry the signal.
.tcell_F4_pc_scatter <- function(paths, cfg) {
  if (!requireNamespace("ggh4x", quietly = TRUE)) return(invisible())
  scores_csv <- file.path(paths$results_tables, "pca_subject_scores.csv")
  sig_csv    <- file.path(paths$results_tables, "pca_pc1_significance.csv")
  if (!file.exists(scores_csv) || !file.exists(sig_csv)) {
    log_message("  F4_D: PCA score / significance CSVs missing for tcell; ",
                "skipping panel D.")
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
                                  function(id) get_substate_display(cfg, "tcell", id),
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
    labs(title = "T cell per-substate pseudobulk PCA",
         subtitle = "PC1 oriented so Viral centroid is positive",
         x = "PC1", y = "PC2") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          aspect.ratio  = 1,
          panel.grid.minor = element_blank())

  out_dir <- file.path(viz_subdir(paths, "pca_coupling"), "F4_tcell")
  ensure_dir(out_dir)
  n_facets <- length(unique(scores$substate))
  save_pdf_png(p, file.path(out_dir, "F4_D_pc12_disease"),
               w = 8.5,
               h = 3.4 * ceiling(n_facets / 2) + 2)
}

# Panel E — Lollipop of per-program PC1 loadings (TCR_signal, CD4_CD8,
# Costim_receptor, Checkpoint, Cytokine_receptor, Effector_cytokine,
# Chemokine_receptor, Tissue_residence, Transcription_factor). Restricted
# to separating substates by construction of pc1_loadings_by_program.csv.
.tcell_F4_pc1_loadings_by_program <- function(paths, cfg) {
  prog_csv <- file.path(paths$results_tables, "pc1_loadings_by_program.csv")
  if (!file.exists(prog_csv)) {
    log_message("  F4_E: pc1_loadings_by_program.csv missing for tcell; ",
                "skipping panel E.")
    return(invisible())
  }
  df <- utils::read.csv(prog_csv, stringsAsFactors = FALSE)
  if (nrow(df) == 0L) {
    log_message("  F4_E: pc1_loadings_by_program empty (no separating ",
                "substates?); skipping.")
    return(invisible())
  }
  df$substate_label <- vapply(as.character(df$substate),
                              function(id) get_substate_display(cfg, "tcell", id),
                              character(1))
  df$direction <- ifelse(df$PC1_oriented >= 0, "Viral-driving", "NIU-driving")

  pal <- if (exists("ETIOLOGY_GROUP_COLORS", inherits = TRUE))
            c(`NIU-driving`   = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
              `Viral-driving` = unname(ETIOLOGY_GROUP_COLORS["Viral"]))
        else c(`NIU-driving` = "#E21F26", `Viral-driving` = "#397FB9")

  # Per-(substate, program) gene ordering by signed loading so the lollipop
  # reads top -> bottom within each facet.
  df <- df |>
    dplyr::arrange(.data$substate_label, .data$program, .data$PC1_oriented) |>
    dplyr::mutate(gene_facet = factor(.data$gene,
                                      levels = unique(.data$gene)))

  p <- ggplot(df, aes(.data$PC1_oriented, .data$gene_facet,
                      color = .data$direction)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_segment(aes(x = 0, xend = .data$PC1_oriented,
                     y = .data$gene_facet, yend = .data$gene_facet),
                 linewidth = 0.4) +
    geom_point(size = 2.2) +
    scale_color_manual(values = pal, name = NULL) +
    facet_grid(.data$program ~ .data$substate_label, scales = "free_y",
               space = "free_y") +
    labs(title = "T cell PC1 loadings on canonical receptor / effector programs",
         subtitle = "Restricted to disease-separating substates (q < 0.05)",
         x = "PC1 loading (NIU <-> Viral)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          strip.text.y  = element_text(angle = 0),
          axis.text.y   = element_text(size = 8),
          legend.position = "bottom",
          panel.grid.minor = element_blank())

  out_dir <- file.path(viz_subdir(paths, "pca_coupling"), "F4_tcell")
  ensure_dir(out_dir)
  n_substates <- length(unique(df$substate))
  n_programs  <- length(unique(df$program))
  save_pdf_png(p, file.path(out_dir, "F4_E_pc1_loadings"),
               w = max(8, n_substates * 2.4 + 3),
               h = max(7, n_programs * 1.1 + 2))
}

# ---------------------------------------------------------------------------
# Paired eye<->blood Gini panel (T cell analog of .bcell_fig6_gini).
# Reads outputs/tables/eye/tcell/tcell_paired_eye_blood_metrics.csv produced
# by R/55_tcell_paired_eye_blood_metrics.R. Faceted by Phenotype_2 with
# per-etiology paired Wilcoxon annotated onto each facet.
# ---------------------------------------------------------------------------
.tcell_paired_metrics <- function() {
  csv_path <- "outputs/tables/eye/tcell/tcell_paired_eye_blood_metrics.csv"
  if (!file.exists(csv_path)) return(NULL)
  d <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  d |>
    dplyr::group_by(Subject) |>
    dplyr::filter(all(c("Eye","Blood") %in% Tissue_1)) |>
    dplyr::ungroup()
}

.tcell_fig6_gini <- function(paths, cfg) {
  d <- .tcell_paired_metrics()
  if (is.null(d) || nrow(d) < 4) {
    log_message("  fig6 tcell gini: paired-metrics CSV missing or <2 paired ",
                "subjects; run cfg$steps_fig6$tcell_paired_metrics first.")
    return(invisible())
  }
  .paired_p <- function(sub) {
    tryCatch({
      sub <- sub |>
        dplyr::filter(Tissue_1 %in% c("Eye","Blood"), !is.na(gini)) |>
        dplyr::group_by(Subject) |>
        dplyr::filter(dplyr::n_distinct(Tissue_1) == 2L) |>
        dplyr::ungroup()
      if (dplyr::n_distinct(sub$Subject) < 2L) return(NA_real_)
      sub <- sub |> dplyr::arrange(Subject, Tissue_1)
      stats::wilcox.test(sub$gini[sub$Tissue_1 == "Eye"],
                         sub$gini[sub$Tissue_1 == "Blood"],
                         paired = TRUE)$p.value
    }, error = function(e) NA_real_)
  }
  stats_df <- d |>
    dplyr::group_by(Phenotype_2) |>
    dplyr::group_modify(~ tibble::tibble(
      p      = .paired_p(.x),
      n_subj = dplyr::n_distinct(.x$Subject[.x$Tissue_1 %in% c("Eye","Blood")])
    )) |>
    dplyr::ungroup() |>
    dplyr::mutate(label = sprintf("Wilcoxon p = %.3g\nn = %d subjects",
                                  p, n_subj))
  p <- ggplot(d, aes(x = Tissue_1, y = gini, group = Subject)) +
    geom_line(aes(group = Subject), color = "grey60", alpha = 0.6) +
    geom_point(aes(color = Phenotype_2), size = 2.2) +
    geom_text(data = stats_df,
              aes(x = 1.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.3, size = 3,
              lineheight = 0.95) +
    facet_wrap(~ Phenotype_2, nrow = 1) +
    scale_color_manual(values = ETIOLOGY_GROUP_COLORS) +
    labs(title = "TCR Gini: eye vs blood",
         x = NULL, y = "Gini (clone-size)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "none")
  save_pdf_png(p, file.path(viz_subdir(paths, "repertoire"), "tcell_fig6_gini_eye_vs_blood"),
               w = 4, h = 5)
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
run_visualizations_tcell <- function(cfg) {
  paths <- get_target_paths(cfg, "tcell")
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Tcell IntegratedSeuratObject.rds not found. Skipping tcell viz.")
    return(invisible(TRUE))
  }
  log_message("=== T cell compartment visualizations ===")
  ensure_dir(paths$viz_dir)
  obj <- readRDS(obj_path)

  # --- Shared cross-cutting
  viz_compartment_umap(obj, "tcell", paths, cfg)
  viz_compartment_umap_etiology(obj, "tcell", paths, cfg)
  viz_compartment_dotplot(obj, "tcell", paths, cfg)

  # Canonical T-cell markers on UMAP (12-gene feature plot, 4x3 grid).
  # Rows are coherent biological axes:
  #   row 1 - lineage + naive/CM:   CD4, CD8A, IL7R, TCF7
  #   row 2 - effector to exhaustion: GZMK, GZMB, CXCL13, PDCD1
  #   row 3 - TRM, Treg, cycling:   ITGA1, CXCR6, FOXP3, MKI67
  tcell_markers <- c("CD4",   "CD8A",  "IL7R",  "TCF7",
                     "GZMK",  "GZMB",  "CXCL13","PDCD1",
                     "ITGA1", "CXCR6", "FOXP3", "MKI67")
  tcell_markers <- intersect(tcell_markers, rownames(obj))
  if (length(tcell_markers) > 0) {
    safe_plot(bquote(
      FeatureStatPlot(obj,
                      raster    = FALSE,
                      features  = .(tcell_markers),
                      plot_type = "dim",
                      reduction = "UMAP",
                      ncol      = 4,
                      palette   = "viridis",
                      bg_cutoff = -Inf,
                      hex       = TRUE) &
        theme(axis.title = element_blank(),
              axis.text  = element_blank(),
              axis.ticks = element_blank())
    ), file.path(viz_subdir(paths, "markers"), "tcell_canonical_markers_featureplot.pdf"),
       width = 16, height = 10)
  }

  viz_compartment_milo_da_box("tcell", paths, cfg)
  viz_compartment_composition(obj, "tcell", paths, cfg)
  viz_compartment_volcano("tcell", paths, cfg)
  viz_compartment_volcano_per_substate("tcell", paths, cfg)
  viz_compartment_functional_gestalt_full_heatmap(obj, "tcell", paths, cfg)

  # --- T-cell-specific panels
  tryCatch(.tcell_fig6_gini(paths, cfg),
           error = function(e) log_message("  tcell fig6 gini failed: ",
                                           conditionMessage(e)))
  .tcell_eye_blood_overlap(cfg)

  # --- F4 panels D + E: PC1 vs PC2 ellipses and PC1 loadings by program
  .tcell_F4_pc_scatter(paths, cfg)
  .tcell_F4_pc1_loadings_by_program(paths, cfg)

  # --- Viral vs NIU intraocular TCR comparison (immGLIPH + immLynx)
  if (isTRUE(cfg$steps$viz_tcr_compare)) {
    # scRepertoire descriptive panels (homeostasis, proportion, overlap,
    # StartracDiversity, circlize chord) -- general clonal-structure views
    # that complement the focused Viral-vs-NIU enrichment plots below.
    .tcell_screp_homeostasis(obj, paths, cfg)
    .tcell_screp_proportion(obj, paths, cfg)
    .tcell_screp_overlap(obj, paths, cfg)
    .tcell_screp_startrac(obj, paths, cfg)
    .tcell_screp_circlize(obj, paths, cfg)
    # GLIPH (motif convergence)
    .tcell_gliph_enrichment_volcano(paths, cfg)
    .tcell_gliph_motif_heatmap(paths, cfg)
    .tcell_gliph_motif_logos(paths, cfg)
    .tcell_gliph_subject_heatmap(paths, cfg)
    .tcell_gliph_network(paths, cfg)
    # tcrdist (sequence distance)
    .tcell_tcrdist_umap(paths, cfg)
    .tcell_tcrdist_neighborhood(paths, cfg)
    # OLGA / SoNNia (generation prob + selection)
    .tcell_olga_pgen_ridges(paths, cfg)
    .tcell_pgen_vs_clone_size(paths, cfg)
    .tcell_pgen_vs_clone_size_gliph(paths, cfg)
    .tcell_sonnia_volcano(paths, cfg)
    # Repertoire structure
    .tcell_public_private_clones(paths, cfg)
    .tcell_clonality_curve(paths, cfg)
    # HLA-B27 pathogenic signature
    .tcell_hla_b27_pathogenic_summary(paths, cfg)
    .tcell_hla_b27_trbv_usage(paths, cfg)
    .tcell_gliph_b27_contingency(paths, cfg)
    # UMAP overlays on gene-expression space (main RNA-based UMAP)
    .tcell_umap_hla_b27(obj, paths, cfg)
    .tcell_umap_gliph_clusters(obj, paths, cfg)
    .tcell_umap_pgen(obj, paths, cfg)
    .tcell_umap_clone_size(obj, paths, cfg)
    # Bar plots of the same highlighted populations across substate /
    # phenotype / individual etiology
    .tcell_bar_hla_b27(obj, paths, cfg)
    .tcell_bar_gliph_directional(obj, paths, cfg)
    .tcell_pgen_by_substate(obj, paths, cfg)
    # Cross-method summary table
    .tcell_cross_method_summary(paths, cfg)
  }

  # --- Full MiloR block
  if (exists("viz_milo")) {
    tryCatch(viz_milo(obj, cfg, paths, target = "tcell"),
             error = function(e)
               log_message("  viz_milo failed for tcell: ",
                           conditionMessage(e)))
  }

  log_message("=== T cell visualizations complete ===")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Viral vs NIU intraocular TCR comparison panels (immGLIPH + immLynx)
# Reads outputs produced by R/64_immgliph.R, R/65_immlynx_tcrdist.R,
# R/66_immlynx_olga_sonnia.R. Headline plots use Phenotype_2 (Viral/NIU)
# with ETIOLOGY_GROUP_COLORS; secondary facets break out individual
# Etiology subtypes using ETIOLOGY_SUBTYPE_COLORS so the phenotype family
# (red = NIU, blue = Viral) stays readable at a glance.
#
# Convention: theme_classic(base_size = 10), save_pdf_png() for output,
# bold plot.title, italic plot.subtitle. Diverging fills use
# scale_fill_gradient2(low = NIU, mid = "white", high = Viral) when the
# signed metric is log2(Viral / NIU) so positive => Viral-enriched (blue).
# ---------------------------------------------------------------------------

# Internal: dispatch a per-etiology color scale that falls back gracefully
# if a subtype isn't covered by ETIOLOGY_SUBTYPE_COLORS.
.etiology_subtype_scale <- function(levels, aesthetic = c("color", "fill")) {
  aesthetic <- match.arg(aesthetic)
  pal <- ETIOLOGY_SUBTYPE_COLORS
  missing <- setdiff(levels, names(pal))
  if (length(missing) > 0) {
    extras <- setNames(rep("grey60", length(missing)), missing)
    pal <- c(pal, extras)
  }
  if (aesthetic == "color")
    ggplot2::scale_color_manual(values = pal, name = "Etiology")
  else
    ggplot2::scale_fill_manual(values = pal, name = "Etiology")
}

.tcr_compare_theme <- function() {
  theme_classic(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic", size = 9,
                                       color = "grey25"),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text       = element_text(face = "bold", size = 9))
}

# ---- Panel 1: GLIPH enrichment volcano -----------------------------------
# log2 OR vs -log10 FDR. Replaces the earlier bar-style "heatmap" because a
# volcano makes effect size + significance visible in one plot. Top hits on
# each side labeled.
.tcell_gliph_enrichment_volcano <- function(paths, cfg) {
  enrich_csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                          "repertoire", "gliph_enrichment_viral_vs_niu.csv")
  if (!file.exists(enrich_csv)) {
    log_message("  TCR-compare: gliph_enrichment CSV missing; skipping volcano.")
    return(invisible())
  }
  df <- utils::read.csv(enrich_csv, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible())

  df <- df |>
    dplyr::mutate(
      log2_OR  = log2(pmax(median_OR, 1e-3)),
      nlog_FDR = -log10(pmax(FDR, 1e-300)),
      direction = dplyr::case_when(
        FDR < 0.05 & log2_OR >  0 ~ "Viral_up",
        FDR < 0.05 & log2_OR <  0 ~ "NIU_up",
        TRUE                       ~ "ns"
      ),
      label_txt = ifelse(!is.na(motif) & nchar(motif) > 0,
                         sprintf("%s [%s]", cluster_id, motif),
                         as.character(cluster_id))
    )

  top_each_side <- df |>
    dplyr::filter(direction != "ns") |>
    dplyr::group_by(direction) |>
    dplyr::slice_min(FDR, n = 8, with_ties = FALSE) |>
    dplyr::ungroup()

  p <- ggplot(df, aes(x = log2_OR, y = nlog_FDR, color = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(alpha = 0.8, size = 1.8) +
    ggrepel::geom_text_repel(data = top_each_side,
                             aes(label = label_txt),
                             size = 2.6, max.overlaps = Inf,
                             min.segment.length = 0.1,
                             show.legend = FALSE) +
    scale_color_manual(values = ETIOLOGY_DIRECTION_COLORS,
                       labels = c(NIU_up = "NIU-enriched",
                                  Viral_up = "Viral-enriched",
                                  ns = "n.s."),
                       name = NULL) +
    labs(title = "GLIPH convergence groups: Viral vs NIU enrichment",
         subtitle = sprintf("Bootstrap median OR; FDR < 0.05 dashed; n = %d groups",
                            nrow(df)),
         x = "log2 OR (Viral / NIU)",
         y = "-log10 FDR") +
    .tcr_compare_theme() +
    theme(legend.position = "top")
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_gliph_enrichment_volcano"),
               w = 8, h = 6.5)
}

# ---- Panel 2: GLIPH significant-motif lollipop (companion to the volcano)
# Filtered to FDR < 0.05 only. Lollipop point size = -log10 FDR, color =
# direction (Viral-enriched = blue, NIU-enriched = red), x = log2 OR.
.tcell_gliph_motif_heatmap <- function(paths, cfg) {
  enrich_csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                          "repertoire", "gliph_enrichment_viral_vs_niu.csv")
  if (!file.exists(enrich_csv)) return(invisible())
  df <- utils::read.csv(enrich_csv, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible())

  sig_alpha <- cfg$gliph$fdr_cutoff %||% 0.05
  sig <- df |>
    dplyr::filter(!is.na(FDR), FDR < sig_alpha) |>
    dplyr::mutate(
      label     = ifelse(!is.na(motif) & nchar(motif) > 0,
                         sprintf("%s [%s]", cluster_id, motif),
                         as.character(cluster_id)),
      log2_OR   = log2(pmax(median_OR, 1e-3)),
      nlog_FDR  = -log10(pmax(FDR, 1e-300)),
      direction = ifelse(log2_OR > 0, "Viral_up", "NIU_up")
    ) |>
    dplyr::arrange(log2_OR)
  if (nrow(sig) == 0) {
    log_message("  GLIPH lollipop: no groups passed FDR < ", sig_alpha,
                "; skipping.")
    return(invisible())
  }
  sig$label <- factor(sig$label, levels = sig$label)

  p <- ggplot(sig, aes(x = log2_OR, y = label, color = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_segment(aes(x = 0, xend = log2_OR,
                     y = label, yend = label),
                 linewidth = 0.45) +
    geom_point(aes(size = nlog_FDR)) +
    scale_color_manual(values = ETIOLOGY_DIRECTION_COLORS,
                       labels = c(NIU_up   = "NIU-enriched",
                                  Viral_up = "Viral-enriched"),
                       name = NULL) +
    scale_size_continuous(range = c(2, 7),
                          name  = expression(-log[10]~FDR)) +
    labs(title = sprintf("Significant GLIPH groups (FDR < %.2f)", sig_alpha),
         subtitle = sprintf("%d groups; point size = -log10 FDR; x = bootstrap log2 OR",
                            nrow(sig)),
         x = "log2 OR (Viral / NIU)", y = NULL) +
    .tcr_compare_theme() +
    theme(axis.text.y     = element_text(size = 7.5),
          legend.position = "right")
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_gliph_significant_lollipop"),
               w = 8.5, h = max(5, 0.26 * nrow(sig) + 2))
}

# ---- Panel 3: Motif logos -------------------------------------------------
.tcell_gliph_motif_logos <- function(paths, cfg) {
  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    log_message("  TCR-compare: ggseqlogo not available; skipping logos.")
    return(invisible())
  }
  res_path <- file.path(get_target_paths(cfg, "all")$results_objects,
                        "ImmGLIPHResults.rds")
  if (!file.exists(res_path)) return(invisible())
  res <- readRDS(res_path)

  clusters <- res$clusters
  enrich   <- res$enrich
  if (is.null(clusters) || is.null(enrich)) return(invisible())

  top_viral <- enrich |>
    dplyr::filter(median_OR > 1, FDR < 0.1) |>
    dplyr::slice_min(FDR, n = 6, with_ties = FALSE) |>
    dplyr::pull(cluster_id)
  top_niu <- enrich |>
    dplyr::filter(median_OR < 1, FDR < 0.1) |>
    dplyr::slice_min(FDR, n = 6, with_ties = FALSE) |>
    dplyr::pull(cluster_id)

  picks <- c(top_viral, top_niu)
  if (length(picks) == 0) return(invisible())

  seq_lists <- lapply(picks, function(cid) {
    s <- clusters$CDR3b[clusters$cluster_id == cid]
    s <- s[!is.na(s) & nchar(s) > 0]
    if (length(s) < 3) return(NULL)
    L <- as.integer(names(sort(table(nchar(s)), decreasing = TRUE))[1])
    s <- s[nchar(s) == L]
    if (length(s) < 3) return(NULL)
    s
  })
  names(seq_lists) <- vapply(picks, function(cid) {
    side <- if (cid %in% top_viral) "Viral" else "NIU"
    sprintf("%s | %s", side, cid)
  }, character(1))
  seq_lists <- seq_lists[!vapply(seq_lists, is.null, logical(1))]
  if (length(seq_lists) == 0) return(invisible())

  p <- ggseqlogo::ggseqlogo(seq_lists, ncol = 3) +
    labs(title = "Top GLIPH motifs: Viral- and NIU-enriched groups",
         subtitle = "Length-matched CDR3-beta within each convergence group") +
    .tcr_compare_theme() +
    theme(strip.text = element_text(size = 8))
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_gliph_motif_logos"),
               w = 12, h = max(6, length(seq_lists) * 0.9))
}

# ---- Panel 4: Per-subject x GLIPH cluster heatmap -------------------------
# Catches the case where one outlier subject drives the entire signal.
# Rows = top GLIPH groups, columns = subjects (ordered NIU then Viral),
# fill = log1p(n CDR3s from this subject in this group).
.tcell_gliph_subject_heatmap <- function(paths, cfg) {
  res_path <- file.path(get_target_paths(cfg, "all")$results_objects,
                        "ImmGLIPHResults.rds")
  if (!file.exists(res_path)) return(invisible())
  res <- readRDS(res_path)
  if (is.null(res$clusters) || is.null(res$trb) || is.null(res$enrich)) {
    return(invisible())
  }

  top_groups <- res$enrich |>
    dplyr::arrange(FDR) |>
    dplyr::slice_head(n = 30) |>
    dplyr::pull(cluster_id)

  clust_long <- res$clusters |>
    dplyr::filter(cluster_id %in% top_groups) |>
    dplyr::left_join(res$trb[, c("CDR3b", "Subject", "Phenotype_2")],
                     by = "CDR3b") |>
    dplyr::filter(!is.na(Subject)) |>
    dplyr::count(cluster_id, Subject, Phenotype_2, name = "n_cdr3") |>
    dplyr::mutate(value = log1p(n_cdr3))

  if (nrow(clust_long) == 0) return(invisible())

  subj_order <- clust_long |>
    dplyr::distinct(Subject, Phenotype_2) |>
    dplyr::arrange(factor(Phenotype_2, levels = c("NIU", "Viral")), Subject) |>
    dplyr::pull(Subject)
  clust_long$Subject <- factor(clust_long$Subject, levels = subj_order)

  clust_order <- res$enrich |>
    dplyr::filter(cluster_id %in% top_groups) |>
    dplyr::arrange(median_OR) |>
    dplyr::pull(cluster_id)
  clust_long$cluster_id <- factor(clust_long$cluster_id, levels = clust_order)

  anno_df <- clust_long |>
    dplyr::distinct(Subject, Phenotype_2)

  p_anno <- ggplot(anno_df, aes(x = Subject, y = 1, fill = Phenotype_2)) +
    geom_tile() +
    scale_fill_manual(values = ETIOLOGY_GROUP_COLORS, name = NULL) +
    theme_void() +
    theme(legend.position = "top",
          legend.text  = element_text(size = 9),
          plot.margin  = margin(0, 0, 0, 0)) +
    guides(fill = guide_legend(override.aes = list(size = 4)))

  p_main <- ggplot(clust_long,
                   aes(x = Subject, y = cluster_id, fill = value)) +
    geom_tile() +
    viridis::scale_fill_viridis(option = "viridis",
                                name = "log1p(n CDR3)") +
    labs(title = "GLIPH group membership per subject",
         subtitle = sprintf("Top %d groups by FDR; subjects ordered NIU then Viral",
                            length(top_groups)),
         x = "Subject", y = "GLIPH group") +
    .tcr_compare_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- patchwork::wrap_plots(p_anno, p_main, ncol = 1,
                               heights = c(0.05, 1))
  } else p <- p_main

  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_gliph_subject_heatmap"),
               w = max(8, length(subj_order) * 0.35 + 3),
               h = max(7, length(top_groups) * 0.22 + 3))
}

# ---- Panel 5: GLIPH group network of top hits ----------------------------
# Network of CDR3 sequences within the top FDR groups. Edges connect CDR3s
# within the same group; nodes colored by Phenotype_2 of the subject(s).
.tcell_gliph_network <- function(paths, cfg) {
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("tidygraph", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE)) {
    log_message("  TCR-compare: ggraph/tidygraph not available; skipping network.")
    return(invisible())
  }
  res_path <- file.path(get_target_paths(cfg, "all")$results_objects,
                        "ImmGLIPHResults.rds")
  if (!file.exists(res_path)) return(invisible())
  res <- readRDS(res_path)
  if (is.null(res$clusters) || is.null(res$trb) || is.null(res$enrich)) {
    return(invisible())
  }

  top_groups <- res$enrich |>
    dplyr::filter(FDR < 0.1) |>
    dplyr::arrange(FDR) |>
    dplyr::slice_head(n = 12) |>
    dplyr::pull(cluster_id)
  if (length(top_groups) == 0) return(invisible())

  clust <- res$clusters |>
    dplyr::filter(cluster_id %in% top_groups) |>
    dplyr::distinct(cluster_id, CDR3b)

  trb_lookup <- res$trb |>
    dplyr::group_by(CDR3b) |>
    dplyr::summarise(Phenotype_2 = paste(sort(unique(Phenotype_2)),
                                         collapse = "/"),
                     n_cells = dplyr::n(), .groups = "drop")

  nodes <- clust |>
    dplyr::distinct(CDR3b) |>
    dplyr::left_join(trb_lookup, by = "CDR3b") |>
    dplyr::mutate(Phenotype_2 = dplyr::coalesce(Phenotype_2, "n.s."))

  # Build edges: connect all CDR3 pairs within the same group.
  edges <- clust |>
    dplyr::group_by(cluster_id) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::reframe(from = CDR3b[utils::combn(dplyr::n(), 2)[1, ]],
                   to   = CDR3b[utils::combn(dplyr::n(), 2)[2, ]]) |>
    dplyr::ungroup()

  if (nrow(edges) == 0) return(invisible())

  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

  pal_full <- c(ETIOLOGY_GROUP_COLORS,
                `NIU/Viral` = "grey50", `n.s.` = "grey80")

  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.18, color = "grey60", width = 0.3) +
    ggraph::geom_node_point(aes(color = Phenotype_2,
                                size  = n_cells)) +
    scale_color_manual(values = pal_full, name = "Phenotype") +
    scale_size_continuous(trans = "log10", range = c(1, 5),
                          name = "Cells") +
    labs(title = "GLIPH convergence-group network",
         subtitle = sprintf("Top %d groups by FDR; node = CDR3-beta, edge = same group",
                            length(top_groups))) +
    .tcr_compare_theme() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.line  = element_blank(),
          axis.ticks = element_blank())
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_gliph_network"),
               w = 10, h = 8)
}

# ---- Panel 6: tcrdist UMAP -----------------------------------------------
.tcell_tcrdist_umap <- function(paths, cfg) {
  res_path <- file.path(get_target_paths(cfg, "all")$results_objects,
                        "ImmLynxTcrdistResults.rds")
  if (!file.exists(res_path)) return(invisible())
  res <- readRDS(res_path)
  if (is.null(res$umap) || is.null(res$trb)) return(invisible())

  umap <- as.data.frame(res$umap)
  colnames(umap)[1:2] <- c("UMAP1", "UMAP2")
  umap$barcode     <- rownames(umap)
  trb              <- res$trb
  join_cols <- c("barcode", "Subject", "Etiology", "Phenotype_2", "clone_id")
  if ("HLA_B27_pathogenic" %in% names(trb))
    join_cols <- c(join_cols, "HLA_B27_pathogenic")
  df <- dplyr::left_join(umap, trb[, join_cols], by = "barcode")
  df <- df |>
    dplyr::group_by(clone_id) |>
    dplyr::mutate(clone_size = dplyr::n()) |>
    dplyr::ungroup()

  b27 <- if ("HLA_B27_pathogenic" %in% names(df))
           df[which(df$HLA_B27_pathogenic), , drop = FALSE]
         else df[FALSE, , drop = FALSE]

  p_pheno <- ggplot(df, aes(UMAP1, UMAP2,
                            color = Phenotype_2, size = clone_size)) +
    geom_point(alpha = 0.55) +
    geom_point(data = b27,
               aes(UMAP1, UMAP2, size = clone_size),
               inherit.aes = FALSE,
               shape = 21, color = "black", fill = NA, stroke = 0.7) +
    scale_color_manual(values = ETIOLOGY_GROUP_COLORS, name = "Phenotype") +
    scale_size_continuous(trans = "log10", range = c(0.4, 4),
                          name = "Clone size") +
    labs(title = "tcrdist + ESM-2 UMAP of intraocular TRB",
         subtitle = sprintf("Color: Viral vs NIU; size: log10 clone size; black ring = HLA-B27 pathogenic (n=%d)",
                            nrow(b27)),
         x = "UMAP-1", y = "UMAP-2") +
    .tcr_compare_theme()
  save_pdf_png(p_pheno, file.path(viz_subdir(paths, "tcr_motif"), "tcell_tcrdist_umap_phenotype"),
               w = 8, h = 7)

  # Per-etiology facet: color by subtype palette so Viral subfamilies (blue
  # hues) and NIU subfamilies (red hues) remain distinguishable.
  df$Etiology <- factor(df$Etiology,
                        levels = intersect(names(ETIOLOGY_SUBTYPE_COLORS),
                                           unique(df$Etiology)))
  p_et <- ggplot(df, aes(UMAP1, UMAP2,
                         color = Etiology, size = clone_size)) +
    geom_point(alpha = 0.6) +
    .etiology_subtype_scale(levels(df$Etiology), "color") +
    scale_size_continuous(trans = "log10", range = c(0.3, 3)) +
    facet_wrap(~ Etiology, ncol = 4) +
    labs(title = "tcrdist UMAP by individual etiology",
         x = "UMAP-1", y = "UMAP-2") +
    .tcr_compare_theme() +
    theme(legend.position = "bottom")
  save_pdf_png(p_et,
               file.path(viz_subdir(paths, "tcr_motif"), "tcell_tcrdist_umap_by_etiology"),
               w = 12, h = 8)
}

# ---- Panel 7: KNN neighborhood density per subject -----------------------
.tcell_tcrdist_neighborhood <- function(paths, cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "tcrdist_per_subject.csv")
  if (!file.exists(csv)) return(invisible())
  ps <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(ps) == 0) return(invisible())

  glm_csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                       "repertoire", "tcrdist_glm_results.csv")
  subtitle <- "Per-subject mean of (same-Phenotype_2 neighbors / K)"
  if (file.exists(glm_csv)) {
    g <- utils::read.csv(glm_csv, stringsAsFactors = FALSE)
    pe <- g[grepl("Phenotype_2", g$term), , drop = FALSE]
    if (nrow(pe) == 1)
      subtitle <- sprintf("GLM Phenotype_2 effect: beta = %.3f, p = %.3g (cell-count adjusted)",
                          pe$estimate, pe$p)
  }

  ps$Etiology <- factor(ps$Etiology,
                        levels = intersect(names(ETIOLOGY_SUBTYPE_COLORS),
                                           unique(ps$Etiology)))

  p <- ggplot(ps, aes(x = Phenotype_2, y = mean_same_group_frac,
                      fill = Phenotype_2)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
    geom_jitter(width = 0.12, alpha = 0.95, size = 2.2,
                aes(color = Etiology)) +
    scale_fill_manual(values = ETIOLOGY_GROUP_COLORS, guide = "none") +
    .etiology_subtype_scale(levels(ps$Etiology), "color") +
    labs(title = "tcrdist KNN: within-phenotype neighborhood density",
         subtitle = subtitle,
         x = "Phenotype", y = "Mean same-group neighbor fraction") +
    .tcr_compare_theme()
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_tcrdist_neighborhood_density"),
               w = 8, h = 6)
}

# ---- Panel 8: OLGA log10 Pgen ridges -------------------------------------
.tcell_olga_pgen_ridges <- function(paths, cfg) {
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    log_message("  TCR-compare: ggridges missing; skipping Pgen ridges.")
    return(invisible())
  }
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "olga_pgen_per_clone.csv")
  if (!file.exists(csv)) return(invisible())
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())

  # OLGA returns NA or Pgen = 0 (-> log10 = -Inf) for sequences it can't
  # score. Drop those so density() can compute a bandwidth.
  n_in <- nrow(d)
  d <- d[is.finite(d$log10_Pgen), , drop = FALSE]
  if (nrow(d) < n_in)
    log_message("  Pgen ridges: dropped ", n_in - nrow(d),
                " rows with non-finite log10 Pgen.")
  if (nrow(d) < 10) {
    log_message("  Pgen ridges: too few finite rows; skipping.")
    return(invisible())
  }

  glm_csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                       "repertoire", "olga_glm_results.csv")
  subtitle <- "log10 OLGA generation probability per CDR3-beta"
  if (file.exists(glm_csv)) {
    g <- utils::read.csv(glm_csv, stringsAsFactors = FALSE)
    pe <- g[grepl("Phenotype_2", g$term), , drop = FALSE]
    if (nrow(pe) == 1)
      subtitle <- sprintf("GLM Phenotype_2 effect on per-subject median: beta = %.3f, p = %.3g (cell-count adjusted)",
                          pe$estimate, pe$p)
  }

  p <- ggplot(d, aes(x = log10_Pgen, y = Phenotype_2, fill = Phenotype_2)) +
    ggridges::geom_density_ridges(alpha = 0.75, scale = 0.95,
                                  quantile_lines = TRUE, quantiles = 0.5,
                                  na.rm = TRUE) +
    scale_fill_manual(values = ETIOLOGY_GROUP_COLORS, guide = "none") +
    labs(title = "OLGA generation probability: intraocular Viral vs NIU",
         subtitle = subtitle,
         x = "log10 Pgen", y = NULL) +
    .tcr_compare_theme()
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_olga_pgen_ridges"),
               w = 8, h = 4.5)

  d$Etiology <- factor(d$Etiology,
                       levels = intersect(names(ETIOLOGY_SUBTYPE_COLORS),
                                          unique(d$Etiology)))
  # Drop etiology levels that have <10 finite Pgen values -- density()
  # also blows up on too-thin distributions.
  et_n <- table(d$Etiology)
  keep_et <- names(et_n)[et_n >= 10]
  d_et <- d[as.character(d$Etiology) %in% keep_et, , drop = FALSE]
  d_et$Etiology <- droplevels(factor(d_et$Etiology, levels = keep_et))
  if (nrow(d_et) >= 10 && nlevels(d_et$Etiology) >= 1) {
    p_et <- ggplot(d_et, aes(x = log10_Pgen, y = Etiology, fill = Etiology)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 0.95,
                                    quantile_lines = TRUE, quantiles = 0.5,
                                    na.rm = TRUE) +
      .etiology_subtype_scale(levels(d_et$Etiology), "fill") +
      labs(title = "OLGA Pgen by individual etiology",
           x = "log10 Pgen", y = NULL) +
      .tcr_compare_theme() +
      theme(legend.position = "none")
    save_pdf_png(p_et,
                 file.path(viz_subdir(paths, "tcr_motif"), "tcell_olga_pgen_ridges_by_etiology"),
                 w = 8, h = 7)
  } else {
    log_message("  Pgen ridges by etiology: insufficient data; skipping.")
  }
}

# ---- Panel 9: Pgen x clone-size scatter ----------------------------------
# Mechanistic: antigen-driven public clones should be in the high-abundance,
# low-Pgen quadrant. NIU clones should distribute closer to background.
.tcell_pgen_vs_clone_size <- function(paths, cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "olga_pgen_per_clone.csv")
  if (!file.exists(csv)) return(invisible())
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())

  # Drop non-finite Pgen (NA or Pgen == 0 -> log10 = -Inf) so the
  # scatter axis is well-defined.
  d <- d[is.finite(d$log10_Pgen), , drop = FALSE]
  if (nrow(d) < 10) return(invisible())

  df <- d |>
    dplyr::group_by(CDR3b, Phenotype_2) |>
    dplyr::summarise(n_cells     = dplyr::n(),
                     log10_Pgen  = stats::median(log10_Pgen),
                     n_subjects  = dplyr::n_distinct(Subject),
                     .groups = "drop") |>
    dplyr::filter(is.finite(log10_Pgen)) |>
    dplyr::mutate(sharing = ifelse(n_subjects > 1, "public", "private"))

  # Label clones with >10 cells, capped per phenotype.
  label_df <- df |>
    dplyr::filter(n_cells > 10) |>
    dplyr::group_by(Phenotype_2) |>
    dplyr::slice_max(n_cells, n = 10, with_ties = FALSE) |>
    dplyr::ungroup()
  log_message("  Pgen-vs-size: labeling ", nrow(label_df),
              " clones (n_cells > 10, top 10 per phenotype).")

  size_vals <- c(public = 4.8, private = 1.2)

  p <- ggplot(df, aes(x = log10_Pgen, y = n_cells, color = Phenotype_2)) +
    geom_point(aes(size = sharing), alpha = 0.55) +
    ggrepel::geom_text_repel(data = label_df,
                             aes(label = CDR3b),
                             size = 2.4, max.overlaps = Inf,
                             min.segment.length = 0.1,
                             box.padding = 0.35,
                             force = 2,
                             show.legend = FALSE) +
    scale_y_log10() +
    scale_color_manual(values = ETIOLOGY_GROUP_COLORS, name = "Phenotype") +
    scale_size_manual(values = size_vals, name = "Sharing") +
    labs(title = "Generation probability vs clonal expansion",
         subtitle = "Antigen-driven candidates: low Pgen + high abundance",
         x = "log10 Pgen", y = "Cells per CDR3 (log10)") +
    .tcr_compare_theme() +
    theme(aspect.ratio = 1)
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_pgen_vs_clone_size"),
               w = 8, h = 8)
}

# Panel 9b -- Same scatter, colored by GLIPH directional assignment.
# Viral-enriched GLIPH = blue, NIU-enriched GLIPH = red, no significant
# GLIPH group = grey.
.tcell_pgen_vs_clone_size_gliph <- function(paths, cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "olga_pgen_per_clone.csv")
  if (!file.exists(csv)) return(invisible())
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())
  d <- d[is.finite(d$log10_Pgen), , drop = FALSE]
  if (nrow(d) < 10) return(invisible())

  lk <- .tcr_compare_per_cell_lookup(cfg)
  v_cdr3 <- lk$gliph_cdr3_viral %||% character(0)
  n_cdr3 <- lk$gliph_cdr3_niu   %||% character(0)
  if (length(v_cdr3) + length(n_cdr3) == 0) {
    log_message("  Pgen-vs-size (GLIPH): no enriched GLIPH CDR3s; skipping.")
    return(invisible())
  }

  df <- d |>
    dplyr::group_by(CDR3b) |>
    dplyr::summarise(n_cells     = dplyr::n(),
                     log10_Pgen  = stats::median(log10_Pgen),
                     n_subjects  = dplyr::n_distinct(Subject),
                     .groups = "drop") |>
    dplyr::filter(is.finite(log10_Pgen)) |>
    dplyr::mutate(
      sharing = ifelse(n_subjects > 1, "public", "private"),
      gliph_direction = dplyr::case_when(
        CDR3b %in% v_cdr3 & CDR3b %in% n_cdr3 ~ "Both (ambiguous)",
        CDR3b %in% v_cdr3                     ~ "Viral GLIPH",
        CDR3b %in% n_cdr3                     ~ "NIU GLIPH",
        TRUE                                  ~ "No GLIPH"
      )
    )

  df$gliph_direction <- factor(df$gliph_direction,
    levels = c("No GLIPH", "Viral GLIPH", "NIU GLIPH", "Both (ambiguous)"))

  # Plot order: No GLIPH (grey) first/bottom, then enriched on top.
  df <- df[order(df$gliph_direction), , drop = FALSE]

  pal <- c(`No GLIPH`         = "grey70",
           `Viral GLIPH`      = unname(ETIOLOGY_GROUP_COLORS["Viral"]),
           `NIU GLIPH`        = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
           `Both (ambiguous)` = "#6A0DAD")

  size_vals <- c(public = 4.8, private = 1.2)

  # Label only GLIPH-enriched clones with >10 cells, capped per direction.
  label_df <- df |>
    dplyr::filter(gliph_direction != "No GLIPH", n_cells > 10) |>
    dplyr::group_by(gliph_direction) |>
    dplyr::slice_max(n_cells, n = 10, with_ties = FALSE) |>
    dplyr::ungroup()
  log_message("  Pgen-vs-size (GLIPH): labeling ", nrow(label_df),
              " enriched clones (n_cells > 10, top 10 per direction).")

  p <- ggplot(df, aes(x = log10_Pgen, y = n_cells,
                      color = gliph_direction)) +
    geom_point(aes(size = sharing), alpha = 0.55) +
    ggrepel::geom_text_repel(data = label_df,
                             aes(label = CDR3b),
                             size = 2.4, max.overlaps = Inf,
                             min.segment.length = 0.1,
                             box.padding = 0.35,
                             force = 2,
                             show.legend = FALSE) +
    scale_y_log10() +
    scale_color_manual(values = pal, name = "GLIPH assignment") +
    scale_size_manual(values = size_vals, name = "Sharing") +
    labs(title = "Generation probability vs clonal expansion",
         subtitle = "Colored by GLIPH directional enrichment",
         x = "log10 Pgen", y = "Cells per CDR3 (log10)") +
    .tcr_compare_theme() +
    theme(aspect.ratio = 1)
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"),
                            "tcell_pgen_vs_clone_size_gliph"),
               w = 8, h = 8)
}

# ---- Panel 10: SoNNia selection volcano ----------------------------------
.tcell_sonnia_volcano <- function(paths, cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "sonnia_selection_features.csv")
  if (!file.exists(csv)) return(invisible())
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())

  q_col <- intersect(c("Q", "selection_factor", "q"), colnames(d))[1]
  p_col <- intersect(c("p", "p_value", "pvalue"), colnames(d))[1]
  if (is.na(q_col) || is.na(p_col)) {
    log_message("  TCR-compare: SoNNia table missing Q/p columns; skipping volcano.")
    return(invisible())
  }
  d$log2_Q <- log2(pmax(d[[q_col]], .Machine$double.xmin))
  d$nlogp  <- -log10(pmax(d[[p_col]], .Machine$double.xmin))
  d$direction <- dplyr::case_when(
    d[[p_col]] < 0.05 & d$log2_Q >  0 ~ "Viral_up",
    d[[p_col]] < 0.05 & d$log2_Q <  0 ~ "NIU_up",
    TRUE                              ~ "ns"
  )

  vfam <- if ("V" %in% colnames(d)) d$V
          else if ("TRBV" %in% colnames(d)) d$TRBV
          else NA
  d$V_family <- stringr::str_extract(as.character(vfam), "TRBV[0-9]+")
  feature_label <- d$feature %||%
                   ifelse(!is.na(d$V_family), d$V_family, as.character(seq_len(nrow(d))))

  top_each_side <- d |>
    dplyr::filter(direction != "ns") |>
    dplyr::group_by(direction) |>
    dplyr::slice_max(abs(log2_Q), n = 8, with_ties = FALSE) |>
    dplyr::ungroup()

  p <- ggplot(d, aes(x = log2_Q, y = nlogp, color = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(alpha = 0.75, size = 1.7) +
    ggrepel::geom_text_repel(data = top_each_side,
                             aes(label = feature),
                             size = 2.5, max.overlaps = Inf,
                             show.legend = FALSE) +
    scale_color_manual(values = ETIOLOGY_DIRECTION_COLORS,
                       labels = c(NIU_up = "NIU-favored",
                                  Viral_up = "Viral-favored",
                                  ns = "n.s."),
                       name = NULL) +
    labs(title = "SoNNia selection factors on TRB features",
         subtitle = "log2 Q > 0: enriched relative to OLGA background",
         x = "log2 Q (selection)", y = "-log10 p") +
    .tcr_compare_theme() +
    theme(legend.position = "top")
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_sonnia_selection_volcano"),
               w = 8.5, h = 6.5)
}

# ---- Panel 11: Public vs private clone proportions -----------------------
# Public = CDR3 present in >= 2 subjects within the same Phenotype_2 group.
# Reports both per-clone (counts) and per-cell (abundance-weighted) views.
.tcell_public_private_clones <- function(paths, cfg) {
  pgen_csv  <- file.path(get_target_paths(cfg, "all")$results_tables,
                         "repertoire", "olga_pgen_per_clone.csv")
  if (!file.exists(pgen_csv)) return(invisible())
  d <- utils::read.csv(pgen_csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())

  per_clone <- d |>
    dplyr::group_by(Phenotype_2, CDR3b) |>
    dplyr::summarise(n_subjects = dplyr::n_distinct(Subject),
                     n_cells    = dplyr::n(),
                     .groups = "drop") |>
    dplyr::mutate(sharing = ifelse(n_subjects >= 2, "Public", "Private"))

  by_clone <- per_clone |>
    dplyr::count(Phenotype_2, sharing, name = "n_clones") |>
    dplyr::group_by(Phenotype_2) |>
    dplyr::mutate(frac = n_clones / sum(n_clones)) |>
    dplyr::ungroup() |>
    dplyr::mutate(view = "Per-clone")

  by_cell <- per_clone |>
    dplyr::group_by(Phenotype_2, sharing) |>
    dplyr::summarise(n = sum(n_cells), .groups = "drop") |>
    dplyr::group_by(Phenotype_2) |>
    dplyr::mutate(frac = n / sum(n), n_clones = n) |>
    dplyr::ungroup() |>
    dplyr::select(Phenotype_2, sharing, n_clones, frac) |>
    dplyr::mutate(view = "Per-cell")

  combined <- dplyr::bind_rows(by_clone, by_cell)
  combined$view <- factor(combined$view, levels = c("Per-clone", "Per-cell"))
  combined$sharing <- factor(combined$sharing, levels = c("Private", "Public"))

  p <- ggplot(combined, aes(x = Phenotype_2, y = frac, fill = sharing)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.4) +
    geom_text(aes(label = scales::percent(frac, accuracy = 1)),
              position = position_stack(vjust = 0.5),
              color = "white", size = 3, fontface = "bold") +
    scale_fill_manual(values = c(Private = "grey75",
                                 Public  = unname(ETIOLOGY_GROUP_COLORS["NIU"])),
                      name = NULL) +
    facet_wrap(~ view) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Public vs private CDR3-beta",
         subtitle = "Public: CDR3 found in 2+ subjects of the same phenotype",
         x = "Phenotype", y = "Fraction") +
    .tcr_compare_theme()
  save_pdf_png(p, file.path(viz_subdir(paths, "repertoire"), "tcell_public_private_clones"),
               w = 7, h = 5)
}

# ---- Panel 12: Lorenz / cumulative top-N clone curve ---------------------
# Cumulative fraction of cells covered by the top N CDR3s, per Phenotype_2.
# Steeper rise = more oligoclonal repertoire. Per-subject curves shown as
# thin lines, group medians as thick lines. Aspect ratio is locked to 1:1
# (both axes are proportions on [0,1]) and the legend reports the per-
# subject Gini coefficient + a Wilcoxon test of Gini Viral vs NIU.
.tcell_clonality_curve <- function(paths, cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "olga_pgen_per_clone.csv")
  if (!file.exists(csv)) return(invisible())
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())

  per_subj <- d |>
    dplyr::count(Subject, Phenotype_2, CDR3b, name = "n_cells") |>
    dplyr::group_by(Subject, Phenotype_2) |>
    dplyr::arrange(dplyr::desc(n_cells), .by_group = TRUE) |>
    dplyr::mutate(rank      = dplyr::row_number(),
                  total     = sum(n_cells),
                  cum_frac  = cumsum(n_cells) / total,
                  rank_frac = rank / dplyr::n()) |>
    dplyr::ungroup()

  group_summary <- per_subj |>
    dplyr::mutate(rank_bin = round(rank_frac * 100) / 100) |>
    dplyr::group_by(Phenotype_2, rank_bin) |>
    dplyr::summarise(median_cum = stats::median(cum_frac), .groups = "drop")

  # Per-subject Gini coefficient on clone-count vector.
  .gini <- function(x) {
    x <- sort(x[!is.na(x) & x > 0])
    n <- length(x)
    if (n < 2) return(NA_real_)
    (2 * sum(seq_len(n) * x) - (n + 1) * sum(x)) / (n * sum(x))
  }
  gini_by_subj <- per_subj |>
    dplyr::group_by(Subject, Phenotype_2) |>
    dplyr::summarise(gini = .gini(n_cells),
                     n_clones = dplyr::n(),
                     .groups = "drop") |>
    dplyr::filter(!is.na(gini))
  utils::write.csv(gini_by_subj,
                   file.path(get_target_paths(cfg, "all")$results_tables,
                             "repertoire", "clonality_gini_per_subject.csv"),
                   row.names = FALSE)

  # Wilcoxon Viral vs NIU on Gini.
  stat_line <- if (sum(gini_by_subj$Phenotype_2 == "Viral") >= 3 &&
                   sum(gini_by_subj$Phenotype_2 == "NIU")   >= 3) {
    w <- tryCatch(
      stats::wilcox.test(gini ~ Phenotype_2, data = gini_by_subj),
      error = function(e) NULL
    )
    med <- gini_by_subj |>
      dplyr::group_by(Phenotype_2) |>
      dplyr::summarise(med = stats::median(gini), .groups = "drop")
    medV <- med$med[med$Phenotype_2 == "Viral"]
    medN <- med$med[med$Phenotype_2 == "NIU"]
    if (!is.null(w))
      sprintf("Gini median: Viral=%.2f, NIU=%.2f; Wilcoxon p = %.3g",
              medV, medN, w$p.value)
    else
      sprintf("Gini median: Viral=%.2f, NIU=%.2f", medV, medN)
  } else "Gini test n.a. (insufficient subjects per group)"

  p <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey70", linewidth = 0.4) +
    geom_line(data = per_subj,
              aes(x = rank_frac, y = cum_frac,
                  group = Subject, color = Phenotype_2),
              alpha = 0.35, linewidth = 0.4) +
    geom_line(data = group_summary,
              aes(x = rank_bin, y = median_cum, color = Phenotype_2),
              linewidth = 1.4) +
    scale_color_manual(values = ETIOLOGY_GROUP_COLORS, name = "Phenotype") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1), expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    labs(title = "Cumulative clonality curve (Lorenz)",
         subtitle = paste0(
           "Steeper = more oligoclonal. Thin = per subject, thick = group median\n",
           stat_line),
         x = "Top fraction of unique CDR3s (ranked by abundance)",
         y = "Cumulative fraction of cells") +
    .tcr_compare_theme()
  save_pdf_png(p, file.path(viz_subdir(paths, "repertoire"), "tcell_clonality_curve"),
               w = 7, h = 7)

  # Companion Gini boxplot -- discrete view of the same statistic.
  if (nrow(gini_by_subj) >= 3) {
    # Attach per-subject Etiology so jitter can be colored by subtype.
    if ("Etiology" %in% colnames(d)) {
      gini_by_subj <- gini_by_subj |>
        dplyr::left_join(dplyr::distinct(d[, c("Subject", "Etiology")]),
                         by = "Subject")
      gini_by_subj$Etiology <- factor(gini_by_subj$Etiology,
        levels = intersect(names(ETIOLOGY_SUBTYPE_COLORS),
                           unique(gini_by_subj$Etiology)))
      p_gini <- ggplot(gini_by_subj,
                       aes(x = Phenotype_2, y = gini, fill = Phenotype_2)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
        geom_jitter(width = 0.12, alpha = 0.95, size = 2.2,
                    aes(color = Etiology)) +
        scale_fill_manual(values = ETIOLOGY_GROUP_COLORS, guide = "none") +
        .etiology_subtype_scale(levels(gini_by_subj$Etiology), "color")
    } else {
      p_gini <- ggplot(gini_by_subj,
                       aes(x = Phenotype_2, y = gini, fill = Phenotype_2)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
        geom_jitter(width = 0.12, alpha = 0.95, size = 2.2) +
        scale_fill_manual(values = ETIOLOGY_GROUP_COLORS, guide = "none")
    }
    p_gini <- p_gini +
      labs(title = "Per-subject Gini coefficient (clone-size inequality)",
           subtitle = stat_line,
           x = "Phenotype", y = "Gini") +
      .tcr_compare_theme()
    save_pdf_png(p_gini,
                 file.path(viz_subdir(paths, "repertoire"), "tcell_clonality_gini_box"),
                 w = 7, h = 6)
  }
}

# ---- Panel 13: UMAP overlays driven by TCR analysis ----------------------
# These run on the T cell compartment Seurat object (main gene-expression
# UMAP). They paint TCR-derived labels onto the canonical embedding so
# clinical / mechanistic features can be read off the same coordinates the
# rest of the manuscript uses.

# Internal helper: load lookup tables (b27 cdr3s, gliph membership by
# direction, per-cell Pgen) and return a tibble keyed by barcode.
.tcr_compare_per_cell_lookup <- function(cfg) {
  tables <- file.path(get_target_paths(cfg, "all")$results_tables, "repertoire")
  res <- list()

  b27_csv <- file.path(tables, "hla_b27_pathogenic_clones.csv")
  if (file.exists(b27_csv)) {
    b27 <- utils::read.csv(b27_csv, stringsAsFactors = FALSE)
    res$b27_barcodes <- b27$barcode
    res$b27_cdr3     <- unique(b27$CDR3b)
  }

  gliph_rds <- file.path(get_target_paths(cfg, "all")$results_objects,
                         "ImmGLIPHResults.rds")
  if (file.exists(gliph_rds)) {
    g <- readRDS(gliph_rds)
    if (!is.null(g$enrich) && !is.null(g$clusters)) {
      sig <- g$enrich |>
        dplyr::filter(!is.na(FDR), FDR < 0.1) |>
        dplyr::mutate(direction = ifelse(median_OR > 1, "Viral", "NIU"))
      vir_ids <- sig$cluster_id[sig$direction == "Viral"]
      niu_ids <- sig$cluster_id[sig$direction == "NIU"]
      res$gliph_cdr3_viral <- unique(g$clusters$CDR3b[g$clusters$cluster_id
                                                      %in% vir_ids])
      res$gliph_cdr3_niu   <- unique(g$clusters$CDR3b[g$clusters$cluster_id
                                                      %in% niu_ids])
    }
  }

  pgen_csv <- file.path(tables, "olga_pgen_per_clone.csv")
  if (file.exists(pgen_csv)) {
    p <- utils::read.csv(pgen_csv, stringsAsFactors = FALSE)
    res$pgen_by_bc <- setNames(p$log10_Pgen, p$barcode)
  }

  res
}

# Internal helper: stamp a "tcr_label" column on the Seurat object so a
# given highlight set is plotable via CellDimPlot(group_by = "tcr_label").
.stamp_highlight_col <- function(obj, col_name, highlight_mask,
                                 label_yes, label_no = "Other") {
  vec <- rep(label_no, ncol(obj))
  vec[highlight_mask] <- label_yes
  obj[[col_name]] <- factor(vec, levels = c(label_no, label_yes))
  obj
}

# Internal helper: extract per-cell TRB CDR3 from CTaa for matching.
.cell_trb_cdr3 <- function(obj) {
  ctaa <- obj@meta.data$CTaa
  trb  <- stringr::str_split(ctaa, "_", simplify = TRUE)
  if (ncol(trb) >= 2) trb[, 2] else rep(NA_character_, ncol(obj))
}

# Panel 13a -- HLA-B27 pathogenic clones on the main UMAP.
# Uses scplotter highlight = '<expr>' so target cells draw on top with a
# stroke and full opacity, like the Mudd-lab TFH UMAP spike highlight.
.tcell_umap_hla_b27 <- function(obj, paths, cfg) {
  lk <- .tcr_compare_per_cell_lookup(cfg)
  if (is.null(lk$b27_cdr3) || length(lk$b27_cdr3) == 0) {
    log_message("  UMAP B27: no HLA-B27 pathogenic CDR3s available; skipping.")
    return(invisible())
  }
  trb_aa <- .cell_trb_cdr3(obj)
  mask   <- !is.na(trb_aa) & trb_aa %in% lk$b27_cdr3
  if (sum(mask) == 0) {
    log_message("  UMAP B27: no cells matched pathogenic CDR3s; skipping.")
    return(invisible())
  }
  log_message("  UMAP B27: highlighting ", sum(mask),
              " cells (", length(unique(trb_aa[mask])), " unique CDR3-beta).")
  obj <- .stamp_highlight_col(obj, "tcr_hla_b27",
                              mask, label_yes = "HLA-B27 pathogenic")

  pal <- c(Other                = "grey85",
           `HLA-B27 pathogenic` = "#9D0208")
  umap_sz <- cfg$visualization$umap_size %||% 8
  base    <- file.path(viz_subdir(paths, "tcr_motif"), "tcell_umap_hla_b27_pathogenic")
  dual_save_plot(bquote(
    scplotter::CellDimPlot(obj, group_by = "tcr_hla_b27",
                           reduction = "UMAP",
                           highlight = 'tcr_hla_b27 == "HLA-B27 pathogenic"',
                           palcolor = .(pal),
                           bg_color = "grey92",
                           pt_alpha = 0.5,
                           highlight_size = 2.4,
                           highlight_stroke = 0.6,
                           raster = FALSE,
                           show_stat = FALSE) +
      ggplot2::ggtitle(.(sprintf("HLA-B27 pathogenic clones (n=%d cells)",
                                 sum(mask))))
  ), base, width = umap_sz, height = umap_sz)
}

# Panel 13b -- GLIPH Viral-enriched and NIU-enriched CDR3s on main UMAP.
# Same highlight idiom; both directions get full-opacity points colored by
# Phenotype_2 family, background cells fade to light grey.
.tcell_umap_gliph_clusters <- function(obj, paths, cfg) {
  lk <- .tcr_compare_per_cell_lookup(cfg)
  if (is.null(lk$gliph_cdr3_viral) && is.null(lk$gliph_cdr3_niu)) {
    log_message("  UMAP GLIPH: no enriched cluster CDR3s available; skipping.")
    return(invisible())
  }
  trb_aa <- .cell_trb_cdr3(obj)
  in_v <- !is.na(trb_aa) & trb_aa %in% (lk$gliph_cdr3_viral %||% character(0))
  in_n <- !is.na(trb_aa) & trb_aa %in% (lk$gliph_cdr3_niu   %||% character(0))
  cls <- rep("Other", ncol(obj))
  cls[in_v]        <- "Viral GLIPH-enriched"
  cls[in_n]        <- "NIU GLIPH-enriched"
  cls[in_v & in_n] <- "Both (ambiguous)"
  if (sum(cls != "Other") == 0) {
    log_message("  UMAP GLIPH: no cells in enriched clusters; skipping.")
    return(invisible())
  }
  log_message("  UMAP GLIPH: Viral=", sum(cls == "Viral GLIPH-enriched"),
              ", NIU=", sum(cls == "NIU GLIPH-enriched"),
              ", ambiguous=", sum(cls == "Both (ambiguous)"))
  obj$tcr_gliph_dir <- factor(cls,
    levels = c("Other", "Viral GLIPH-enriched",
               "NIU GLIPH-enriched", "Both (ambiguous)"))

  pal <- c(Other                  = "grey85",
           `Viral GLIPH-enriched` = unname(ETIOLOGY_GROUP_COLORS["Viral"]),
           `NIU GLIPH-enriched`   = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
           `Both (ambiguous)`     = "#6A0DAD")
  umap_sz <- cfg$visualization$umap_size %||% 8
  base    <- file.path(viz_subdir(paths, "tcr_motif"), "tcell_umap_gliph_directional")
  dual_save_plot(bquote(
    scplotter::CellDimPlot(obj, group_by = "tcr_gliph_dir",
                           reduction = "UMAP",
                           highlight = 'tcr_gliph_dir != "Other"',
                           palcolor = .(pal),
                           bg_color = "grey92",
                           pt_alpha = 0.5,
                           highlight_size = 1.8,
                           highlight_stroke = 0.4,
                           raster = FALSE,
                           show_stat = FALSE) +
      ggplot2::ggtitle("GLIPH directional-enriched CDR3-beta")
  ), base, width = umap_sz, height = umap_sz)
}

# Panel 13c-bar -- HLA-B27 pathogenic enrichment bar plots.
# Plots the per-group proportion of cells flagged as HLA-B27 pathogenic
# across (i) T cell substate, (ii) pooled Phenotype_2, (iii) individual
# Etiology subtype. Implemented via FeatureStatPlot(plot_type="bar",
# agg=mean) on a 0/1 numeric column so the bar height is directly the
# fraction of cells in that group carrying the signature.
.tcell_bar_hla_b27 <- function(obj, paths, cfg) {
  lk <- .tcr_compare_per_cell_lookup(cfg)
  if (is.null(lk$b27_cdr3) || length(lk$b27_cdr3) == 0) {
    log_message("  Bar B27: no HLA-B27 pathogenic CDR3s available; skipping.")
    return(invisible())
  }
  trb_aa <- .cell_trb_cdr3(obj)
  obj$is_hla_b27 <- as.numeric(!is.na(trb_aa) & trb_aa %in% lk$b27_cdr3)
  if (sum(obj$is_hla_b27) == 0) {
    log_message("  Bar B27: no matching cells; skipping.")
    return(invisible())
  }
  if ("knn.leiden.cluster" %in% colnames(obj@meta.data))
    obj$substate_label <- substate_labels(cfg, "tcell",
                                          obj$knn.leiden.cluster)

  pct_y <- ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1)
  )

  # By substate
  if ("substate_label" %in% colnames(obj@meta.data)) {
    log_message("  FeatureStatPlot bar HLA-B27 by substate ...")
    p <- tryCatch(
      scplotter::FeatureStatPlot(obj, features = "is_hla_b27",
                                 plot_type = "bar",
                                 ident = "substate_label",
                                 agg = mean,
                                 x_text_angle = 45,
                                 show_stat = FALSE) +
        pct_y +
        ggplot2::ylab("Fraction HLA-B27 pathogenic") +
        ggplot2::ggtitle("HLA-B27 pathogenic frequency per T cell substate"),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p))
      save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"),
                                "tcell_bar_hla_b27_by_substate"),
                   w = max(8, length(unique(obj$substate_label)) * 0.7 + 2),
                   h = 5)
  }

  # By Phenotype_2 (pooled etiology family)
  log_message("  FeatureStatPlot bar HLA-B27 by Phenotype_2 ...")
  p <- tryCatch(
    scplotter::FeatureStatPlot(obj, features = "is_hla_b27",
                               plot_type = "bar",
                               ident = "Phenotype_2",
                               agg = mean,
                               palcolor = ETIOLOGY_GROUP_COLORS,
                               show_stat = FALSE) +
      pct_y +
      ggplot2::ylab("Fraction HLA-B27 pathogenic") +
      ggplot2::ggtitle("HLA-B27 pathogenic frequency per phenotype"),
    error = function(e) {
      log_message("    failed: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(p))
    save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"),
                              "tcell_bar_hla_b27_by_phenotype"),
                 w = 6, h = 5)

  # By individual Etiology, ordered Viral subtypes -> NIU subtypes so the
  # phenotype family is visually grouped.
  if ("Etiology" %in% colnames(obj@meta.data)) {
    et_levels <- intersect(names(ETIOLOGY_SUBTYPE_COLORS),
                           unique(obj$Etiology))
    obj$Etiology <- factor(obj$Etiology, levels = et_levels)
    log_message("  FeatureStatPlot bar HLA-B27 by Etiology ...")
    p <- tryCatch(
      scplotter::FeatureStatPlot(obj, features = "is_hla_b27",
                                 plot_type = "bar",
                                 ident = "Etiology",
                                 agg = mean,
                                 palcolor = ETIOLOGY_SUBTYPE_COLORS,
                                 x_text_angle = 45,
                                 show_stat = FALSE) +
        pct_y +
        ggplot2::ylab("Fraction HLA-B27 pathogenic") +
        ggplot2::ggtitle("HLA-B27 pathogenic frequency per individual etiology"),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p))
      save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"),
                                "tcell_bar_hla_b27_by_etiology"),
                   w = max(7, length(et_levels) * 0.6 + 2), h = 5)
  }
}

# Panel 13c-gliph -- GLIPH directional bar plots.
# Two binary features (is_gliph_viral, is_gliph_niu) plotted side-by-side
# via FeatureStatPlot's multi-feature mode. Same three stratifications as
# the HLA-B27 panel.
.tcell_bar_gliph_directional <- function(obj, paths, cfg) {
  lk <- .tcr_compare_per_cell_lookup(cfg)
  if (is.null(lk$gliph_cdr3_viral) && is.null(lk$gliph_cdr3_niu)) {
    log_message("  Bar GLIPH: no enriched CDR3s available; skipping.")
    return(invisible())
  }
  trb_aa <- .cell_trb_cdr3(obj)
  obj$is_gliph_viral <- as.numeric(
    !is.na(trb_aa) & trb_aa %in% (lk$gliph_cdr3_viral %||% character(0))
  )
  obj$is_gliph_niu   <- as.numeric(
    !is.na(trb_aa) & trb_aa %in% (lk$gliph_cdr3_niu   %||% character(0))
  )
  if (sum(obj$is_gliph_viral) + sum(obj$is_gliph_niu) == 0) {
    log_message("  Bar GLIPH: no matching cells; skipping.")
    return(invisible())
  }
  if ("knn.leiden.cluster" %in% colnames(obj@meta.data))
    obj$substate_label <- substate_labels(cfg, "tcell",
                                          obj$knn.leiden.cluster)

  feats <- c("is_gliph_viral", "is_gliph_niu")
  pct_y <- ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1)
  )

  if ("substate_label" %in% colnames(obj@meta.data)) {
    log_message("  FeatureStatPlot bar GLIPH-directional by substate ...")
    p <- tryCatch(
      scplotter::FeatureStatPlot(obj, features = feats,
                                 plot_type = "bar",
                                 ident = "substate_label",
                                 agg = mean,
                                 x_text_angle = 45,
                                 show_stat = FALSE) +
        pct_y +
        ggplot2::ylab("Fraction in GLIPH cluster") +
        ggplot2::ggtitle("GLIPH directional CDR3-beta per T cell substate"),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p))
      save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"),
                                "tcell_bar_gliph_directional_by_substate"),
                   w = max(10, length(unique(obj$substate_label)) * 0.7 + 3),
                   h = 6)
  }

  log_message("  FeatureStatPlot bar GLIPH-directional by Phenotype_2 ...")
  p <- tryCatch(
    scplotter::FeatureStatPlot(obj, features = feats,
                               plot_type = "bar",
                               ident = "Phenotype_2",
                               agg = mean,
                               palcolor = ETIOLOGY_GROUP_COLORS,
                               show_stat = FALSE) +
      pct_y +
      ggplot2::ylab("Fraction in GLIPH cluster") +
      ggplot2::ggtitle("GLIPH directional CDR3-beta per phenotype"),
    error = function(e) {
      log_message("    failed: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(p))
    save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"),
                              "tcell_bar_gliph_directional_by_phenotype"),
                 w = 8, h = 5)

  if ("Etiology" %in% colnames(obj@meta.data)) {
    et_levels <- intersect(names(ETIOLOGY_SUBTYPE_COLORS),
                           unique(obj$Etiology))
    obj$Etiology <- factor(obj$Etiology, levels = et_levels)
    log_message("  FeatureStatPlot bar GLIPH-directional by Etiology ...")
    p <- tryCatch(
      scplotter::FeatureStatPlot(obj, features = feats,
                                 plot_type = "bar",
                                 ident = "Etiology",
                                 agg = mean,
                                 palcolor = ETIOLOGY_SUBTYPE_COLORS,
                                 x_text_angle = 45,
                                 show_stat = FALSE) +
        pct_y +
        ggplot2::ylab("Fraction in GLIPH cluster") +
        ggplot2::ggtitle("GLIPH directional CDR3-beta per individual etiology"),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p))
      save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"),
                                "tcell_bar_gliph_directional_by_etiology"),
                   w = max(10, length(et_levels) * 0.7 + 3), h = 6)
  }
}

# Panel 13d -- scRepertoire cloneSize on the main UMAP.
# cloneSize is an ordered factor (Single -> ... -> Hyperexpanded). Viridis
# default direction maps the last level (Hyperexpanded) to yellow so the
# most expanded clones stand out.
.tcell_umap_clone_size <- function(obj, paths, cfg) {
  meta <- obj@meta.data
  if (!"cloneSize" %in% colnames(meta)) {
    log_message("  UMAP cloneSize: cloneSize column missing in metadata; skipping.")
    return(invisible())
  }
  cs <- meta$cloneSize
  if (all(is.na(cs))) {
    log_message("  UMAP cloneSize: all NA; skipping.")
    return(invisible())
  }
  # Ensure factor ordering runs Single -> Hyperexpanded so viridis yellow
  # lands on Hyperexpanded.
  if (!is.factor(cs)) cs <- factor(cs)
  lvls <- levels(cs)
  rank_key <- c("Rare" = 0, "Single" = 1, "Small" = 2, "Medium" = 3,
                "Large" = 4, "Hyperexpanded" = 5)
  prefix <- sub("[^A-Za-z].*$", "", lvls)
  ord <- order(unname(rank_key[prefix] %||% seq_along(lvls)))
  if (!any(is.na(rank_key[prefix])) && length(unique(rank_key[prefix])) > 1)
    cs <- factor(cs, levels = lvls[ord], ordered = TRUE)
  obj$cloneSize <- cs

  log_message("  UMAP cloneSize: ", sum(!is.na(cs)),
              " cells with assignment; levels = ",
              paste(levels(cs), collapse = " | "))

  umap_sz <- cfg$visualization$umap_size %||% 8
  base    <- file.path(viz_subdir(paths, "tcr_motif"), "tcell_umap_clone_size")
  dual_save_plot(quote(
    scplotter::CellDimPlot(obj, group_by = "cloneSize",
                           reduction = "UMAP",
                           palette = "viridis",
                           bg_color = "grey92",
                           raster = FALSE,
                           show_stat = FALSE) +
      ggplot2::ggtitle("scRepertoire clone size (yellow = hyperexpanded)")
  ), base, width = umap_sz, height = umap_sz)
}

# Panel 13c -- log10 Pgen on UMAP, overall + split by Phenotype_2
.tcell_umap_pgen <- function(obj, paths, cfg) {
  lk <- .tcr_compare_per_cell_lookup(cfg)
  if (is.null(lk$pgen_by_bc) || length(lk$pgen_by_bc) == 0) {
    log_message("  UMAP Pgen: olga_pgen_per_clone.csv not available; skipping.")
    return(invisible())
  }
  pgen <- lk$pgen_by_bc[colnames(obj)]
  if (sum(is.finite(pgen)) < 50) {
    log_message("  UMAP Pgen: <50 cells with finite Pgen; skipping.")
    return(invisible())
  }
  obj$log10_Pgen <- as.numeric(pgen)

  umap_sz <- cfg$visualization$umap_size %||% 8
  stripped_bg <- cfg$visualization$stripped_bg %||% "black"

  base_all <- file.path(viz_subdir(paths, "tcr_motif"), "tcell_umap_pgen")
  dual_save_plot(bquote(
    scplotter::FeatureStatPlot(obj, features = "log10_Pgen",
                               plot_type = "dim",
                               reduction = "UMAP",
                               palette = "viridis",
                               raster = FALSE,
                               show_stat = FALSE,
                               bg_cutoff = -Inf) +
      ggplot2::ggtitle("log10 OLGA Pgen (TRB)")
  ), base_all, width = umap_sz, height = umap_sz,
     stripped_bg = stripped_bg)

  base_split <- file.path(viz_subdir(paths, "tcr_motif"), "tcell_umap_pgen_split_phenotype")
  dual_save_plot(bquote(
    scplotter::FeatureStatPlot(obj, features = "log10_Pgen",
                               plot_type = "dim",
                               reduction = "UMAP",
                               split_by = "Phenotype_2",
                               palette = "viridis",
                               raster = FALSE,
                               show_stat = FALSE,
                               bg_cutoff = -Inf) +
      ggplot2::ggtitle("log10 OLGA Pgen split by Phenotype_2")
  ), base_split, width = umap_sz * 1.8, height = umap_sz,
     stripped_bg = stripped_bg)
}

# ---------------------------------------------------------------------------
# scRepertoire descriptive panels
#
# The T cell compartment Seurat object already carries scRepertoire metadata
# (CTaa, CTgene, CTstrict, cloneSize) from combineExpression() earlier in the
# pipeline, so the visualization functions accept obj directly.
# Filename stems: tcell_screp_<function>_<grouping>.{pdf,png}.
# ---------------------------------------------------------------------------

# Internal: subset the T cell object to cells with valid TCR data and a
# resolved substate label. Returns NULL if nothing usable.
.tcr_tcell_subset_for_screp <- function(obj, cfg, drop_healthy = TRUE) {
  if (!requireNamespace("scRepertoire", quietly = TRUE)) return(NULL)
  meta <- obj@meta.data
  has_tcr <- !is.na(meta$CTstrict %||% meta$CTaa)
  in_groups <- if (drop_healthy)
                 meta$Phenotype_2 %in% c("Viral", "NIU")
               else
                 !is.na(meta$Phenotype_2)
  keep <- has_tcr & in_groups
  if (sum(keep) < 50) return(NULL)
  o <- subset(obj, cells = colnames(obj)[keep])
  if ("knn.leiden.cluster" %in% colnames(o@meta.data)) {
    o$substate_label <- substate_labels(cfg, "tcell",
                                        o$knn.leiden.cluster)
  }
  o
}

.tcell_screp_homeostasis <- function(obj, paths, cfg) {
  o <- .tcr_tcell_subset_for_screp(obj, cfg)
  if (is.null(o)) {
    log_message("  scRepertoire: skipping homeostasis (no TCR-carrying T cells)."); return(invisible())
  }
  cc <- cfg$repertoire$clone_definition %||% "strict"

  for (grp in c("Phenotype_2", "Etiology", "substate_label")) {
    if (!(grp %in% colnames(o@meta.data))) next
    log_message("  scRepertoire clonalHomeostasis by ", grp, " ...")
    p <- tryCatch(
      scRepertoire::clonalHomeostasis(o,
                                      cloneCall = cc,
                                      chain = "TRB",
                                      group.by = grp,
                                      palette = "viridis"),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p)) {
      n_grp <- length(unique(o@meta.data[[grp]]))
      save_pdf_png(p,
                   file.path(viz_subdir(paths, "repertoire"),
                             paste0("tcell_screp_homeostasis_", grp)),
                   w = max(6, n_grp * 0.6 + 3), h = 5)
    }
  }
}

.tcell_screp_proportion <- function(obj, paths, cfg) {
  o <- .tcr_tcell_subset_for_screp(obj, cfg)
  if (is.null(o)) {
    log_message("  scRepertoire: skipping proportion."); return(invisible())
  }
  cc <- cfg$repertoire$clone_definition %||% "strict"

  for (grp in c("Phenotype_2", "Etiology", "substate_label")) {
    if (!(grp %in% colnames(o@meta.data))) next
    log_message("  scRepertoire clonalProportion by ", grp, " ...")
    p <- tryCatch(
      scRepertoire::clonalProportion(o,
                                     cloneCall = cc,
                                     chain = "TRB",
                                     group.by = grp,
                                     palette = "viridis"),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p)) {
      n_grp <- length(unique(o@meta.data[[grp]]))
      save_pdf_png(p,
                   file.path(viz_subdir(paths, "repertoire"),
                             paste0("tcell_screp_proportion_", grp)),
                   w = max(6, n_grp * 0.6 + 3), h = 5)
    }
  }
}

.tcell_screp_overlap <- function(obj, paths, cfg) {
  o <- .tcr_tcell_subset_for_screp(obj, cfg)
  if (is.null(o)) {
    log_message("  scRepertoire: skipping overlap."); return(invisible())
  }
  cc     <- cfg$repertoire$clone_definition %||% "strict"
  method <- cfg$repertoire$overlap_method   %||% "morisita"

  # 1) Overall: Etiology x Etiology
  log_message("  scRepertoire clonalOverlap (Etiology x Etiology, ",
              method, ") ...")
  p <- tryCatch(
    scRepertoire::clonalOverlap(o,
                                cloneCall = cc,
                                chain = "TRB",
                                method = method,
                                group.by = "Etiology"),
    error = function(e) {
      log_message("    failed: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(p))
    save_pdf_png(p, file.path(viz_subdir(paths, "repertoire"),
                              "tcell_screp_overlap_etiology"),
                 w = 8, h = 7)

  # 2 + 3) Subject x Subject within Viral and within NIU
  for (ph in c("Viral", "NIU")) {
    sub_cells <- colnames(o)[o$Phenotype_2 == ph]
    if (length(sub_cells) < 50) {
      log_message("  scRepertoire overlap (", ph, "): <50 cells; skipping.")
      next
    }
    s <- subset(o, cells = sub_cells)
    n_sub <- length(unique(s$Subject))
    if (n_sub < 2) {
      log_message("  scRepertoire overlap (", ph, "): only ", n_sub,
                  " subject; skipping.")
      next
    }
    log_message("  scRepertoire clonalOverlap (Subject x Subject within ",
                ph, ", ", method, ") ...")
    p <- tryCatch(
      scRepertoire::clonalOverlap(s,
                                  cloneCall = cc,
                                  chain = "TRB",
                                  method = method,
                                  group.by = "Subject"),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p))
      save_pdf_png(p, file.path(viz_subdir(paths, "repertoire"),
                                paste0("tcell_screp_overlap_subjects_",
                                       tolower(ph))),
                   w = max(7, n_sub * 0.35 + 3),
                   h = max(6, n_sub * 0.35 + 3))
  }
}

.tcell_screp_startrac <- function(obj, paths, cfg) {
  o <- .tcr_tcell_subset_for_screp(obj, cfg)
  if (is.null(o)) {
    log_message("  scRepertoire: skipping StartracDiversity."); return(invisible())
  }
  if (!"Subject" %in% colnames(o@meta.data)) {
    log_message("  StartracDiversity: Subject column missing; skipping.")
    return(invisible())
  }
  cc <- cfg$repertoire$clone_definition %||% "strict"

  # type = patient/subject column (one row per subject in the output);
  # group.by = the stratifier (Phenotype_2 / Etiology) to color points.
  for (grp in c("Phenotype_2", "Etiology")) {
    if (!(grp %in% colnames(o@meta.data))) next
    log_message("  scRepertoire StartracDiversity (type=Subject, group.by=",
                grp, ") ...")
    p <- tryCatch(
      scRepertoire::StartracDiversity(o,
                                      cloneCall = cc,
                                      chain     = "TRB",
                                      type      = "Subject",
                                      group.by  = grp),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p))
      save_pdf_png(p,
                   file.path(viz_subdir(paths, "repertoire"),
                             paste0("tcell_screp_startrac_", grp)),
                   w = 9, h = 5)
  }
}

.tcell_screp_circlize <- function(obj, paths, cfg) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    log_message("  scRepertoire circlize: 'circlize' package missing; skipping.")
    return(invisible())
  }
  o <- .tcr_tcell_subset_for_screp(obj, cfg)
  if (is.null(o)) {
    log_message("  scRepertoire: skipping circlize."); return(invisible())
  }
  cc <- cfg$repertoire$clone_definition %||% "strict"

  for (grp in c("Phenotype_2", "Etiology", "substate_label")) {
    if (!(grp %in% colnames(o@meta.data))) next

    log_message("  scRepertoire getCirclize (group.by=", grp, ") ...")
    cmat <- tryCatch(
      scRepertoire::getCirclize(o,
                                cloneCall    = cc,
                                group.by     = grp,
                                proportion   = FALSE,
                                include.self = FALSE),
      error = function(e) {
        log_message("    failed: ", conditionMessage(e)); NULL
      }
    )
    if (is.null(cmat) || nrow(cmat) == 0) next

    # Sector color palette per grouping.
    sectors <- unique(c(as.character(cmat$from), as.character(cmat$to)))
    grid_cols <- if (grp == "Phenotype_2") {
      setNames(ifelse(sectors %in% names(ETIOLOGY_GROUP_COLORS),
                      ETIOLOGY_GROUP_COLORS[sectors], "grey50"),
               sectors)
    } else if (grp == "Etiology") {
      setNames(ifelse(sectors %in% names(ETIOLOGY_SUBTYPE_COLORS),
                      ETIOLOGY_SUBTYPE_COLORS[sectors], "grey50"),
               sectors)
    } else {
      # substate_label -- no global palette in the project, so generate a
      # viridis ramp keyed on natural sector ordering.
      setNames(viridis::viridis(length(sectors), option = "turbo"),
               sectors)
    }

    base <- file.path(viz_subdir(paths, "repertoire"),
                      paste0("tcell_screp_circlize_", grp))
    .draw_chord <- function() {
      circlize::circos.clear()
      circlize::circos.par(start.degree = 90, gap.degree = 4)
      circlize::chordDiagram(
        cmat,
        grid.col        = grid_cols,
        annotationTrack = c("name", "grid"),
        transparency    = 0.25
      )
      title(main = sprintf("Shared TCR clones across %s", grp),
            cex.main = 0.95)
      circlize::circos.clear()
    }
    # Substate panels have many short sector labels -- give them a slightly
    # larger canvas so the names don't collide.
    dim_in <- if (grp == "substate_label") 9 else 8
    ensure_dir(dirname(base))  # raw grDevices::pdf/png need the bucket to exist
    grDevices::pdf(paste0(base, ".pdf"), width = dim_in, height = dim_in)
    .draw_chord(); grDevices::dev.off()
    grDevices::png(paste0(base, ".png"),
                   width = dim_in, height = dim_in,
                   units = "in", res = 300)
    .draw_chord(); grDevices::dev.off()
    log_message("  Saved: tcell_screp_circlize_", grp, ".{pdf,png}")
  }
}

# ---- Panel 14: TRBV usage of HLA-B27 pathogenic clones -------------------
# Even though TRAV21 is required for the pathogenic label, the paired TRBV
# usage is informative -- are the paired beta chains drawn from a narrow
# subset? Compares HLA-B27 pathogenic vs the background intraocular pool.
.tcell_hla_b27_trbv_usage <- function(paths, cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "hla_b27_pathogenic_clones.csv")
  if (!file.exists(csv)) return(invisible())
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())

  # Background = all intraocular TRB rows (from GLIPH results if available,
  # otherwise from the b27 file's complement -- approximate).
  gliph_rds <- file.path(get_target_paths(cfg, "all")$results_objects,
                         "ImmGLIPHResults.rds")
  bg <- if (file.exists(gliph_rds)) {
          g <- readRDS(gliph_rds)
          if (!is.null(g$trb)) g$trb else NULL
        } else NULL
  if (is.null(bg)) {
    log_message("  HLA-B27 TRBV: background not available; skipping.")
    return(invisible())
  }

  bg <- bg[!is.na(bg$TRBV) & nzchar(bg$TRBV), , drop = FALSE]
  bg_freq <- bg |>
    dplyr::count(TRBV, name = "n") |>
    dplyr::mutate(set = "Background", frac = n / sum(n))
  b27_freq <- d |>
    dplyr::filter(!is.na(TRBV), nzchar(TRBV)) |>
    dplyr::count(TRBV, name = "n") |>
    dplyr::mutate(set = "HLA-B27 pathogenic", frac = n / sum(n))

  top_v <- union(b27_freq$TRBV,
                 bg_freq |>
                   dplyr::slice_max(n, n = 25, with_ties = FALSE) |>
                   dplyr::pull(TRBV))
  combo <- dplyr::bind_rows(bg_freq, b27_freq) |>
    dplyr::filter(TRBV %in% top_v)
  combo$TRBV <- factor(combo$TRBV, levels = top_v)

  p <- ggplot(combo, aes(x = TRBV, y = frac, fill = set)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    scale_fill_manual(values = c(Background = "grey60",
                                 `HLA-B27 pathogenic` = "#9D0208"),
                      name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Paired TRBV usage in HLA-B27 pathogenic clones",
         subtitle = "All have TRAV21 alpha-chain; this shows the paired beta",
         x = NULL, y = "Frequency within set") +
    .tcr_compare_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "top")
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_hla_b27_trbv_usage"),
               w = max(8, length(top_v) * 0.35 + 3), h = 5)
}

# ---- Panel 16: log10 Pgen by T cell substate -----------------------------
# Boxplot of per-cell log10 Pgen across the curated T cell substates,
# split by Phenotype_2. Tests the prediction that certain substates
# (effector, GZMK-high, exhausted) concentrate low-Pgen antigen-driven clones.
.tcell_pgen_by_substate <- function(obj, paths, cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "olga_pgen_per_clone.csv")
  if (!file.exists(csv)) return(invisible())
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  d <- d[is.finite(d$log10_Pgen), , drop = FALSE]
  if (nrow(d) < 20) return(invisible())

  meta <- obj@meta.data
  cb_col <- resolve_celltype_broad(meta) %||% "knn.leiden.cluster"
  if (!"knn.leiden.cluster" %in% colnames(meta)) return(invisible())

  meta$barcode <- rownames(meta)
  meta$substate_label <- substate_labels(cfg, "tcell",
                                         meta$knn.leiden.cluster)
  j <- dplyr::inner_join(d, meta[, c("barcode", "substate_label",
                                     "Phenotype_2")],
                         by = c("barcode", "Phenotype_2"))
  if (nrow(j) < 20) return(invisible())

  p <- ggplot(j, aes(x = substate_label, y = log10_Pgen, fill = Phenotype_2)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
    scale_fill_manual(values = ETIOLOGY_GROUP_COLORS, name = "Phenotype") +
    labs(title = "OLGA Pgen by T cell substate",
         subtitle = "Lower Pgen suggests antigen-driven selection",
         x = NULL, y = "log10 Pgen") +
    .tcr_compare_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          legend.position = "top")
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_pgen_by_substate"),
               w = max(8, length(unique(j$substate_label)) * 0.7 + 3),
               h = 5.5)
}

# ---- Panel 17: GLIPH x HLA-B27 contingency -------------------------------
# Are HLA-B27 pathogenic CDR3s preferentially captured by NIU-enriched
# GLIPH groups? Fisher OR + 95% CI; barplot of overlap counts.
.tcell_gliph_b27_contingency <- function(paths, cfg) {
  gliph_rds <- file.path(get_target_paths(cfg, "all")$results_objects,
                         "ImmGLIPHResults.rds")
  b27_csv   <- file.path(get_target_paths(cfg, "all")$results_tables,
                         "repertoire", "hla_b27_pathogenic_clones.csv")
  if (!file.exists(gliph_rds) || !file.exists(b27_csv)) return(invisible())
  g <- readRDS(gliph_rds)
  b27 <- utils::read.csv(b27_csv, stringsAsFactors = FALSE)
  if (is.null(g$clusters) || is.null(g$enrich) || nrow(b27) == 0)
    return(invisible())

  b27_cdr3 <- unique(b27$CDR3b)
  cls <- g$clusters |> dplyr::distinct(cluster_id, CDR3b)
  enr <- g$enrich |>
    dplyr::mutate(direction = ifelse(median_OR > 1, "Viral", "NIU"),
                  sig = !is.na(FDR) & FDR < 0.1)

  # Tag each unique CDR3 with the strongest-FDR cluster direction it belongs
  # to (Viral-enriched > NIU-enriched > unenriched > none).
  cls <- cls |>
    dplyr::left_join(enr[, c("cluster_id", "direction", "FDR", "sig")],
                     by = "cluster_id")
  cdr3_dir <- cls |>
    dplyr::group_by(CDR3b) |>
    dplyr::summarise(
      assignment = dplyr::case_when(
        any(sig & direction == "Viral", na.rm = TRUE) ~ "Viral GLIPH",
        any(sig & direction == "NIU",   na.rm = TRUE) ~ "NIU GLIPH",
        TRUE                                          ~ "Unenriched GLIPH"
      ), .groups = "drop")

  cdr3_dir$is_b27 <- cdr3_dir$CDR3b %in% b27_cdr3

  # Build a 2-column table even when one side is empty -- e.g. when none
  # of the B27 pathogenic CDR3s landed in any GLIPH cluster.
  cdr3_dir$b27_lab <- factor(ifelse(cdr3_dir$is_b27, "B27+", "B27-"),
                             levels = c("B27+", "B27-"))
  cdr3_dir$assignment <- factor(cdr3_dir$assignment,
    levels = c("Viral GLIPH", "NIU GLIPH", "Unenriched GLIPH"))
  tbl <- table(cdr3_dir$assignment, cdr3_dir$b27_lab)

  n_b27_total <- sum(tbl[, "B27+"])
  if (n_b27_total == 0) {
    log_message("  GLIPH x B27: 0 HLA-B27 pathogenic CDR3-beta found in any ",
                "GLIPH cluster (",
                length(b27_cdr3), " unique B27+ CDR3s, ",
                length(unique(cls$CDR3b)), " GLIPH-assigned CDR3s). ",
                "Skipping contingency plot.")
    # Still write the empty contingency table for the record.
    utils::write.csv(as.data.frame.matrix(tbl),
                     file.path(get_target_paths(cfg, "all")$results_tables,
                               "repertoire",
                               "gliph_x_hla_b27_contingency.csv"))
    return(invisible())
  }

  # Fisher OR for each direction vs not.
  rows <- lapply(rownames(tbl), function(r) {
    a <- tbl[r, "B27+"]; b <- tbl[r, "B27-"]
    c <- sum(tbl[, "B27+"]) - a; d <- sum(tbl[, "B27-"]) - b
    ft <- tryCatch(stats::fisher.test(matrix(c(a, c, b, d), nrow = 2)),
                   error = function(e) NULL)
    data.frame(assignment = r,
               n_B27_pos = a, n_B27_neg = b,
               OR    = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
               OR_lo = if (!is.null(ft)) ft$conf.int[1]      else NA_real_,
               OR_hi = if (!is.null(ft)) ft$conf.int[2]      else NA_real_,
               p     = if (!is.null(ft)) ft$p.value           else NA_real_)
  })
  ct <- do.call(rbind, rows)
  utils::write.csv(ct,
                   file.path(get_target_paths(cfg, "all")$results_tables,
                             "repertoire", "gliph_x_hla_b27_contingency.csv"),
                   row.names = FALSE)

  pal <- c(`Viral GLIPH`       = unname(ETIOLOGY_GROUP_COLORS["Viral"]),
           `NIU GLIPH`         = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
           `Unenriched GLIPH`  = "grey70")

  p <- ggplot(ct, aes(x = assignment, y = OR, color = assignment)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_pointrange(aes(ymin = OR_lo, ymax = OR_hi), size = 0.7) +
    geom_text(aes(label = sprintf("OR=%.1f, p=%.2g\n(B27+:%d, B27-:%d)",
                                  OR, p, n_B27_pos, n_B27_neg)),
              vjust = -0.7, size = 2.7, color = "grey20") +
    scale_color_manual(values = pal, guide = "none") +
    scale_y_log10() +
    labs(title = "GLIPH directional enrichment of HLA-B27 pathogenic CDR3s",
         subtitle = "OR > 1 means the assignment contains more B27+ CDR3-beta than expected",
         x = NULL, y = "Odds ratio (B27+ vs B27-)") +
    .tcr_compare_theme()
  save_pdf_png(p, file.path(viz_subdir(paths, "tcr_motif"), "tcell_gliph_x_hla_b27"),
               w = 7.5, h = 5.5)
}

# ---- Panel 14: HLA-B27 pathogenic clone summary --------------------------
# Cells flagged by flag_hla_b27_pathogenic_tcr() (TRAV21 + CDR3b [YF]S[TS]
# motif). Two views in one figure: per-subject count colored by Etiology,
# and per-Phenotype_2 stacked summary. Also writes a per-clone CSV.
.tcell_hla_b27_pathogenic_summary <- function(paths, cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "hla_b27_pathogenic_clones.csv")
  if (!file.exists(csv)) {
    log_message("  TCR-compare: hla_b27_pathogenic_clones.csv missing; skipping.")
    return(invisible())
  }
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())

  d$Etiology <- factor(d$Etiology,
                       levels = intersect(names(ETIOLOGY_SUBTYPE_COLORS),
                                          unique(d$Etiology)))

  per_subject <- d |>
    dplyr::group_by(Subject, Etiology, Phenotype_2) |>
    dplyr::summarise(n_cells       = dplyr::n(),
                     n_unique_cdr3 = dplyr::n_distinct(CDR3b),
                     .groups = "drop") |>
    dplyr::arrange(factor(Phenotype_2, levels = c("NIU", "Viral")),
                   dplyr::desc(n_cells))
  per_subject$Subject <- factor(per_subject$Subject,
                                levels = per_subject$Subject)

  p_subj <- ggplot(per_subject,
                   aes(x = Subject, y = n_cells, fill = Etiology)) +
    geom_col(width = 0.85, color = "white", linewidth = 0.2) +
    geom_text(aes(label = n_unique_cdr3), vjust = -0.3, size = 2.6) +
    .etiology_subtype_scale(levels(per_subject$Etiology), "fill") +
    facet_grid(. ~ Phenotype_2, scales = "free_x", space = "free_x") +
    labs(title = "HLA-B27 pathogenic clones per subject",
         subtitle = "Bar = cells with TRAV21 + CDR3b [YF]S[TS] motif; numerals = unique CDR3-beta",
         x = NULL, y = "Cells") +
    .tcr_compare_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "right")
  save_pdf_png(p_subj,
               file.path(viz_subdir(paths, "tcr_motif"), "tcell_hla_b27_pathogenic_per_subject"),
               w = max(8, nrow(per_subject) * 0.35 + 3), h = 6)

  # Etiology subtype stack within each Phenotype_2
  by_et <- d |>
    dplyr::count(Phenotype_2, Etiology, name = "n_cells")
  p_stack <- ggplot(by_et, aes(x = Phenotype_2, y = n_cells, fill = Etiology)) +
    geom_col(width = 0.55, color = "white", linewidth = 0.3) +
    .etiology_subtype_scale(levels(d$Etiology), "fill") +
    labs(title = "HLA-B27 pathogenic clone yield by phenotype",
         subtitle = "Stack: contribution from each individual etiology",
         x = "Phenotype", y = "Cells") +
    .tcr_compare_theme() +
    theme(legend.position = "right")
  save_pdf_png(p_stack,
               file.path(viz_subdir(paths, "tcr_motif"), "tcell_hla_b27_pathogenic_by_phenotype"),
               w = 7, h = 5)

  # Concise CDR3-level table for the methods/supplement
  per_clone <- d |>
    dplyr::group_by(CDR3b, TRBV, TRBJ) |>
    dplyr::summarise(n_cells      = dplyr::n(),
                     n_subjects   = dplyr::n_distinct(Subject),
                     n_etiologies = dplyr::n_distinct(Etiology),
                     subjects     = paste(sort(unique(Subject)),
                                          collapse = ";"),
                     etiologies   = paste(sort(unique(Etiology)),
                                          collapse = ";"),
                     phenotypes   = paste(sort(unique(Phenotype_2)),
                                          collapse = ";"),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(n_cells))
  utils::write.csv(per_clone,
                   file.path(get_target_paths(cfg, "all")$results_tables,
                             "repertoire", "hla_b27_pathogenic_per_cdr3.csv"),
                   row.names = FALSE)
  log_message("  Saved: hla_b27_pathogenic_per_cdr3.csv (",
              nrow(per_clone), " unique CDR3-beta)")
}

# ---- Panel 15: Cross-method top-50 summary table -------------------------
.tcell_cross_method_summary <- function(paths, cfg) {
  trb_tables <- file.path(get_target_paths(cfg, "all")$results_tables, "repertoire")
  out_csv <- file.path(trb_tables, "tcr_summary_top50.csv")

  pgen_csv  <- file.path(trb_tables, "olga_pgen_per_clone.csv")
  neigh_csv <- file.path(trb_tables, "tcrdist_neighborhoods.csv")
  cl_csv    <- file.path(trb_tables, "gliph_clusters.csv")

  if (!file.exists(pgen_csv) && !file.exists(neigh_csv) && !file.exists(cl_csv)) {
    return(invisible())
  }

  pg <- if (file.exists(pgen_csv))
          utils::read.csv(pgen_csv, stringsAsFactors = FALSE) else NULL
  nb <- if (file.exists(neigh_csv))
          utils::read.csv(neigh_csv, stringsAsFactors = FALSE) else NULL
  cl <- if (file.exists(cl_csv))
          utils::read.csv(cl_csv,    stringsAsFactors = FALSE) else NULL

  base <- pg %||% nb
  if (is.null(base)) return(invisible())
  has_pgen <- !is.null(pg) && "log10_Pgen" %in% names(pg)
  base <- base |>
    dplyr::group_by(CDR3b) |>
    dplyr::summarise(
      n_cells     = dplyr::n(),
      n_subjects  = dplyr::n_distinct(Subject),
      Phenotype_2 = paste(sort(unique(Phenotype_2)), collapse = "/"),
      log10_Pgen  = if (has_pgen) stats::median(log10_Pgen, na.rm = TRUE)
                    else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_cells)) |>
    dplyr::slice_head(n = 50)

  if (!is.null(nb)) {
    nb_summ <- nb |>
      dplyr::group_by(CDR3b) |>
      dplyr::summarise(mean_same_group_frac = mean(same_group_frac),
                       .groups = "drop")
    base <- dplyr::left_join(base, nb_summ, by = "CDR3b")
  }
  if (!is.null(cl)) {
    cl_summ <- cl |>
      dplyr::group_by(CDR3b) |>
      dplyr::summarise(gliph_group = paste(unique(cluster_id), collapse = ";"),
                       .groups = "drop")
    base <- dplyr::left_join(base, cl_summ, by = "CDR3b")
  }

  utils::write.csv(base, out_csv, row.names = FALSE)
  log_message("  Saved: tcr_summary_top50.csv (", nrow(base), " clones)")
}

# run_visualizations_tcell renders the T cell compartment figure block:
# cross-cutting cluster panels (substate UMAP labeled+stripped, UMAP by
# etiology, marker dot plot, NIU-vs-Viral composition, compartment-global
# and per-substate volcanoes, GSEA pathway-by-substate heatmap, NIU-vs-Viral
# top pathway bar, pathway-by-etiology heatmap, NIU sub-contrast heatmap)
# plus adaptive-receptor panels (CDR3-beta length NIU vs Viral, TRBV usage
# heatmap NIU vs Viral, public-clone count at multiple sharing thresholds),
# and T-cell-specific blocks: top-pathway UCell UMAP, eye-blood TCR overlap
# boxplot, HLA-B27-stratified TCR features (CDR3 length and TRBV usage split
# into NIU_B27+, NIU_other, Viral), and top-expanded-clone alluvial eye ->
# blood. The full MiloR viz block runs last so the tcell milo nhood-graph,
# beeswarm, boxplot, summary, composition, volcano, and nhood-size hist all
# land under outputs/viz/eye/tcell/07_milo/. All filenames are snake_case
# descriptions; no F5 panel-letter scheme is used.
