# R/85_viz_myeloid.R
# Myeloid compartment figure block. Uses the shared helpers in
# R/81_viz_compartment_helpers.R (substate UMAP, marker dot plot,
# NIU-vs-Viral composition / volcano / GSEA heatmap / pathway bar,
# per-substate volcano grid, pathway-by-etiology heatmap, NIU sub-contrast
# heatmap) and adds myeloid-specific panels:
#   - Lineage-marker UMAP featureplots
#   - Complement module score by substate and NIU vs Viral
#   - APC MHC-I/II module dotplot and Panel D / G cross-compartment views
# It also wires viz_milo (the full F1/F2 milo viz block) so myeloid gets the
# same beeswarm + boxplot + summary + composition + volcano + nhoodgraph +
# nhood-size hist as the eye level.
#
# All filenames use snake_case descriptions (no F3_panelA). Titles describe
# the panel content rather than referencing a figure letter.

suppressPackageStartupMessages({
  library(Seurat)
  library(scplotter)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# Top differential pathways from compartment GSEA, used as UCell overlay set
# ---------------------------------------------------------------------------
.top_axis_pathways_myeloid <- function(paths, n_per_axis = 2) {
  gsea_path <- file.path(paths$results_tables, "GSEA_Autoimmune_vs_Viral.csv")
  if (!file.exists(gsea_path)) return(character(0))
  g <- read.csv(gsea_path, stringsAsFactors = FALSE)
  if (!"axis_tag" %in% colnames(g)) return(character(0))
  g_curated <- g |>
    dplyr::filter(!is.na(axis_tag), stratum != "global", !is.na(p_val_adj))
  if (nrow(g_curated) == 0) return(character(0))
  g_curated |>
    dplyr::group_by(axis_tag) |>
    dplyr::slice_min(p_val_adj, n = n_per_axis, with_ties = FALSE) |>
    dplyr::pull(pathway) |>
    unique()
}

# ---------------------------------------------------------------------------
# Marker featureplot grid (lineage markers + key activation markers).
# ---------------------------------------------------------------------------
.myeloid_marker_featureplots <- function(obj, paths, cfg) {
  genes <- c("CD1A", "CD1C", "RORC", "IFI27", "CXCL10",
             "GPNMB", "LYVE1", "FCN1", "CLEC4C",
             "LAMP3", "CCR7", "IDO1", "S100A8", "S100A9")
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(invisible())
  umap_sz <- cfg$visualization$umap_size %||% 8
  stripped_bg <- cfg$visualization$stripped_bg %||% "black"
  dual_save_plot(bquote(
    FeatureStatPlot(obj, features = .(genes),
                    plot_type = "dim", reduction = "UMAP",
                    palette = "viridis", raster = FALSE, show_stat = FALSE)
  ), file.path(viz_subdir(paths, "markers"), "myeloid_lineage_marker_umap"),
    width = umap_sz * 1.2, height = umap_sz, stripped_bg = stripped_bg)
}

# ---------------------------------------------------------------------------
# Panel D variants: three sibling views of myeloid substate transcriptional
# response to Autoimmune (NIU) vs Viral. All three pull from data the existing
# DGE step already produced (one long-format pseudobulk CSV) plus the
# compartment Seurat object. The three views are saved side by side so the
# survivor can be picked on visual inspection before locking the figure.
#   v1: DEG counts per substate, split by direction (likely winner on clarity)
#   v2: 1 - Spearman r between AI and Viral mean pseudobulks per substate
#   v3: per-substate pseudobulk PCA small multiples with etiology ellipses
# Output: paths$viz_dir/11_pca_coupling/myeloid_etiology_distance_variants/{pdf,png}
# ---------------------------------------------------------------------------
.myeloid_panelD_paths <- function(paths) {
  out_dir <- file.path(viz_subdir(paths, "pca_coupling"), "myeloid_etiology_distance_variants")
  ensure_dir(out_dir)
  out_dir
}

# Curated gene categories used to color the panel D v3 loadings axis labels
# so the figure can serve as a launchpad for T cell analysis. HLA / antigen
# processing genes get one color, T-cell-interaction ligands and receptors
# another, everything else stays neutral. Extend either list with no other
# code changes; the labeller resolves matches via fixed string and prefix
# tests.
.PANEL_D_HLA_GENES <- c(
  # Class I
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
  "B2M", "TAP1", "TAP2", "TAPBP", "ERAP1", "ERAP2",
  "PSMB8", "PSMB9", "PSMB10", "NLRC5",
  # Class II + assembly
  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
  "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  "CD74", "CIITA", "IFI30",
  # Cathepsins involved in MHC-II processing
  "CTSS", "CTSB", "CTSL", "CTSH")

.PANEL_D_TCELL_INTERACTION_GENES <- c(
  # Cytokines
  "IL12A", "IL12B", "IL15", "IL18", "IL23A", "IL27", "IL10",
  "IL6", "IL1B", "TNF", "TGFB1", "IFNB1", "IFNG",
  # Chemokines that recruit / shape T cell behavior
  "CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL2", "CCL3", "CCL4",
  "CCL5", "CCL17", "CCL19", "CCL21", "CCL22",
  # Costim and checkpoint
  "CD80", "CD86", "CD40", "ICOSLG", "CD70",
  "TNFSF4", "TNFSF8", "TNFSF9", "TNFSF10", "TNFSF13B", "TNFSF14",
  "CD274", "PDCD1LG2", "VSIR", "LGALS9", "HVEM", "CD200",
  "IDO1", "IDO2", "CD58", "CADM1", "SLAMF8")

.gene_category <- function(genes) {
  cat <- rep("other", length(genes))
  cat[genes %in% .PANEL_D_HLA_GENES] <- "HLA"
  cat[genes %in% .PANEL_D_TCELL_INTERACTION_GENES] <- "Tcell"
  cat
}

# Fine-grained categories used by the curated panels (v4 heatmap, v5
# lollipop, v6 dumbbell). Splits the coarse HLA / Tcell labels into the
# six functional axes so the curated panels can be grouped by biology
# rather than by "in-list vs not".
.PANEL_D_CATEGORY_HLA_I <- c(
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
  "B2M", "TAP1", "TAP2", "TAPBP", "ERAP1", "ERAP2",
  "PSMB8", "PSMB9", "PSMB10", "NLRC5")
.PANEL_D_CATEGORY_HLA_II <- c(
  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
  "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  "CD74", "CIITA", "IFI30",
  "CTSS", "CTSB", "CTSL", "CTSH")
.PANEL_D_CATEGORY_CYTOKINE <- c(
  "IL12A", "IL12B", "IL15", "IL18", "IL23A", "IL27", "IL10",
  "IL6", "IL1B", "TNF", "TGFB1", "IFNB1", "IFNG")
.PANEL_D_CATEGORY_CHEMOKINE <- c(
  "CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL2", "CCL3", "CCL4",
  "CCL5", "CCL17", "CCL19", "CCL21", "CCL22")
.PANEL_D_CATEGORY_COSTIM <- c(
  "CD80", "CD86", "CD40", "ICOSLG", "CD70",
  "TNFSF4", "TNFSF8", "TNFSF9", "TNFSF10", "TNFSF13B", "TNFSF14",
  "CD58", "CADM1", "SLAMF8")
.PANEL_D_CATEGORY_CHECKPOINT <- c(
  "CD274", "PDCD1LG2", "VSIR", "LGALS9", "HVEM", "CD200",
  "IDO1", "IDO2")
.PANEL_D_CATEGORY_LEVELS <- c(
  "HLA-I",
  "HLA-II",
  "Cytokine",
  "Chemokine",
  "Costim",
  "Checkpoint")

.gene_category_fine <- function(genes) {
  out <- rep(NA_character_, length(genes))
  out[genes %in% .PANEL_D_CATEGORY_HLA_I]      <- .PANEL_D_CATEGORY_LEVELS[1]
  out[genes %in% .PANEL_D_CATEGORY_HLA_II]     <- .PANEL_D_CATEGORY_LEVELS[2]
  out[genes %in% .PANEL_D_CATEGORY_CYTOKINE]   <- .PANEL_D_CATEGORY_LEVELS[3]
  out[genes %in% .PANEL_D_CATEGORY_CHEMOKINE]  <- .PANEL_D_CATEGORY_LEVELS[4]
  out[genes %in% .PANEL_D_CATEGORY_COSTIM]     <- .PANEL_D_CATEGORY_LEVELS[5]
  out[genes %in% .PANEL_D_CATEGORY_CHECKPOINT] <- .PANEL_D_CATEGORY_LEVELS[6]
  out
}

# Long-format dataframe of curated HLA / T cell loadings across the focus
# substates, annotated with within-substate |loading| rank against the full
# variable-gene set. One row per (substate, PC, gene). Used by v4/v5/v6.
.myeloid_panelD_curated_long <- function(per_sub_focus, cfg) {
  rows <- lapply(per_sub_focus, function(x) {
    df <- x$loadings
    df$rank_PC1 <- rank(-abs(df$PC1), ties.method = "min")
    df$rank_PC2 <- rank(-abs(df$PC2), ties.method = "min")
    cat <- .gene_category_fine(df$gene)
    keep <- !is.na(cat)
    df  <- df[keep, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    sub_lab <- substate_labels(cfg, "myeloid", x$substate)
    rbind(
      data.frame(substate_label = sub_lab,
                 gene = df$gene,
                 category = factor(cat[keep],
                                   levels = .PANEL_D_CATEGORY_LEVELS),
                 pc = "PC1",
                 loading = df$PC1,
                 rank = df$rank_PC1,
                 stringsAsFactors = FALSE),
      data.frame(substate_label = sub_lab,
                 gene = df$gene,
                 category = factor(cat[keep],
                                   levels = .PANEL_D_CATEGORY_LEVELS),
                 pc = "PC2",
                 loading = df$PC2,
                 rank = df$rank_PC2,
                 stringsAsFactors = FALSE))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# Build the markdown-formatted axis label for a gene_key vector. Color +
# bold for HLA (gold) and T cell interaction (teal); plain for everything
# else. The gene_key encodes "<substate>||<pc>||<gene>" so we strip the
# prefix before classifying.
.panelD_loading_axis_labels <- function(keys,
                                        hla_color = "#B9770E",
                                        tcell_color = "#117864",
                                        other_color = "grey25") {
  genes <- sub(".*\\|\\|", "", as.character(keys))
  cat <- .gene_category(genes)
  fmt <- ifelse(
    cat == "HLA",
    sprintf("<span style='color:%s;font-weight:bold'>%s</span>",
            hla_color, genes),
    ifelse(cat == "Tcell",
           sprintf("<span style='color:%s;font-weight:bold'>%s</span>",
                   tcell_color, genes),
           sprintf("<span style='color:%s'>%s</span>",
                   other_color, genes)))
  fmt
}

# ---------------------------------------------------------------------------
# v4 curated heatmap. Rows = curated HLA / T-cell-interaction genes (grouped
# by category), cols = focus substate, faceted by PC. Fill = signed loading
# on a NIU<->Viral diverging palette. Cell text = |loading| rank in that
# substate against the full ~2000 variable gene set (blank when rank > 50)
# so the calibration against the rest of the transcriptome stays visible.
# ---------------------------------------------------------------------------
.myeloid_panelD_v4_curated_heatmap <- function(per_sub_focus, paths, cfg,
                                                focus_labels) {
  long <- .myeloid_panelD_curated_long(per_sub_focus, cfg)
  if (is.null(long) || nrow(long) == 0) {
    log_message("  panel D v4 (curated heatmap): no curated genes; skipping.")
    return(invisible())
  }
  long$substate_label <- factor(long$substate_label, levels = focus_labels)
  long$pc <- factor(long$pc, levels = c("PC1", "PC2"))
  gene_order <- long |>
    dplyr::filter(.data$pc == "PC1") |>
    dplyr::group_by(.data$category, .data$gene) |>
    dplyr::summarise(mu = mean(.data$loading, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(.data$category, .data$mu) |>
    dplyr::pull(.data$gene)
  long$gene <- factor(long$gene, levels = unique(gene_order))
  long$rank_lbl <- ifelse(long$rank <= 50,
                          as.character(long$rank), "")
  lim <- max(abs(long$loading), na.rm = TRUE)

  p <- ggplot2::ggplot(long,
                       ggplot2::aes(.data$substate_label, .data$gene,
                                    fill = .data$loading)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$rank_lbl),
                       size = 2.3, color = "grey15") +
    ggplot2::scale_fill_gradient2(
      low      = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
      mid      = "white",
      high     = unname(ETIOLOGY_GROUP_COLORS["Viral"]),
      midpoint = 0, limits = c(-lim, lim),
      name     = "Signed loading\nneg = AI, pos = Viral") +
    ggh4x::facet_grid2(rows = ggplot2::vars(.data$category),
                       cols = ggplot2::vars(.data$pc),
                       scales = "free_y", space = "free_y",
                       strip = ggh4x::strip_themed(
                         text_y = ggh4x::elem_list_text(face = "bold",
                                                         size = 9))) +
    ggplot2::labs(
      title    = "Curated HLA / T cell PCA loadings by substate",
      subtitle = paste0("Cell text = |loading| rank within that substate's ",
                        "PCA (blank if > 50). Rows sorted by mean PC1 ",
                        "loading within category."),
      x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title       = ggplot2::element_text(face = "bold"),
                   plot.subtitle    = ggplot2::element_text(size = 9),
                   panel.grid       = ggplot2::element_blank(),
                   axis.text.x      = ggplot2::element_text(angle = 30,
                                                            hjust = 1),
                   axis.text.y      = ggplot2::element_text(size = 8),
                   strip.text.y     = ggplot2::element_text(angle = 0,
                                                            hjust = 0,
                                                            face = "bold"),
                   legend.position  = "right")

  n_genes <- nlevels(long$gene)
  n_subs  <- length(focus_labels)
  save_pdf_png(p, file.path(.myeloid_panelD_paths(paths),
                            "v4_pseudobulk_pca_curated_heatmap"),
               w = max(8, 1.0 * n_subs + 6),
               h = max(8, 0.20 * n_genes + 2.5))
}

# ---------------------------------------------------------------------------
# v5 curated lollipop. One file per PC. Rows = category, x = signed loading,
# y = gene (ordered by mean loading within category), color = substate. Each
# row of the facet stack is one functional axis of T cell engagement, so the
# reader can see which axes each substate uses to separate AI from Viral.
# ---------------------------------------------------------------------------
.myeloid_panelD_v5_curated_lollipop <- function(per_sub_focus, paths, cfg,
                                                 focus_labels, strip_pal) {
  long <- .myeloid_panelD_curated_long(per_sub_focus, cfg)
  if (is.null(long) || nrow(long) == 0) {
    log_message("  panel D v5 (curated lollipop): no curated genes; skipping.")
    return(invisible())
  }
  long$substate_label <- factor(long$substate_label, levels = focus_labels)
  sub_pal <- strip_pal[focus_labels]
  if (any(is.na(sub_pal))) {
    sub_pal <- grDevices::hcl.colors(length(focus_labels), palette = "Dark 3")
    names(sub_pal) <- focus_labels
  }

  for (this_pc in c("PC1", "PC2")) {
    df <- long[long$pc == this_pc, , drop = FALSE]
    if (nrow(df) == 0) next
    gene_order <- df |>
      dplyr::group_by(.data$category, .data$gene) |>
      dplyr::summarise(mu = mean(.data$loading, na.rm = TRUE),
                       .groups = "drop") |>
      dplyr::arrange(.data$category, .data$mu) |>
      dplyr::pull(.data$gene)
    df$gene <- factor(df$gene, levels = unique(gene_order))

    p <- ggplot2::ggplot(df,
                         ggplot2::aes(.data$loading, .data$gene,
                                      color = .data$substate_label)) +
      ggplot2::geom_vline(xintercept = 0,
                          linewidth = 0.2, color = "grey60") +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$loading,
                                         yend = .data$gene),
                            linewidth = 0.3, alpha = 0.6) +
      ggplot2::geom_point(size = 0.9) +
      ggh4x::facet_grid2(rows = ggplot2::vars(.data$category),
                         scales = "free_y", space = "free_y",
                         strip = ggh4x::strip_themed(
                           text_y = ggh4x::elem_list_text(face = "bold",
                                                           size = 4,
                                                           angle = 270))) +
      ggplot2::scale_color_manual(values = sub_pal, name = NULL) +
      ggplot2::labs(x = paste0(this_pc, " loading"), y = NULL) +
      ggplot2::theme_bw(base_size = 6) +
      ggplot2::theme(panel.grid.minor   = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     axis.text          = ggplot2::element_text(size = 5),
                     axis.title         = ggplot2::element_text(size = 6),
                     strip.text.y       = ggplot2::element_text(angle = 270,
                                                                hjust = 0.5,
                                                                vjust = 0.5,
                                                                face = "bold",
                                                                size = 4),
                     legend.position    = "bottom",
                     legend.key.size    = grid::unit(0.25, "lines"),
                     legend.text        = ggplot2::element_text(size = 5),
                     legend.margin      = ggplot2::margin(0, 0, 0, 0),
                     legend.box.margin  = ggplot2::margin(-4, 0, 0, 0),
                     plot.margin        = ggplot2::margin(2, 4, 2, 2))

    save_pdf_png(p, file.path(.myeloid_panelD_paths(paths),
                              paste0("v5_pseudobulk_curated_lollipop_",
                                     this_pc)),
                 w = 2, h = 2.67)
  }
}

# ---------------------------------------------------------------------------
# v6 curated dumbbell. PC1 only. y = curated gene, x = signed PC1 loading;
# one point per focus substate connected by a segment from min to max. Red
# segment = the gene flips sign across substates (i.e. substate-specific
# direction along the AI<->Viral axis), grey = consistent direction. Designed
# to surface which T cell / HLA genes are substate-specific vs pan-myeloid.
# ---------------------------------------------------------------------------
.myeloid_panelD_v6_curated_dumbbell <- function(per_sub_focus, paths, cfg,
                                                 focus_labels, strip_pal) {
  long <- .myeloid_panelD_curated_long(per_sub_focus, cfg)
  if (is.null(long) || nrow(long) == 0) {
    log_message("  panel D v6 (curated dumbbell): no curated genes; skipping.")
    return(invisible())
  }
  long$substate_label <- factor(long$substate_label, levels = focus_labels)
  df <- long[long$pc == "PC1", , drop = FALSE]
  if (nrow(df) == 0) return(invisible())

  summ <- df |>
    dplyr::group_by(.data$category, .data$gene) |>
    dplyr::summarise(lo   = min(.data$loading, na.rm = TRUE),
                     hi   = max(.data$loading, na.rm = TRUE),
                     mu   = mean(.data$loading, na.rm = TRUE),
                     flip = .data$lo < 0 & .data$hi > 0,
                     n_sub = dplyr::n(),
                     .groups = "drop") |>
    dplyr::filter(.data$n_sub >= 2)
  if (nrow(summ) == 0) {
    log_message("  panel D v6: no curated gene present in >=2 focus substates.")
    return(invisible())
  }
  summ <- summ |> dplyr::arrange(.data$category, .data$mu)
  summ$gene <- factor(summ$gene, levels = unique(summ$gene))
  df <- df[df$gene %in% levels(summ$gene), , drop = FALSE]
  df$gene <- factor(df$gene, levels = levels(summ$gene))

  sub_pal <- strip_pal[focus_labels]
  if (any(is.na(sub_pal))) {
    sub_pal <- grDevices::hcl.colors(length(focus_labels), palette = "Dark 3")
    names(sub_pal) <- focus_labels
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0,
                        linewidth = 0.3, color = "grey60") +
    ggplot2::geom_segment(data = summ,
                          ggplot2::aes(x = .data$lo, xend = .data$hi,
                                       y = .data$gene, yend = .data$gene,
                                       color = .data$flip),
                          linewidth = 0.7) +
    ggplot2::geom_point(data = df,
                        ggplot2::aes(.data$loading, .data$gene,
                                     fill = .data$substate_label),
                        shape = 21, size = 2.7, color = "grey20",
                        stroke = 0.3) +
    ggplot2::scale_color_manual(
      values = c(`TRUE` = "#C0392B", `FALSE` = "grey55"),
      labels = c(`TRUE` = "Yes",     `FALSE` = "No"),
      name   = "Sign flips across substates") +
    ggplot2::scale_fill_manual(values = sub_pal, name = "Substate") +
    ggh4x::facet_grid2(rows = ggplot2::vars(.data$category),
                       scales = "free_y", space = "free_y",
                       strip = ggh4x::strip_themed(
                         text_y = ggh4x::elem_list_text(face = "bold",
                                                         size = 9))) +
    ggplot2::labs(
      title    = "Curated PC1 loadings across focus substates",
      subtitle = paste0("Each row is one curated gene. Bar spans the min/max ",
                        "PC1 loading across substates; red = sign flips (the ",
                        "gene drives different directions in different ",
                        "substates). Positive = Viral-driving."),
      x = "PC1 loading", y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title         = ggplot2::element_text(face = "bold"),
                   plot.subtitle      = ggplot2::element_text(size = 9),
                   panel.grid.minor   = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   strip.text.y       = ggplot2::element_text(angle = 0,
                                                              hjust = 0,
                                                              face = "bold"),
                   legend.position    = "bottom",
                   legend.box         = "vertical")

  n_genes <- nlevels(summ$gene)
  save_pdf_png(p, file.path(.myeloid_panelD_paths(paths),
                            "v6_pseudobulk_curated_dumbbell_PC1"),
               w = 9, h = max(6, 0.24 * n_genes + 2.5))
}

.myeloid_panelD_deg_counts <- function(paths, cfg) {
  dge <- read_per_substate_dge(paths, contrast_name = "Autoimmune_vs_Viral")
  if (is.null(dge) || nrow(dge) == 0) {
    log_message("  panel D v1 (DEG counts): DGE table not found; skipping.")
    return(invisible())
  }
  padj_thr <- cfg$dge$padj_threshold %||% 0.05
  lfc_thr  <- cfg$dge$lfc_threshold  %||% 0.5
  pf <- "padj" %in% colnames(dge)
  if (!pf || !"log2FoldChange" %in% colnames(dge)) {
    log_message("  panel D v1: DGE table missing padj/log2FoldChange.")
    return(invisible())
  }
  df <- dge |>
    dplyr::filter(!is.na(.data$padj), !is.na(.data$log2FoldChange),
                  .data$padj < padj_thr,
                  abs(.data$log2FoldChange) > lfc_thr) |>
    dplyr::mutate(direction = dplyr::if_else(.data$log2FoldChange > 0,
                                             "AI-up", "Viral-up")) |>
    dplyr::count(.data$substate, .data$direction, name = "n_sig")
  if (nrow(df) == 0) {
    log_message("  panel D v1: no substate passed thresholds.")
    return(invisible())
  }
  df$substate_label <- substate_labels(cfg, "myeloid", df$substate)
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(.data$n_sig,
                                    forcats::fct_reorder(.data$substate_label,
                                                          .data$n_sig, sum),
                                    fill = .data$direction)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(values = c("AI-up"    = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
                                          "Viral-up" = unname(ETIOLOGY_GROUP_COLORS["Viral"]))) +
    ggplot2::labs(title = "Myeloid substate DEG counts: AI (NIU) vs Viral",
                  subtitle = sprintf("padj < %.2g, |log2FC| > %.2g",
                                     padj_thr, lfc_thr),
                  x = "Significant genes", y = NULL, fill = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  save_pdf_png(p, file.path(.myeloid_panelD_paths(paths),
                            "v1_deg_counts_by_substate"),
               w = 7, h = max(4, 0.35 * length(unique(df$substate)) + 2))
}

.myeloid_panelD_cor_distance <- function(obj, paths, cfg) {
  pbs <- build_per_substate_pseudobulks(obj)
  if (length(pbs) == 0) {
    log_message("  panel D v2 (cor distance): no pseudobulks; skipping.")
    return(invisible())
  }
  rows <- lapply(names(pbs), function(ck) {
    pb <- pbs[[ck]]
    if (inherits(pb, "SummarizedExperiment")) {
      cd <- SummarizedExperiment::colData(pb)
      m  <- SummarizedExperiment::assay(pb, "counts")
    } else {
      cd <- pb$coldata; m <- pb$counts
    }
    grp <- as.character(cd$group)
    ai_cols  <- which(grp == "NIU")
    vir_cols <- which(grp == "Viral")
    if (length(ai_cols) < 2 || length(vir_cols) < 2) return(NULL)
    # Log-normalize sample-wise so the correlation isn't dominated by depth.
    norm <- function(x) {
      x <- as.matrix(x)
      sf <- pmax(colSums(x), 1)
      log1p(t(t(x) / sf) * 1e4)
    }
    nm <- norm(m)
    ai  <- rowMeans(nm[, ai_cols,  drop = FALSE])
    vir <- rowMeans(nm[, vir_cols, drop = FALSE])
    use <- ai > 0 | vir > 0
    if (sum(use) < 200) return(NULL)
    data.frame(substate = ck,
               cor_dist = 1 - stats::cor(ai[use], vir[use],
                                          method = "spearman"),
               n_ai = length(ai_cols), n_vir = length(vir_cols),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(df) || nrow(df) == 0) {
    log_message("  panel D v2: no substate had >=2 samples per arm.")
    return(invisible())
  }
  df$substate_label <- substate_labels(cfg, "myeloid", df$substate)
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(.data$cor_dist,
                                    forcats::fct_reorder(.data$substate_label,
                                                          .data$cor_dist))) +
    ggplot2::geom_col(fill = "grey40") +
    ggplot2::labs(title = "Myeloid AI vs Viral transcriptional distance",
                  subtitle = "1 - Spearman r on mean pseudobulk per substate",
                  x = "1 - Spearman r (higher = more AI/Viral divergence)",
                  y = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  save_pdf_png(p, file.path(.myeloid_panelD_paths(paths),
                            "v2_cor_distance_by_substate"),
               w = 7, h = max(4, 0.35 * nrow(df) + 2))
}

.myeloid_panelD_pca_facets <- function(obj, paths, cfg) {
  pbs <- build_per_substate_pseudobulks(obj)
  if (length(pbs) == 0) {
    log_message("  panel D v3 (PCA facets): no pseudobulks; skipping.")
    return(invisible())
  }
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    log_message("  panel D v3: DESeq2 not available; skipping.")
    return(invisible())
  }
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    log_message("  panel D v3: ggh4x not installed; installing for colored strips.")
    utils::install.packages("ggh4x")
  }

  # Per-substate PCA: refactored 2026-05-19 to call compute_per_substate_pca()
  # from R/45_compartment_pca.R so F3 panel F and F4 panels D/E share a single
  # implementation. The function is sign-flip-equivalent to the prior inline
  # block (PC1_oriented = pc$x[,1] * flip), so downstream plots reproduce.
  # Side-effect: also persists pca_subject_scores.csv / pca_gene_loadings.csv /
  # pca_variance_explained.csv / pca_pc1_significance.csv under
  # outputs/tables/eye/myeloid/ — same CSVs the F4 cross-compartment bridge
  # reads, so running F3 viz alone is now sufficient to refresh the F4 input.
  pca_res <- compute_per_substate_pca(
    obj,
    min_cells_per_pb = 10L,            # match historical F3 floor
    min_gene_count   = 10L,
    hvg_n            = 2000L,
    n_pcs            = 5L,
    vst_blind        = FALSE,
    pc1_split_fdr    = 0.05
  )
  if (is.null(pca_res)) {
    log_message("  panel D v3 (PCA facets): compute_per_substate_pca returned NULL.")
    return(invisible())
  }

  # Persist CSVs (the canonical F4 input). Idempotent; safe to rerun.
  ensure_dir(paths$results_tables)
  utils::write.csv(pca_res$scores,
                   file.path(paths$results_tables, "pca_subject_scores.csv"),
                   row.names = FALSE)
  utils::write.csv(pca_res$loadings,
                   file.path(paths$results_tables, "pca_gene_loadings.csv"),
                   row.names = FALSE)
  utils::write.csv(pca_res$variance,
                   file.path(paths$results_tables, "pca_variance_explained.csv"),
                   row.names = FALSE)
  utils::write.csv(pca_res$significance,
                   file.path(paths$results_tables, "pca_pc1_significance.csv"),
                   row.names = FALSE)

  # Rebuild the per_sub list-of-lists shape the downstream scores + loadings
  # plots expect. PC1 == PC1_oriented (post-flip) per the historical contract.
  scores_by_sub   <- split(pca_res$scores,   pca_res$scores$substate)
  loadings_by_sub <- split(pca_res$loadings, pca_res$loadings$substate)
  variance_by_sub <- split(pca_res$variance, pca_res$variance$substate)
  per_sub <- lapply(names(scores_by_sub), function(ck) {
    sc <- scores_by_sub[[ck]]
    ld <- loadings_by_sub[[ck]]
    vr <- variance_by_sub[[ck]]
    if (is.null(sc) || is.null(ld) || is.null(vr) || nrow(sc) == 0L) return(NULL)
    list(
      substate = ck,
      scores   = data.frame(sample   = sc$sample,
                            PC1      = sc$PC1_oriented,
                            PC2      = sc$PC2,
                            etiology = sc$Phenotype_2,
                            stringsAsFactors = FALSE),
      loadings = data.frame(gene = ld$gene,
                            PC1  = ld$PC1_oriented,
                            PC2  = ld$PC2,
                            stringsAsFactors = FALSE),
      var_pct  = vr$var_explained
    )
  })
  per_sub <- per_sub[!vapply(per_sub, is.null, logical(1))]
  if (length(per_sub) == 0) {
    log_message("  panel D v3: no substate produced a PCA.")
    return(invisible())
  }

  # Per-substate palette for facet strips. hcl.colors palette is colorblind
  # friendlier than the default; alpha-mute slightly so dark text on the
  # strip stays legible.
  sub_ids        <- vapply(per_sub, `[[`, character(1), "substate")
  sub_label_vec  <- substate_labels(cfg, "myeloid", sub_ids)
  ord <- order(sub_label_vec)
  sub_ids       <- sub_ids[ord]
  sub_label_vec <- sub_label_vec[ord]
  per_sub       <- per_sub[ord]
  strip_pal     <- grDevices::hcl.colors(length(sub_label_vec),
                                          palette = "Dark 3")
  names(strip_pal) <- sub_label_vec
  fill_pal <- vapply(strip_pal,
                     function(h) {
                       rgb <- grDevices::col2rgb(h) / 255
                       grDevices::rgb(rgb[1], rgb[2], rgb[3], alpha = 0.55)
                     }, character(1))

  facet_strip <- ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = fill_pal),
    text_x       = ggh4x::elem_list_text(face = "bold", size = 9))

  # ---- 1. Scores plot: per-substate sample positions in PC space ----------
  scores_df <- do.call(rbind, lapply(per_sub, function(x) {
    df <- x$scores
    df$substate       <- x$substate
    df$substate_label <- substate_labels(cfg, "myeloid", x$substate)
    df$var_pc1        <- x$var_pct[1]
    df$var_pc2        <- x$var_pct[2]
    df
  }))
  scores_df$substate_label <- factor(scores_df$substate_label,
                                     levels = sub_label_vec)
  scores_df$etiology_label <- ifelse(scores_df$etiology == "NIU",
                                     "Autoimmune (NIU)", "Viral")

  p_scores <- ggplot2::ggplot(scores_df,
                              ggplot2::aes(.data$PC1, .data$PC2,
                                           color = .data$etiology_label)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        linewidth = 0.25, color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        linewidth = 0.25, color = "grey70") +
    ggplot2::geom_point(size = 2.2, alpha = 0.9) +
    ggplot2::stat_ellipse(level = 0.7,
                          ggplot2::aes(group = .data$etiology_label),
                          linewidth = 0.5) +
    ggplot2::scale_color_manual(values = c(
      "Autoimmune (NIU)" = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
      "Viral"            = unname(ETIOLOGY_GROUP_COLORS["Viral"]))) +
    ggh4x::facet_wrap2(~ .data$substate_label, scales = "free", ncol = 2,
                       strip = facet_strip) +
    ggplot2::labs(title = "Per-substate pseudobulk PCA (DESeq2 vst)",
                  subtitle = "PC1 oriented so Viral centroid is positive",
                  x = "PC1", y = "PC2", color = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title       = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank(),
                   aspect.ratio     = 1)

  n_facets <- length(per_sub)
  # Match canvas to a 2-col grid of square panels so theme(aspect.ratio = 1)
  # isn't stretched by the device.
  n_cols   <- 2
  panel_in <- 3.4
  n_rows   <- ceiling(n_facets / n_cols)
  save_pdf_png(p_scores, file.path(.myeloid_panelD_paths(paths),
                                   "v3_pseudobulk_pca_facets"),
               w = n_cols * panel_in + 2.5,
               h = n_rows * panel_in + 2)

  # ---- 2. Loadings plot: signed bars of top genes on PC1 and PC2 ----------
  # Restrict to the substates with clean AI vs Viral separation on PC1
  # (FOLR2+ TRM mac, Macrophage, Langerhans/cDC2, cDC1) so the panel is
  # interpretable rather than dominated by underpowered substates. Drop
  # genes that are stress / ribosomal / mito / pseudo-gene noise rather
  # than biology — these are common scRNA loading artifacts on small
  # pseudobulks. Override either list via cfg$visualization$panelD_loadings
  # when you want a different focus.
  focus_subs <- cfg$visualization$panelD_loadings$substates_keep %||%
                c("1", "3", "5", "7")
  drop_regex <- cfg$visualization$panelD_loadings$gene_drop_regex %||%
                paste0("^(HSP|DNAJ|LINC|RPS|RPL|MRPS|MRPL|MT-|MTRNR|",
                       "AC\\d|AL\\d|AP\\d|AF\\d|ENSG|XIST|MALAT1|NEAT1)")
  per_sub_focus <- per_sub[vapply(per_sub,
                                   function(x) x$substate %in% focus_subs,
                                   logical(1))]
  if (length(per_sub_focus) == 0) {
    log_message("  panel D v3 loadings: no focus substates ('",
                paste(focus_subs, collapse = ","),
                "') survived pseudobulk filter; skipping loadings plot.")
    return(invisible())
  }

  # Pick top n positive + n negative by signed loading, then also pull in any
  # HLA or T-cell-interaction gene that's within `priority_pool` of the top
  # by |loading|. Caps priority pickups at `priority_cap` so the panel
  # height stays bounded.
  pick_top <- function(df, col, n = 6, priority_pool = 60,
                       priority_cap = 4) {
    pos <- df |> dplyr::slice_max(.data[[col]], n = n, with_ties = FALSE)
    neg <- df |> dplyr::slice_min(.data[[col]], n = n, with_ties = FALSE)
    base <- dplyr::bind_rows(pos, neg)
    pool <- df |>
      dplyr::mutate(.absload = abs(.data[[col]])) |>
      dplyr::slice_max(.data$.absload, n = priority_pool, with_ties = FALSE)
    flagged <- pool[.gene_category(pool$gene) %in% c("HLA", "Tcell"), ]
    flagged <- flagged |>
      dplyr::slice_max(.data$.absload, n = priority_cap, with_ties = FALSE) |>
      dplyr::select(-".absload")
    dplyr::bind_rows(base, flagged) |>
      dplyr::distinct(.data$gene, .keep_all = TRUE)
  }
  top_df <- do.call(rbind, lapply(per_sub_focus, function(x) {
    df <- x$loadings
    df <- df[!grepl(drop_regex, df$gene), , drop = FALSE]
    df$substate_label <- substate_labels(cfg, "myeloid", x$substate)
    pc1 <- pick_top(df, "PC1") |>
      dplyr::mutate(pc = "PC1 (AI <-> Viral)", loading = .data$PC1)
    pc2 <- pick_top(df, "PC2") |>
      dplyr::mutate(pc = "PC2 (within-arm spread)", loading = .data$PC2)
    dplyr::bind_rows(pc1, pc2) |>
      dplyr::select("substate_label", "pc", "gene", "loading")
  }))
  focus_labels <- substate_labels(cfg, "myeloid", focus_subs)
  focus_labels <- focus_labels[focus_labels %in% top_df$substate_label]
  top_df$substate_label <- factor(top_df$substate_label,
                                  levels = focus_labels)
  top_df$pc             <- factor(top_df$pc,
                                  levels = c("PC1 (AI <-> Viral)",
                                             "PC2 (within-arm spread)"))
  # Per-facet gene ordering: encode the panel into the y-value so each
  # (substate x pc) panel sorts independently. Strip the prefix at axis-
  # rendering time so the visible label is just the gene symbol.
  top_df <- top_df |>
    dplyr::group_by(.data$substate_label, .data$pc) |>
    dplyr::arrange(.data$loading, .by_group = TRUE) |>
    dplyr::mutate(gene_key = paste(.data$substate_label, .data$pc,
                                   .data$gene, sep = "||")) |>
    dplyr::ungroup()
  top_df$gene_key <- factor(top_df$gene_key, levels = top_df$gene_key)
  top_df$bar_fill <- ifelse(
    top_df$pc == "PC1 (AI <-> Viral)",
    ifelse(top_df$loading > 0, "Viral-driving",      "AI-driving"),
    ifelse(top_df$loading > 0, "PC2 positive end",   "PC2 negative end"))
  bar_pal <- c(
    "AI-driving"        = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
    "Viral-driving"     = unname(ETIOLOGY_GROUP_COLORS["Viral"]),
    "PC2 positive end"  = "#5D6D7E",
    "PC2 negative end"  = "#AEB6BF")

  if (!requireNamespace("ggtext", quietly = TRUE)) {
    log_message("  panel D v3 loadings: ggtext not installed; installing for ",
                "colored gene labels.")
    utils::install.packages("ggtext")
  }

  hla_color   <- "#B9770E"
  tcell_color <- "#117864"
  subtitle_md <- paste0(
    "Top 6 +/- 6 genes per (substate, PC); positive PC1 = Viral-driving. ",
    "Axis labels colored: ",
    "<span style='color:", hla_color,   ";font-weight:bold'>HLA / antigen processing</span>, ",
    "<span style='color:", tcell_color, ";font-weight:bold'>T cell interaction (cytokine/chemokine/costim/checkpoint)</span>. ",
    "Filtered: HSP / DNAJ / LINC / RPS / RPL / MRP / MT- / Ensembl clone IDs / XIST / MALAT1 / NEAT1.")

  p_load <- ggplot2::ggplot(top_df,
                            ggplot2::aes(.data$loading, .data$gene_key,
                                         fill = .data$bar_fill)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
    ggplot2::geom_col() +
    ggplot2::scale_y_discrete(labels = .panelD_loading_axis_labels) +
    ggplot2::scale_fill_manual(values = bar_pal, name = NULL) +
    ggh4x::facet_grid2(rows = ggplot2::vars(.data$substate_label),
                        cols = ggplot2::vars(.data$pc),
                        scales = "free", independent = "all",
                        strip = ggh4x::strip_themed(
                          background_y = ggh4x::elem_list_rect(
                            fill = unname(fill_pal[focus_labels])),
                          text_y       = ggh4x::elem_list_text(face = "bold",
                                                                size = 10))) +
    ggplot2::labs(title = "Per-substate PCA gene loadings",
                  subtitle = subtitle_md,
                  x = "Loading", y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title         = ggplot2::element_text(face = "bold"),
                   plot.subtitle      = ggtext::element_markdown(size = 9,
                                                                  lineheight = 1.2),
                   panel.grid.minor   = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   strip.text.y       = ggplot2::element_text(angle = 0),
                   legend.position    = "bottom",
                   axis.text.y        = ggtext::element_markdown(size = 9))

  n_focus <- length(focus_labels)
  save_pdf_png(p_load, file.path(.myeloid_panelD_paths(paths),
                                 "v3_pseudobulk_pca_loadings"),
               w = 10, h = max(5, n_focus * 1.8 + 1.5))

  # ---- 3. Curated HLA / T cell views (v4 heatmap, v5 lollipop, v6 dumbbell)
  # Same focus substates and same per-PCA loadings the v3 stack used; just
  # reshaped to put HLA / antigen-processing / cytokine / chemokine /
  # costim / checkpoint biology on the y-axis instead of letting it compete
  # with generic top-loading genes inside each substate facet.
  .myeloid_panelD_v4_curated_heatmap(per_sub_focus, paths, cfg, focus_labels)
  .myeloid_panelD_v5_curated_lollipop(per_sub_focus, paths, cfg,
                                       focus_labels, strip_pal)
  .myeloid_panelD_v6_curated_dumbbell(per_sub_focus, paths, cfg,
                                       focus_labels, strip_pal)
}

# ---------------------------------------------------------------------------
# Panel F: curated APC UCell modules (escape.APC assay produced by
# run_escape_custom_modules in 32_escape.R). Left = module x substate dotplot
# faceted by etiology (size = fraction above threshold, color = mean UCell).
# Right = subject-mean scatter of MHC-I vs MHC-II module scores with etiology
# ellipses. Both use the same UCell scale so reviewers can't accuse the
# panel of double-scoring. Also writes a Wilcoxon test of subject-mean
# MHC-I and MHC-II between AI and Viral to APC_module_etiology_test.csv.
# ---------------------------------------------------------------------------
.myeloid_apc_get_scores <- function(obj, cfg) {
  spec <- cfg$escape$custom_modules$myeloid_apc
  assay_name <- spec$new_assay %||% "escape.APC"
  if (!assay_name %in% Seurat::Assays(obj)) {
    log_message("  APC: assay '", assay_name,
                "' not on myeloid object. Run run_escape_custom_modules first.")
    return(NULL)
  }
  mat <- as.matrix(Seurat::GetAssayData(obj, assay = assay_name, layer = "data"))
  list(mat = mat, threshold = spec$ucell_dotplot_threshold %||% 0.1)
}

# Seurat's CreateAssayObject converts underscores in feature (module) names to
# dashes. Map a config-style name like "MHC_II_CLASSICAL" to whichever spelling
# is on the assay rownames, returning NA when neither form is present.
.apc_resolve_module <- function(name, available) {
  if (name %in% available) return(name)
  alt <- gsub("_", "-", name, fixed = TRUE)
  if (alt %in% available) return(alt)
  alt2 <- gsub("-", "_", name, fixed = TRUE)
  if (alt2 %in% available) return(alt2)
  NA_character_
}

.myeloid_apc_dotplot <- function(obj, paths, cfg) {
  sc <- .myeloid_apc_get_scores(obj, cfg)
  if (is.null(sc)) return(invisible())
  meta <- obj[[]]
  if (!"Phenotype_2" %in% colnames(meta)) return(invisible())
  meta$substate <- substate_labels(cfg, "myeloid", meta$knn.leiden.cluster)
  cells <- intersect(colnames(sc$mat), rownames(meta))
  if (length(cells) < 50) return(invisible())
  mat  <- sc$mat[, cells, drop = FALSE]
  meta <- meta[cells, , drop = FALSE]
  meta <- meta[meta$Phenotype_2 %in% c("NIU", "Viral"), , drop = FALSE]
  if (nrow(meta) < 50) return(invisible())
  mat  <- mat[, rownames(meta), drop = FALSE]

  long <- lapply(rownames(mat), function(mod) {
    data.frame(module = mod, score = mat[mod, ],
               substate = meta$substate, etiology = meta$Phenotype_2,
               stringsAsFactors = FALSE)
  })
  long <- do.call(rbind, long)
  thr <- sc$threshold
  summary <- long |>
    dplyr::group_by(.data$module, .data$substate, .data$etiology) |>
    dplyr::summarise(mean_score = mean(.data$score, na.rm = TRUE),
                     frac_above = mean(.data$score > thr, na.rm = TRUE),
                     n_cells    = dplyr::n(),
                     .groups = "drop") |>
    dplyr::filter(.data$n_cells >= 10)
  if (nrow(summary) == 0) return(invisible())
  summary$etiology_label <- ifelse(summary$etiology == "NIU",
                                   "Autoimmune (NIU)", "Viral")

  p <- ggplot2::ggplot(summary,
                       ggplot2::aes(.data$substate, .data$module,
                                    size = .data$frac_above,
                                    color = .data$mean_score)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ .data$etiology_label) +
    ggplot2::scale_color_viridis_c(option = "viridis") +
    ggplot2::scale_size_continuous(range = c(0.5, 7),
                                   labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(title = "Myeloid APC modules (escape.APC, UCell)",
                  subtitle = sprintf("dot size = fraction of cells with UCell > %.2f",
                                     thr),
                  x = NULL, y = NULL,
                  color = "Mean UCell", size = "Frac > threshold") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                   plot.title  = ggplot2::element_text(face = "bold"))
  save_pdf_png(p, file.path(viz_subdir(paths, "escape"), "myeloid_apc_modules_dotplot"),
               w = 11, h = 4 + 0.3 * nrow(mat))
}

# ---------------------------------------------------------------------------
# Panel G: myeloid -> T cell signalling.
# Primary: read outputs/tables/cross_compartment/liana_myeloid_to_tcell_combined.csv
# (written by R/47_liana_myeloid_tcell) and plot top 20 LR pairs by
# absolute disease_bias (negative = NIU-driving, positive = Viral-driving).
# Fallback: T-cell-relevant ligand expression heatmap on the myeloid object.
# Migrated from CellChat (R/64_cellchat.R, retired 2026-05-19) to the LIANA
# consensus rank — same panel framing, methods-paper-stronger inference.
# ---------------------------------------------------------------------------
.myeloid_panelG_liana <- function(paths, cfg) {
  cc_paths <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/")
  combined_path <- file.path(cc_paths$tables,
                             "liana_myeloid_to_tcell_combined.csv")
  if (!file.exists(combined_path)) return(FALSE)
  df <- tryCatch(utils::read.csv(combined_path, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L ||
      !all(c("source", "target", "ligand_complex", "receptor_complex",
             "disease_bias", "unique_to") %in% colnames(df))) return(FALSE)

  df <- df |>
    dplyr::filter(.data$unique_to == "both",
                  !is.na(.data$disease_bias)) |>
    dplyr::mutate(abs_bias = abs(.data$disease_bias)) |>
    dplyr::slice_max(.data$abs_bias, n = 20)
  if (nrow(df) == 0L) return(FALSE)

  df$lr_pair <- paste(df$ligand_complex, df$receptor_complex, sep = " -> ")
  df$direction <- ifelse(df$disease_bias < 0, "NIU-driving", "Viral-driving")
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(.data$disease_bias,
                                    forcats::fct_reorder(.data$lr_pair,
                                                         .data$disease_bias),
                                    fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(
      `NIU-driving`   = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
      `Viral-driving` = unname(ETIOLOGY_GROUP_COLORS["Viral"])),
      guide = "none") +
    ggplot2::labs(title = "Myeloid -> T cell signalling (LIANA consensus)",
                  subtitle = "Top 20 LR pairs by |rank_Viral - rank_NIU|",
                  x = "disease_bias  (NIU <-> Viral)", y = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  save_pdf_png(p, file.path(viz_subdir(paths, "pca_coupling"),
                            "myeloid_to_tcell_liana_disease_bias"),
               w = 8, h = 7)
  TRUE
}

.myeloid_panelG_ligand_heatmap <- function(obj, paths, cfg) {
  ligands <- c("IL12B", "IL15", "IL18", "IL23A", "IL27",
               "CXCL9", "CXCL10", "CXCL11",
               "TNFSF10", "TNFSF13B", "CD70", "CD80", "CD86", "ICOSLG")
  ligands <- intersect(ligands, rownames(obj))
  if (length(ligands) < 3) {
    log_message("  panel G fallback: <3 T-cell ligands detected; skipping.")
    return(invisible())
  }
  meta <- obj[[]]
  if (!"Phenotype_2" %in% colnames(meta)) return(invisible())
  meta <- meta[meta$Phenotype_2 %in% c("NIU", "Viral"), , drop = FALSE]
  meta$substate <- substate_labels(cfg, "myeloid", meta$knn.leiden.cluster)
  meta$col      <- paste(meta$substate,
                         ifelse(meta$Phenotype_2 == "NIU",
                                "AI", "Viral"), sep = " | ")
  cells <- rownames(meta)
  # Keep the expression matrix sparse — only subset rows (ligands are a
  # handful) and let Matrix do the column-group means. Coercing to dense
  # over the full 12k cells x 30k genes object allocates several GiB.
  expr <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  expr <- expr[ligands, cells, drop = FALSE]
  col_grp <- factor(meta$col[match(colnames(expr), rownames(meta))])
  # Sparse mean per group = sparse sum %*% indicator / n
  groups <- levels(col_grp)
  n_per  <- as.integer(table(col_grp)[groups])
  ind <- Matrix::sparseMatrix(i = seq_along(col_grp),
                              j = as.integer(col_grp),
                              x = 1,
                              dims = c(length(col_grp), length(groups)))
  mat <- as.matrix((expr %*% ind) %*% Matrix::Diagonal(x = 1 / n_per))
  colnames(mat) <- groups
  rownames(mat) <- ligands

  summary <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(summary) <- c("gene", "col", "mean_expr")
  long_z <- summary |>
    dplyr::group_by(.data$gene) |>
    dplyr::mutate(z = (.data$mean_expr - mean(.data$mean_expr, na.rm = TRUE)) /
                       (stats::sd(.data$mean_expr, na.rm = TRUE) + 1e-6)) |>
    dplyr::ungroup()
  p <- ggplot2::ggplot(long_z,
                       ggplot2::aes(.data$col, .data$gene, fill = .data$z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "viridis") +
    ggplot2::labs(title = "Myeloid T-cell ligand expression (panel G fallback)",
                  subtitle = "Mean log-normalized expression, z-scored per gene",
                  x = NULL, y = NULL, fill = "z") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                   plot.title  = ggplot2::element_text(face = "bold"))
  save_pdf_png(p, file.path(viz_subdir(paths, "pca_coupling"),
                            "myeloid_to_tcell_ligand_heatmap_fallback"),
               w = max(8, ncol(mat) * 0.35 + 3),
               h = 4 + length(ligands) * 0.25)
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
run_visualizations_myeloid <- function(cfg) {
  paths <- get_target_paths(cfg, "myeloid")
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Myeloid IntegratedSeuratObject.rds not found at ", obj_path,
                ". Skipping myeloid viz.")
    return(invisible(TRUE))
  }
  log_message("=== Myeloid compartment visualizations ===")
  ensure_dir(paths$viz_dir)
  obj <- readRDS(obj_path)

  # --- Shared cross-cutting panels (helpers in 81_viz_compartment_helpers.R)
  viz_compartment_umap(obj, "myeloid", paths, cfg)
  viz_compartment_umap_etiology(obj, "myeloid", paths, cfg)
  viz_compartment_dotplot(obj, "myeloid", paths, cfg)
  viz_compartment_milo_da_box("myeloid", paths, cfg)
  viz_compartment_composition(obj, "myeloid", paths, cfg)
  viz_compartment_volcano("myeloid", paths, cfg)
  viz_compartment_volcano_per_substate("myeloid", paths, cfg)
  viz_compartment_niu_subcontrast_heatmap("myeloid", paths, cfg)
  viz_compartment_functional_gestalt_full_heatmap(obj, "myeloid", paths, cfg)

  # --- Myeloid-specific module score: Complement
  myeloid_modules <- list(
    Complement = c("C1QA", "C1QB", "C1QC", "C3", "C2", "CFB", "C3AR1",
                   "C5AR1", "SERPING1")
  )
  viz_compartment_module_score_niu_vs_viral(obj, "myeloid", paths, cfg,
                                            myeloid_modules)

  # --- Myeloid-specific UMAP overlays
  .myeloid_marker_featureplots(obj, paths, cfg)

  # --- Panel D variants (DEG counts / cor distance / per-substate PCA)
  .myeloid_panelD_deg_counts(paths, cfg)
  .myeloid_panelD_cor_distance(obj, paths, cfg)
  .myeloid_panelD_pca_facets(obj, paths, cfg)

  # --- Panel F: curated APC UCell modules (dotplot + scatter + Wilcoxon)
  .myeloid_apc_dotplot(obj, paths, cfg)

  # --- Panel G: myeloid -> T cell signalling. Prefers LIANA combined CSV
  # (R/47_liana_myeloid_tcell); falls back to the curated T-cell-relevant
  # ligand heatmap if the CSV is absent so the panel slot has content
  # either way.
  if (!isTRUE(.myeloid_panelG_liana(paths, cfg))) {
    .myeloid_panelG_ligand_heatmap(obj, paths, cfg)
  }

  # --- Full MiloR viz block (beeswarm/boxplot/summary/composition/volcano/
  # nhoodgraph/nhood-size hist) for myeloid
  if (exists("viz_milo")) {
    tryCatch(viz_milo(obj, cfg, paths, target = "myeloid"),
             error = function(e)
               log_message("  viz_milo failed for myeloid: ",
                           conditionMessage(e)))
  }

  log_message("=== Myeloid visualizations complete ===")
  invisible(TRUE)
}

# run_visualizations_myeloid renders the myeloid compartment figure block.
# Layout: cross-cutting cluster panels (substate UMAP labeled+stripped, UMAP
# by etiology, marker dot plot, NIU-vs-Viral composition, compartment-global
# and per-substate volcanoes, GSEA pathway-by-substate heatmap, NIU-vs-Viral
# top pathway bar, pathway-by-etiology heatmap, NIU sub-contrast heatmap),
# then myeloid-specific overlays (lineage marker featureplots, HLA class I
# and class II expression, top-pathway UCell UMAP, ISG/Type17/Complement
# module scores split by NIU vs Viral, dedicated NIU vs VZV vs other-viral
# ISG check). The full MiloR viz block (drawn from R/82_viz_dispatch.R) is
# run last so the myeloid milo nhood-graph, beeswarm, boxplot, summary,
# composition, volcano, and nhood-size hist all land under
# outputs/viz/eye/myeloid/07_milo/. All filenames are snake_case
# descriptions; no F3 panel-letter scheme is used.
