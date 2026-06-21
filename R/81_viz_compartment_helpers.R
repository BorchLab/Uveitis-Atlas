# R/81_viz_compartment_helpers.R
# Shared visualization helpers used by R/85_viz_myeloid.R, R/86_viz_bcell.R,
# and R/88_viz_tcell.R. Pulled out of the per-compartment files so the cross-
# cutting panels (composition, volcano, GSEA heatmap, NIU-vs-Viral pathway
# bar, NIU sub-contrast heatmap, per-substate volcano grid, repertoire panels)
# stay in one place.
#
# Aesthetic conventions:
#   - Continuous color: viridis everywhere
#   - NIU red (#E21F26), Viral blue (#397FB9), Healthy teal (#A8DADC)
#   - dual_save_plot for UMAP-style panels (labeled + stripped)
#   - save_pdf_png for ggplot panels (one labeled PDF + PNG)
#   - File stems are descriptive: "<compartment>_<what>_<grouping>" (snake_case)

suppressPackageStartupMessages({
  library(Seurat)
  library(scplotter)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(viridis)
})

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
suppressPackageStartupMessages(library(ggrepel))

ETIOLOGY_GROUP_COLORS <- c(NIU = "#E21F26", Viral = "#397FB9",
                           Healthy = "#A8DADC")

# Shared palette for VDJdb antigen species so every pathogen-keyed plot
# (UMAP highlights, fraction bars, motif logos) uses a single color per
# species. Hexes chosen for distinguishability on greyed-out backgrounds
# and at small dot sizes; "Other" falls through to grey85 by convention.
# Add new species here rather than ad-hoc in any caller.
PATHOGEN_COLORS <- c(
  CMV           = "#1B7837",  # green
  EBV           = "#762A83",  # purple
  InfluenzaA    = "#E08214",  # orange
  "SARS-CoV-2"  = "#B2182B",  # red
  HCV           = "#2166AC",  # blue
  "HIV-1"       = "#8C510A",  # brown
  YFV           = "#FFB300",  # amber
  MTuberculosis = "#01665E",  # teal
  "VZV"         = "#5AAE61",  # green-teal (HHV-3)
  "HSV-1"       = "#9970AB",  # muted purple (HHV-1)
  "HSV-2"       = "#C2A5CF",  # paler purple (HHV-2)
  "HTLV-1"      = "#3690C0",  # blue (retrovirus)
  HomoSapiens   = "#878787",  # grey (kept for back-compat; usually excluded)
  Other         = "#BABABA"
)

# Direction-keyed colours, derived from the canonical palette so volcanoes /
# pathway bars stay in sync if the hexes change.
ETIOLOGY_DIRECTION_COLORS <- c(NIU_up   = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
                               Viral_up = unname(ETIOLOGY_GROUP_COLORS["Viral"]),
                               ns       = "grey75")

# Per-subtype etiology palette. Viral subtypes share the cool family of
# ETIOLOGY_GROUP_COLORS["Viral"]; NIU subtypes share the warm family of
# ETIOLOGY_GROUP_COLORS["NIU"]. Use whenever a plot encodes individual
# etiologies (per-etiology facets, alluvials, subject jitter colored by
# Etiology) so the Viral/NIU phenotype is still readable at a glance.
ETIOLOGY_SUBTYPE_COLORS <- c(
  # Viral (blue family)
  VZV_ARN     = "#397FB9",
  HSV1        = "#1F5C92",
  HSV2        = "#4A93C7",
  CMV_CRN     = "#2E6FA1",
  HTLV1       = "#62A2D4",
  # NIU (red family)
  Idiopathic  = "#E21F26",
  HLA_B27     = "#9D0208",
  VKH         = "#C61D27",
  BSCR        = "#F25A60",
  JIA         = "#B11820",
  SLE         = "#ED8186",
  Phakic_lens = "#7C0509",
  # Healthy
  Healthy     = "#A8DADC"
)

# ---------------------------------------------------------------------------
# Compartment marker dot plot (scplotter, viridis) — replaces the old
# F<n>_panelB_dotplot. Filename: <cmp>_substate_marker_dotplot.pdf
# ---------------------------------------------------------------------------
viz_compartment_dotplot <- function(obj, cmp, paths, cfg) {
  top_path <- file.path(paths$results_tables,
                        paste0(cmp, "_top_markers.csv"))
  if (!file.exists(top_path)) {
    log_message("  marker dotplot: ", basename(top_path), " not found.")
    return(invisible())
  }
  top <- read.csv(top_path, stringsAsFactors = FALSE)
  if (nrow(top) == 0) return(invisible())
  genes <- intersect(unique(top$gene), rownames(obj))
  if (length(genes) == 0) return(invisible())
  obj$substate_label <- substate_labels(cfg, cmp, obj$knn.leiden.cluster)
  obj$substate_label <- factor(obj$substate_label, levels = sort(unique(obj$substate_label)))
  out <- file.path(paths$viz_dir,
                   paste0(cmp, "_substate_marker_dotplot.pdf"))
  safe_plot(bquote(
    FeatureStatPlot(obj, features = .(genes),
                    ident = "substate_label",
                    plot_type = "dot",
                    palette = "viridis",
                    row_name_annotation = FALSE)
  ), out, height = max(8, 0.4 * length(genes)), width = 6)
}

# ---------------------------------------------------------------------------
# UMAP colored by Leiden substate (one labeled + one stripped variant)
# ---------------------------------------------------------------------------
viz_compartment_umap <- function(obj, cmp, paths, cfg) {
  obj$substate_label <- substate_labels(cfg, cmp, obj$knn.leiden.cluster)
  obj$substate_label <- factor(obj$substate_label, levels = sort(unique(obj$substate_label)))
  umap_sz <- cfg$visualization$umap_size %||% 8
  stripped_bg <- cfg$visualization$stripped_bg %||% "black"
  base <- file.path(paths$viz_dir, paste0(cmp, "_substate_umap"))
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = "substate_label",
                reduction = "UMAP",
                label = TRUE, label_insitu = TRUE,
                show_stat = FALSE, raster = FALSE)
  ), base, width = umap_sz, height = umap_sz, stripped_bg = stripped_bg)
}

# ---------------------------------------------------------------------------
# UMAP colored by Etiology (per-cell) — quick check that NIU subjects
# don't park in one corner of the embedding (that would be a confounder)
# ---------------------------------------------------------------------------
viz_compartment_umap_etiology <- function(obj, cmp, paths, cfg) {
  if (!"Etiology" %in% colnames(obj[[]])) return(invisible())
  umap_sz <- cfg$visualization$umap_size %||% 8
  stripped_bg <- cfg$visualization$stripped_bg %||% "black"
  base <- file.path(paths$viz_dir, paste0(cmp, "_substate_umap_by_etiology"))
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = "Etiology", reduction = "UMAP",
                show_stat = FALSE, raster = FALSE)
  ), base, width = umap_sz, height = umap_sz, stripped_bg = stripped_bg)
}

# ---------------------------------------------------------------------------
# Composition stacked bar: each substate row is a fill (sums to 1), with the
# segments computed from *phenotype-normalised* counts. That is, for each
# (cluster, phenotype) we use prop = n_in_cluster / total_in_phenotype, then
# the fill geom rescales the (NIU, Viral) pair so it spans 0–1 per row.
# Net effect: each bar shows the share of the cluster's "phenotype-corrected"
# representation that comes from NIU vs Viral — so the 50% line is the
# baseline-no-DA reference (and is independent of the absolute NIU vs Viral
# cell-count imbalance in the compartment).
# Row order: Milo DA call when on disk (NIU-up on top, Viral-up at bottom);
# otherwise by NIU/Viral phenotype-fraction log2 ratio (most NIU-heavy on top).
# Filename: <cmp>_substate_composition_niu_vs_viral.pdf/png
# ---------------------------------------------------------------------------
viz_compartment_composition <- function(obj, cmp, paths, cfg) {
  meta <- obj@meta.data
  meta$substate <- substate_labels(cfg, cmp, meta$knn.leiden.cluster)
  meta <- meta |> dplyr::filter(Phenotype_2 %in% c("NIU", "Viral"))
  if (nrow(meta) == 0) return(invisible())

  # Phenotype-normalised: each (cluster, phenotype) count divided by phenotype
  # total. position="fill" below then rescales (NIU_prop, Viral_prop) within
  # each row to span 0..1, so the 50% line is the baseline-no-DA reference.
  pheno_totals <- meta |>
    dplyr::count(.data$Phenotype_2, name = "pheno_total")
  comp <- meta |>
    dplyr::count(.data$substate, .data$Phenotype_2) |>
    dplyr::left_join(pheno_totals, by = "Phenotype_2") |>
    dplyr::mutate(prop = .data$n / .data$pheno_total)

  per_cluster <- comp |>
    dplyr::select(.data$substate, .data$Phenotype_2, .data$prop) |>
    tidyr::pivot_wider(names_from = .data$Phenotype_2,
                       values_from = .data$prop,
                       values_fill = 0) |>
    dplyr::mutate(log_ratio = log2((.data$NIU + 1e-4) /
                                   (.data$Viral + 1e-4)),
                  n_cells = .data$NIU * pheno_totals$pheno_total[
                                       pheno_totals$Phenotype_2 == "NIU"] +
                            .data$Viral * pheno_totals$pheno_total[
                                       pheno_totals$Phenotype_2 == "Viral"])

  # Row ordering: Milo DA when available, else NIU/Viral log-ratio.
  da <- tryCatch(compartment_milo_da_clusters(cmp, paths, cfg),
                 error = function(e) NULL)
  if (!is.null(da) && nrow(da) > 0) {
    da_ord <- da |>
      dplyr::arrange(factor(.data$direction, levels = c("NIU_up", "Viral_up")),
                     dplyr::desc(.data$median_logFC))
    sub_levels <- c(da_ord$substate_label,
                    setdiff(as.character(per_cluster$substate),
                            da_ord$substate_label))
    order_subtitle <- "Rows ordered by Milo DA direction (NIU-up on top)"
  } else {
    sub_levels <- per_cluster |>
      dplyr::arrange(dplyr::desc(.data$log_ratio)) |>
      dplyr::pull(.data$substate) |>
      as.character()
    order_subtitle <- "Rows ordered by NIU/Viral phenotype-fraction log2 ratio (most NIU-heavy on top)"
  }
  comp$substate <- factor(comp$substate, levels = rev(sub_levels))
  per_cluster$substate <- factor(per_cluster$substate,
                                 levels = rev(sub_levels))

  # Per-cluster propeller q for the right-side annotation.
  prop_path <- file.path(paths$results_tables,
                         "composition_test_Autoimmune_vs_Viral_knn.leiden.cluster.csv")
  q_lookup <- if (file.exists(prop_path)) {
    pp <- read.csv(prop_path, row.names = 1)
    if ("FDR" %in% colnames(pp)) setNames(pp$FDR, rownames(pp)) else NULL
  } else NULL
  ann <- per_cluster |>
    dplyr::mutate(cluster_id = sub(":.*", "", as.character(.data$substate)),
                  q = if (!is.null(q_lookup)) q_lookup[.data$cluster_id]
                      else NA_real_,
                  label = ifelse(is.na(.data$q),
                                 sprintf("log2 NIU/Viral=%.2f", .data$log_ratio),
                                 sprintf("q=%.2g, log2 NIU/Viral=%.2f",
                                         .data$q, .data$log_ratio)))

  p <- ggplot(comp, aes(y = .data$substate, x = .data$prop,
                        fill = .data$Phenotype_2)) +
       geom_col(position = "fill", width = 0.78) +
       geom_vline(xintercept = 0.5, color = "grey25", linewidth = 0.4,
                  linetype = "dashed") +
       geom_text(data = ann,
                 aes(x = 1.02, y = .data$substate, label = .data$label),
                 hjust = 0, size = 3, inherit.aes = FALSE) +
       scale_fill_manual(values = ETIOLOGY_GROUP_COLORS, name = "Phenotype") +
       scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                          labels = scales::percent_format(accuracy = 1),
                          expand = expansion(add = c(0, 0.3))) +
       labs(title = paste0(tools::toTitleCase(cmp),
                           " substate composition: NIU vs Viral"),
            subtitle = paste0(order_subtitle,
                              ". Bars rescale phenotype-normalised fractions to fill; 50% = no DA."),
            x = "Phenotype-normalised share (NIU vs Viral within cluster)",
            y = NULL) +
       theme_classic(base_size = 10) +
       theme(plot.title = element_text(face = "bold"),
             panel.grid.major.x = element_line(color = "grey92"))
  save_pdf_png(p, file.path(paths$viz_dir,
                            paste0(cmp, "_substate_composition_niu_vs_viral")),
               w = 9, h = max(4, 0.45 * length(unique(comp$substate)) + 2))
}

# ---------------------------------------------------------------------------
# Volcano of Autoimmune vs Viral pseudobulk DEGs (compartment-pooled "global")
# Filename: <cmp>_niu_vs_viral_volcano.pdf/png
# ---------------------------------------------------------------------------
viz_compartment_volcano <- function(cmp, paths, cfg) {
  dge_path <- file.path(paths$results_tables,
                        "DGE_pseudobulk_Autoimmune_vs_Viral.csv")
  if (!file.exists(dge_path)) return(invisible())
  dge <- read.csv(dge_path, stringsAsFactors = FALSE)
  if (nrow(dge) == 0 || !"padj" %in% colnames(dge)) return(invisible())
  glob <- if ("cluster" %in% colnames(dge) && any(dge$cluster == "global")) {
    dge |> dplyr::filter(cluster == "global")
  } else dge

  glob$direction <- factor(
    ifelse(glob$padj < 0.05 & glob$log2FoldChange >  1, "NIU_up",
    ifelse(glob$padj < 0.05 & glob$log2FoldChange < -1, "Viral_up", "ns")),
    levels = c("NIU_up", "Viral_up", "ns"))

  top_lab <- glob |>
    dplyr::filter(direction != "ns", !is.na(padj)) |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = 25)

  p <- ggplot(glob, aes(x = log2FoldChange, y = -log10(pmax(padj, 1e-300)),
                        color = direction)) +
       geom_point(alpha = 0.5, size = 1.2) +
       geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
       geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
       scale_color_manual(values = ETIOLOGY_DIRECTION_COLORS) +
       ggrepel::geom_text_repel(data = top_lab, aes(label = gene),
                                size = 3, max.overlaps = 30, color = "black") +
       labs(title = paste0(tools::toTitleCase(cmp),
                           " DEGs: NIU vs Viral (compartment-global)"),
            x = expression(log[2]~FC~"(NIU / Viral)"),
            y = expression(-log[10]~padj),
            color = "Direction") +
       theme_classic() +
       theme(plot.title = element_text(face = "bold"))
  save_pdf_png(p, file.path(paths$viz_dir,
                            paste0(cmp, "_niu_vs_viral_volcano")),
               w = 8, h = 7)
}

# ---------------------------------------------------------------------------
# Per-substate volcano grid (one mini-volcano per Leiden substate)
# Surfaces substate-unique programs that the global volcano hides.
# Filename: <cmp>_niu_vs_viral_volcano_per_substate.pdf/png
# ---------------------------------------------------------------------------
viz_compartment_volcano_per_substate <- function(cmp, paths, cfg) {
  dge_path <- file.path(paths$results_tables,
                        "DGE_pseudobulk_Autoimmune_vs_Viral.csv")
  if (!file.exists(dge_path)) return(invisible())
  dge <- read.csv(dge_path, stringsAsFactors = FALSE)
  if (!"cluster" %in% colnames(dge)) return(invisible())
  per <- dge |>
    dplyr::filter(cluster != "global", !is.na(padj), !is.na(log2FoldChange))
  if (nrow(per) == 0) return(invisible())
  per$substate <- substate_labels(cfg, cmp, per$cluster)
  per$direction <- factor(
    ifelse(per$padj < 0.05 & per$log2FoldChange >  1, "NIU_up",
    ifelse(per$padj < 0.05 & per$log2FoldChange < -1, "Viral_up", "ns")),
    levels = c("NIU_up", "Viral_up", "ns"))

  top_lab <- per |>
    dplyr::filter(direction != "ns") |>
    dplyr::group_by(substate) |>
    dplyr::slice_min(padj, n = 8, with_ties = FALSE) |>
    dplyr::ungroup()

  p <- ggplot(per, aes(x = log2FoldChange, y = -log10(pmax(padj, 1e-300)),
                       color = direction)) +
       geom_point(alpha = 0.55, size = 1.0) +
       geom_vline(xintercept = c(-1, 1), linetype = "dashed",
                  color = "grey60", linewidth = 0.3) +
       geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                  color = "grey60", linewidth = 0.3) +
       facet_wrap(~ substate, scales = "free", ncol = 3) +
       scale_color_manual(values = ETIOLOGY_DIRECTION_COLORS) +
       ggrepel::geom_text_repel(data = top_lab, aes(label = gene),
                                size = 2.5, max.overlaps = 15, color = "black") +
       labs(title = paste0(tools::toTitleCase(cmp),
                           " DEGs: NIU vs Viral by substate"),
            x = expression(log[2]~FC~"(NIU / Viral)"),
            y = expression(-log[10]~padj),
            color = "Direction") +
       theme_classic(base_size = 9) +
       theme(strip.background = element_blank(),
             strip.text = element_text(face = "bold"))
  n_sub <- length(unique(per$substate))
  save_pdf_png(p, file.path(paths$viz_dir,
                            paste0(cmp, "_niu_vs_viral_volcano_per_substate")),
               w = 12, h = max(6, ceiling(n_sub / 3) * 3.5))
}

# ---------------------------------------------------------------------------
# NIU sub-contrast heatmap: each NIU etiology vs other-NIU-pooled. Cell is
# signed -log10 q from the limma fit so HLA_B27-specific (vs other NIU)
# pathways pop out separately from VKH-specific etc.
# Filename: <cmp>_niu_subcontrast_heatmap.pdf/png
# ---------------------------------------------------------------------------
viz_compartment_niu_subcontrast_heatmap <- function(cmp, paths, cfg) {
  tab_path <- file.path(paths$results_tables,
                        paste0("Pathway_NIU_subcontrast_", cmp, ".csv"))
  if (!file.exists(tab_path)) {
    log_message("  NIU sub-contrast heatmap: ", basename(tab_path),
                " not found.")
    return(invisible())
  }
  d <- read.csv(tab_path, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(invisible())

  d$signed_nlp <- sign(d$logFC) * -log10(pmax(d$FDR, 1e-50))

  # Keep only pathways that are sig (FDR<0.1) in at least one NIU etiology
  sig_paths <- d |>
    dplyr::filter(FDR < 0.1) |>
    dplyr::pull(pathway) |>
    unique()
  if (length(sig_paths) == 0) {
    log_message("  NIU sub-contrast heatmap: no FDR<0.1 hits to plot.")
    return(invisible())
  }
  d <- d |> dplyr::filter(pathway %in% sig_paths)

  # Cluster rows
  wide <- d |>
    dplyr::select(pathway, etiology, signed_nlp) |>
    tidyr::pivot_wider(names_from = etiology, values_from = signed_nlp,
                       values_fill = 0) |>
    as.data.frame()
  rownames(wide) <- wide$pathway; wide$pathway <- NULL
  if (nrow(wide) >= 3) {
    ord <- tryCatch(stats::hclust(stats::dist(as.matrix(wide)))$order,
                    error = function(e) seq_len(nrow(wide)))
    d$pathway <- factor(d$pathway, levels = rownames(wide)[ord])
  }

  p <- ggplot(d, aes(x = etiology, y = pathway, fill = signed_nlp)) +
       geom_tile(color = "grey90", linewidth = 0.2) +
       viridis::scale_fill_viridis(option = "viridis",
                                   name = "sign(logFC) x -log10(FDR)") +
       labs(title = paste0(tools::toTitleCase(cmp),
                           " NIU sub-contrast: each etiology vs other-NIU pooled"),
            subtitle = "Positive = up in this NIU etiology vs other NIU",
            x = "NIU etiology", y = NULL) +
       theme_minimal(base_size = 10) +
       theme(axis.text.x = element_text(angle = 30, hjust = 1),
             panel.grid = element_blank())
  save_pdf_png(p, file.path(paths$viz_dir,
                            paste0(cmp, "_niu_subcontrast_heatmap")),
               w = 10, h = max(6, 0.25 * length(unique(d$pathway))))
}

# ---------------------------------------------------------------------------
# Module score helper: AddModuleScore + return augmented obj + score column
# ---------------------------------------------------------------------------
add_module_score_safe <- function(obj, features, name, ctrl = 25) {
  features <- list(intersect(features, rownames(obj)))
  if (length(features[[1]]) < 3) return(list(obj = obj, col = NULL))
  obj <- AddModuleScore(obj, features = features, name = name, ctrl = ctrl)
  list(obj = obj, col = paste0(name, "1"))
}

# ---------------------------------------------------------------------------
# Module-score violin: NIU vs Viral, faceted by substate.
# `score_specs` is a named list: list(Type17 = c("RORC", ...), ISG = c(...))
# Filename: <cmp>_<key>_score_niu_vs_viral.pdf/png
# ---------------------------------------------------------------------------
viz_compartment_module_score_niu_vs_viral <- function(obj, cmp, paths, cfg,
                                                      score_specs) {
  meta_orig <- obj@meta.data
  for (key in names(score_specs)) {
    res <- add_module_score_safe(obj, score_specs[[key]],
                                 paste0(key, "_mod_"))
    obj <- res$obj
    if (is.null(res$col)) next
    obj$substate <- substate_labels(cfg, cmp, obj$knn.leiden.cluster)
    df <- obj@meta.data |>
      dplyr::filter(Phenotype_2 %in% c("NIU", "Viral")) |>
      dplyr::select(substate, Phenotype_2, score = !!res$col)
    if (nrow(df) == 0) next
    p <- ggplot(df, aes(x = Phenotype_2, y = score, fill = Phenotype_2)) +
         geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
         geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                      color = "grey20", alpha = 0.7) +
         facet_wrap(~ substate, scales = "free_y") +
         scale_fill_manual(values = ETIOLOGY_GROUP_COLORS) +
         labs(title = paste0(tools::toTitleCase(cmp), " ", key,
                             " module score: NIU vs Viral"),
              x = "Phenotype", y = paste0(key, " module score")) +
         theme_classic(base_size = 10) +
         theme(legend.position = "none")
    save_pdf_png(p, file.path(paths$viz_dir,
                              paste0(cmp, "_", tolower(key),
                                     "_score_niu_vs_viral")),
                 w = 10, h = 8)
  }
}

# ===========================================================================
# DA-focused panel: Milo DA per-cluster call + figure-D/E analogues + GSEA
# bubble heatmap + per-cluster bar + per-DA-cluster UCell hex + etiology
# resolution within DA clusters (composition, subject-fraction, pathway
# violins/heatmap, Fisher enrichment).
#
# Inputs are produced upstream by R/40_milo.R (per-contrast DA + per-grouping
# summary CSVs).
# ===========================================================================

# Map compartment short-name -> milo contrast name written by R/40_milo.R.
.milo_contrast_name <- function(cmp) paste0("etiology_", cmp, "_eye")

# Cluster-level DA call from the upstream per-grouping summary table.
# Returns a tibble with one row per Leiden cluster:
#   cluster, substate_label, n_neighborhoods, n_sig_DA, frac_sig,
#   median_logFC, direction, da_call (logical).
# Persists the result to outputs/tables/eye/<cmp>/milo_da_cluster_calls.csv.
# Returns NULL invisibly when the summary table is missing or empty.
compartment_milo_da_clusters <- function(cmp, paths, cfg) {
  contrast <- .milo_contrast_name(cmp)
  summary_path <- file.path(paths$results_tables,
                            paste0("milo_summary_", contrast,
                                   "_knn_leiden.csv"))
  if (!file.exists(summary_path)) {
    log_message("  milo DA: ", basename(summary_path), " not found.")
    return(invisible(NULL))
  }
  s <- read.csv(summary_path, stringsAsFactors = FALSE)
  if (nrow(s) == 0) return(invisible(NULL))
  min_sig <- as.integer(cfg$milo$cluster_da_min_sig_nhoods %||% 2L)
  out <- s |>
    dplyr::transmute(cluster = as.character(.data$knn.leiden.cluster),
                     n_neighborhoods = .data$n_neighborhoods,
                     n_sig_DA        = .data$n_sig_DA,
                     n_enriched      = .data$n_enriched,
                     n_depleted      = .data$n_depleted,
                     median_logFC    = .data$median_logFC,
                     mean_logFC      = .data$mean_logFC) |>
    dplyr::mutate(frac_sig  = ifelse(.data$n_neighborhoods > 0,
                                     .data$n_sig_DA / .data$n_neighborhoods,
                                     0),
                  direction = ifelse(.data$n_enriched >= .data$n_depleted,
                                     "NIU_up", "Viral_up"),
                  da_call   = .data$n_sig_DA >= min_sig,
                  substate_label = substate_labels(cfg, cmp, .data$cluster))
  ensure_dir(paths$results_tables)
  write.csv(out,
            file.path(paths$results_tables, "milo_da_cluster_calls.csv"),
            row.names = FALSE)
  log_message("  milo DA: ", sum(out$da_call), "/", nrow(out),
              " clusters called DA (min_sig=", min_sig, ").")
  invisible(out)
}

# Figure-D analogue: per-cluster Milo logFC boxplot with significance marks.
# Reads Milo_<contrast>_DA.csv directly, filters to purity >= 0.7.
viz_compartment_milo_da_box <- function(cmp, paths, cfg) {
  contrast <- .milo_contrast_name(cmp)
  da_path <- file.path(paths$results_tables,
                       paste0("Milo_", contrast, "_DA.csv"))
  if (!file.exists(da_path)) {
    log_message("  milo DA box: ", basename(da_path), " not found.")
    return(invisible())
  }
  res <- read.csv(da_path, stringsAsFactors = FALSE)
  if (!"knn.leiden.cluster" %in% colnames(res) ||
      !"logFC" %in% colnames(res)) return(invisible())
  purity_col <- "knn.leiden.cluster_fraction"
  if (purity_col %in% colnames(res)) {
    res <- res |> dplyr::filter(.data[[purity_col]] >= 0.7,
                                !is.na(.data$knn.leiden.cluster))
  }
  if (nrow(res) == 0) return(invisible())

  da <- compartment_milo_da_clusters(cmp, paths, cfg)
  da_clusters <- if (!is.null(da)) da$cluster[da$da_call] else character(0)

  res$cluster <- as.character(res$knn.leiden.cluster)
  res$substate_label <- substate_labels(cfg, cmp, res$cluster)

  med <- res |>
    dplyr::group_by(.data$cluster, .data$substate_label) |>
    dplyr::summarise(med_lfc = stats::median(.data$logFC, na.rm = TRUE),
                     min_sfdr = suppressWarnings(min(.data$SpatialFDR, na.rm = TRUE)),
                     .groups = "drop") |>
    dplyr::arrange(.data$med_lfc)
  med$is_da <- med$cluster %in% da_clusters
  med$label_marked <- ifelse(med$is_da,
                             paste0(med$substate_label, " **"),
                             med$substate_label)
  res$substate_label <- factor(res$substate_label,
                               levels = med$substate_label)
  med$substate_label <- factor(med$substate_label,
                               levels = med$substate_label)

  res$direction <- factor(ifelse(res$cluster %in% med$cluster[med$med_lfc > 0],
                                 "NIU_up", "Viral_up"),
                          levels = c("NIU_up", "Viral_up"))

  p <- ggplot(res, aes(y = .data$substate_label, x = .data$logFC,
                       fill = .data$direction)) +
       geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
       geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65) +
       geom_jitter(height = 0.15, alpha = 0.35, size = 0.7,
                   color = "grey20") +
       scale_fill_manual(values = ETIOLOGY_DIRECTION_COLORS[c("NIU_up",
                                                              "Viral_up")],
                         labels = c(NIU_up = "Increased in NIU",
                                    Viral_up = "Increased in Viral")) +
       scale_y_discrete(labels = setNames(med$label_marked,
                                          med$substate_label)) +
       labs(title = paste0(tools::toTitleCase(cmp),
                           " Milo DA: logFC per substate"),
            subtitle = "** = cluster has ≥ N nhoods with SpatialFDR<0.1 (cfg$milo$cluster_da_min_sig_nhoods)",
            x = expression(log[2]~FC~"(NIU / Viral)"),
            y = NULL, fill = "Direction") +
       theme_classic(base_size = 10) +
       theme(plot.title = element_text(face = "bold"))
  save_pdf_png(p, file.path(paths$viz_dir,
                            paste0(cmp, "_milo_da_logfc_box")),
               w = 9, h = max(4, 0.5 * length(unique(med$substate_label)) + 2))
}

# ---------------------------------------------------------------------------
# Etiology resolution within DA clusters (7a-7e)
# ---------------------------------------------------------------------------

# Curated immunological keyword list used as a fallback when a pathway has no
# upstream axis_tag (R/32_escape.R::.tag_axis only matches a compartment-specific
# regex; pathways outside that regex are tagged via the keywords below).
# Order matters: first match wins, so put more specific terms before generic.
.FUNCTIONAL_KEYWORDS <- c(
  "INTERFERON_ALPHA", "INTERFERON_BETA", "INTERFERON_GAMMA", "INTERFERON",
  "ANTIGEN_PRESENTATION", "MHC_CLASS_I", "MHC_CLASS_II", "MHC",
  "COMPLEMENT", "COAGULATION",
  "TNFA_SIGNALING", "TNF",
  "IL17", "IL23", "IL10", "IL6", "IL2_STAT5", "IL4",
  "TH17", "TH1", "TH2", "TFH", "TREG",
  "INFLAMMATORY_RESPONSE", "INFLAMMATION", "NEUROINFLAMMATION",
  "CLASS_SWITCH_RECOMBINATION", "SOMATIC_HYPERMUTATION",
  "GERMINAL_CENTER", "MEMORY_B", "PLASMA",
  "TISSUE_RESIDENT", "EXHAUSTION", "CYTOTOXIC", "EFFECTOR", "STEMNESS",
  "APOPTOSIS", "AUTOPHAGY", "HYPOXIA", "GLYCOLYSIS",
  "OXIDATIVE_PHOSPHORYLATION", "FATTY_ACID", "CHOLESTEROL",
  "TGF_BETA", "TGFB", "WNT", "NOTCH", "HEDGEHOG",
  "NF_KAPPA_B", "NFKB", "JAK_STAT", "MAPK",
  "TLR_SIGNALING", "RIG_I", "STING", "INFLAMMASOME", "NLRP3",
  "CELL_CYCLE", "MITOTIC", "DNA_REPAIR", "P53"
)

# Cross-cluster functional-gestalt heatmap. Works from per-cell UCell scores
# (escape.UCell_eye assay if present, else loaded from EscapeChunkManifest)
# rather than the NIU-vs-Viral GSEA logFC. Pathways shown = those that are
# differential *across the cell clusters* (max |z| across cluster means >=
# z_threshold; ANOVA optional via cfg$gsea$cross_cluster_anova=TRUE), so the
# heatmap captures cluster identity / cell-state biology rather than the
# etiology contrast. Tile values = z-score of mean UCell per pathway, so
# yellow = cluster most enriched for that pathway, dark purple = least
# enriched. Functional grouping (rows split by k=8) is unchanged in spirit
# but recomputed on this new matrix.
# Filename: <cmp>_functional_gestalt_full_heatmap.{pdf,png}
viz_compartment_functional_gestalt_full_heatmap <- function(obj, cmp, paths, cfg,
                                                            n_groups = NULL,
                                                            z_threshold = NULL) {
  if (is.null(n_groups))
    n_groups <- as.integer(cfg$gsea$functional_n_groups %||% 8L)
  n_groups <- max(4L, min(8L, n_groups))
  if (is.null(z_threshold))
    z_threshold <- as.numeric(cfg$gsea$cross_cluster_z_threshold %||% 1.5)

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    log_message("  cross-cluster gestalt: ComplexHeatmap/circlize required.")
    return(invisible())
  }

  # 1) Get the full per-cell UCell matrix from the carried-forward assay first,
  # then fall back to the chunk manifest.
  esc_assay <- if ("escape.UCell_eye" %in% Seurat::Assays(obj)) "escape.UCell_eye"
               else if ("escape.UCell"  %in% Seurat::Assays(obj)) "escape.UCell"
               else NA_character_
  esc_mat <- NULL
  if (!is.na(esc_assay)) {
    esc_mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = esc_assay, layer = "data"),
      error = function(e) NULL)
  }
  if (is.null(esc_mat)) {
    manifest_path <- file.path(get_target_paths(cfg, "all")$results_objects,
                               "EscapeChunkManifest.rds")
    if (file.exists(manifest_path) && exists(".load_escape_cells")) {
      manifest <- readRDS(manifest_path)
      esc_mat <- tryCatch(.load_escape_cells(manifest, colnames(obj),
                                             type = "unnorm"),
                          error = function(e) NULL)
    }
  }
  if (is.null(esc_mat) || nrow(esc_mat) == 0) {
    log_message("  cross-cluster gestalt: no escape scores available.")
    return(invisible())
  }

  # 2) Restrict pathways to immune-tagged (axis_regex or .FUNCTIONAL_KEYWORDS).
  all_paths <- rownames(esc_mat)
  resolve_func_tag <- function(pw) {
    nm <- toupper(pw)
    for (kw in .FUNCTIONAL_KEYWORDS) {
      pat <- gsub("_", "[_-]?", kw, fixed = TRUE)
      if (grepl(pat, nm)) return(kw)
    }
    NA_character_
  }
  tag_lookup <- vapply(all_paths, resolve_func_tag, character(1))
  names(tag_lookup) <- all_paths
  immune_paths <- all_paths[!is.na(tag_lookup)]
  if (length(immune_paths) < n_groups + 1L) {
    log_message("  cross-cluster gestalt: only ", length(immune_paths),
                " immune-tagged pathways in escape assay; skipping.")
    return(invisible())
  }
  esc_mat <- esc_mat[immune_paths, , drop = FALSE]
  tag_lookup <- tag_lookup[immune_paths]

  # 3) Mean UCell per (pathway, cluster).
  meta <- obj@meta.data
  if (!"knn.leiden.cluster" %in% colnames(meta)) return(invisible())
  cluster_int <- as.integer(as.character(meta$knn.leiden.cluster))
  names(cluster_int) <- rownames(meta)
  shared <- intersect(colnames(esc_mat), names(cluster_int))
  shared <- shared[!is.na(cluster_int[shared])]
  if (length(shared) == 0) return(invisible())
  esc_mat <- esc_mat[, shared, drop = FALSE]
  cluster_int <- cluster_int[shared]

  uniq_clusters <- sort(unique(cluster_int))
  mean_mat <- matrix(NA_real_,
                     nrow = nrow(esc_mat), ncol = length(uniq_clusters),
                     dimnames = list(rownames(esc_mat),
                                     as.character(uniq_clusters)))
  for (cid in uniq_clusters) {
    cells <- shared[cluster_int == cid]
    if (length(cells) > 0) {
      mat_c <- esc_mat[, cells, drop = FALSE]
      mean_mat[, as.character(cid)] <-
        if (methods::is(mat_c, "sparseMatrix"))
          Matrix::rowMeans(mat_c, na.rm = TRUE)
        else
          rowMeans(as.matrix(mat_c), na.rm = TRUE)
    }
  }

  # 4) Differential filter: z-score per pathway across clusters, keep pathways
  # with at least one cluster whose mean is >= z_threshold SD from the
  # pathway's average.
  z_mat <- t(scale(t(mean_mat)))
  z_max <- suppressWarnings(apply(abs(z_mat), 1, max, na.rm = TRUE))
  diff_paths <- rownames(mean_mat)[is.finite(z_max) & z_max >= z_threshold]
  if (length(diff_paths) < n_groups + 1L) {
    log_message("  cross-cluster gestalt: only ", length(diff_paths),
                " differential pathways (max|z|>=", z_threshold,
                "); skipping.")
    return(invisible())
  }
  log_message("  cross-cluster gestalt: ", length(diff_paths),
              " of ", length(immune_paths),
              " immune-tagged pathways differential (max|z|>=",
              z_threshold, ").")
  z_disp <- z_mat[diff_paths, , drop = FALSE]
  mean_disp <- mean_mat[diff_paths, , drop = FALSE]

  # 5) hclust both axes on the z-score matrix; cut pathway dendrogram.
  row_hc <- stats::hclust(stats::dist(z_disp, method = "euclidean"),
                          method = "ward.D2")
  col_hc <- stats::hclust(stats::dist(t(z_disp), method = "euclidean"),
                          method = "ward.D2")
  pathway_group <- stats::cutree(row_hc, k = n_groups)

  group_labels <- character(n_groups)
  for (gid in seq_len(n_groups)) {
    members <- names(pathway_group)[pathway_group == gid]
    counts  <- sort(table(tag_lookup[members]), decreasing = TRUE)
    top_one <- names(counts)[1]
    if (length(counts) >= 2 && counts[2] >= 0.5 * counts[1]) {
      lbl <- paste(names(counts)[1:2], collapse = "/")
    } else {
      lbl <- top_one
    }
    group_labels[gid] <- paste0(gid, ": ", lbl,
                                " (n=", length(members), ")")
  }
  group_pos <- vapply(seq_len(n_groups), function(gid)
    mean(which(pathway_group[row_hc$order] == gid)), numeric(1))
  group_levels <- group_labels[order(group_pos)]

  # 6) Persist pathway -> functional-group + raw cluster means for inspection.
  ensure_dir(paths$results_tables)
  write.csv(
    data.frame(
      pathway                = names(pathway_group),
      functional_group_id    = unname(pathway_group),
      functional_group_label = group_labels[pathway_group],
      functional_tag         = unname(tag_lookup[names(pathway_group)]),
      max_abs_z              = unname(z_max[names(pathway_group)]),
      stringsAsFactors       = FALSE),
    file.path(paths$results_tables,
              "cross_cluster_gestalt_pathway_assignments.csv"),
    row.names = FALSE)
  write.csv(
    cbind(pathway = rownames(mean_disp),
          as.data.frame(mean_disp)),
    file.path(paths$results_tables,
              "cross_cluster_gestalt_mean_ucell_per_cluster.csv"),
    row.names = FALSE)
  log_message("  cross-cluster gestalt core: ", length(diff_paths),
              " pathways -> ", n_groups, " groups: ",
              paste(group_labels, collapse = "; "))

  # 7) Render heatmap. Fill = z-score of mean UCell per pathway across
  # clusters (yellow = cluster most enriched for that pathway, dark = least),
  # so rows are directly comparable.
  col_labels <- substate_labels(cfg, cmp, colnames(z_disp))
  row_group  <- factor(group_labels[pathway_group[rownames(z_disp)]],
                       levels = group_levels)
  fill_lim <- max(abs(z_disp), na.rm = TRUE)
  if (!is.finite(fill_lim) || fill_lim == 0) fill_lim <- 1
  col_fun <- circlize::colorRamp2(
    seq(-fill_lim, fill_lim, length.out = 11),
    viridis::viridis(11, option = "viridis"))

  ht <- ComplexHeatmap::Heatmap(
    z_disp,
    name = "Z-score\nof mean UCell\n(per pathway,\nacross clusters)",
    col = col_fun,
    row_split = row_group,
    cluster_rows = TRUE,
    cluster_columns = col_hc,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    column_labels = col_labels,
    column_names_rot = 35,
    row_title_side = "left",
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    column_title = paste0(tools::toTitleCase(cmp),
                          " cross-cluster functional gestalt (k=",
                          n_groups, ", max|z|>=", z_threshold, ")"),
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(legend_direction = "horizontal",
                                legend_width = grid::unit(4, "cm")),
    border = TRUE,
    use_raster = TRUE,
    raster_quality = 4)

  base <- file.path(paths$viz_dir,
                    paste0(cmp, "_functional_gestalt_full_heatmap"))
  h_in <- max(8, n_groups * 1.5 + 4)
  w_in <- max(9, ncol(z_disp) * 0.6 + 5)
  grDevices::pdf(paste0(base, ".pdf"), width = w_in, height = h_in)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE)
  grDevices::dev.off()
  grDevices::png(paste0(base, ".png"), width = w_in, height = h_in,
                 units = "in", res = 200)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE)
  grDevices::dev.off()
  log_message("  Saved: ", basename(base), ".{pdf,png}")
}

# ---------------------------------------------------------------------------
# Per-substate DGE / pseudobulk helpers (used by panel D variants in
# 85_viz_myeloid.R; B-cell / T-cell figures can call the same helpers later).
#
# read_per_substate_dge(): filter the long-format pseudobulk DGE CSV that
# run_dge() writes (one row per gene x cluster) into a tibble with an explicit
# `substate` column. Caller decides padj / lfc thresholds.
#
# build_per_substate_pseudobulks(): construct a named list of
# SummarizedExperiment objects keyed by cluster id, with one column per
# Subject_Timepoint (or orig.ident) within each cluster. Counts come from
# Seurat::AggregateExpression on the RNA assay's counts slot; coldata carries
# Phenotype_2 and the sample id. Builds on the fly so D-v2 / D-v3 do not
# require any new on-disk artifact from the DGE step.
# ---------------------------------------------------------------------------
read_per_substate_dge <- function(paths,
                                  contrast_name = "Autoimmune_vs_Viral") {
  f <- file.path(paths$results_tables,
                 paste0("DGE_pseudobulk_", contrast_name, ".csv"))
  if (!file.exists(f)) return(NULL)
  df <- utils::read.csv(f, stringsAsFactors = FALSE)
  if (!"cluster" %in% colnames(df)) return(NULL)
  # The DGE table writes the contrast direction such that positive log2FC
  # corresponds to the FIRST contrast group (NIU == Autoimmune for this
  # contrast). Keep the raw column and let the caller relabel.
  df$substate <- as.character(df$cluster)
  tibble::as_tibble(df)
}

build_per_substate_pseudobulks <- function(obj,
                                           cluster_col = "knn.leiden.cluster",
                                           group_col   = "Phenotype_2",
                                           groups      = c("NIU", "Viral"),
                                           min_cells_per_pb = 10) {
  if (!cluster_col %in% colnames(obj[[]])) return(list())
  meta <- obj[[]]
  sample_col <- if ("Subject_Timepoint" %in% colnames(meta))
                  "Subject_Timepoint" else "orig.ident"
  pb_sample_vec  <- as.character(meta[[sample_col]])
  pb_cluster_vec <- as.character(meta[[cluster_col]])
  pb_group_vec   <- as.character(meta[[group_col]])
  keep <- pb_group_vec %in% groups & !is.na(pb_sample_vec) & !is.na(pb_cluster_vec)
  if (sum(keep) < 50) return(list())

  # Write the helper columns onto the SEURAT object (not just the local meta
  # copy) so AggregateExpression can resolve them. Subset afterwards.
  obj[[".pb_cluster"]] <- pb_cluster_vec
  obj[[".pb_sample"]]  <- pb_sample_vec
  obj[[".pb_group"]]   <- pb_group_vec
  obj_sub <- subset(obj, cells = rownames(meta)[keep])

  agg <- tryCatch(
    Seurat::AggregateExpression(obj_sub, assays = "RNA", slot = "counts",
                                group.by = c(".pb_cluster", ".pb_sample",
                                             ".pb_group"),
                                return.seurat = FALSE),
    error = function(e) NULL)
  if (is.null(agg)) return(list())
  mat <- agg$RNA
  if (is.null(mat) || ncol(mat) == 0) return(list())

  # Seurat's AggregateExpression joins the three group.by fields with `_` AND
  # internally substitutes any underscores in the original values with `-`.
  # So a cell with cluster=0, sample="CMV_CRN2_VIT", group="NIU" lands in
  # column "g0_CMV-CRN2-VIT_NIU". Build the lookup by joining the *converted*
  # values ourselves and matching on the joined string, which is robust to
  # any future change in Seurat's group separator.
  obj_meta <- obj_sub[[]]
  cv <- function(x) gsub("_", "-", as.character(x), fixed = TRUE)
  cell_cluster <- cv(obj_meta$.pb_cluster)
  cell_sample  <- cv(obj_meta$.pb_sample)
  cell_group   <- cv(obj_meta$.pb_group)
  cell_col     <- paste(cell_cluster, cell_sample, cell_group, sep = "_")
  # AggregateExpression prepends "g" to integer cluster ids. Match either form.
  col_lookup_a <- paste(cell_cluster, cell_sample, cell_group, sep = "_")
  col_lookup_b <- paste0("g", cell_cluster, "_", cell_sample, "_", cell_group)

  pb_meta <- data.frame(col = colnames(mat),
                        cluster = NA_character_,
                        sample  = NA_character_,
                        group   = NA_character_,
                        n_cells = 0L,
                        stringsAsFactors = FALSE)

  for (i in seq_len(nrow(pb_meta))) {
    cn <- pb_meta$col[i]
    idx <- which(col_lookup_a == cn | col_lookup_b == cn)
    if (length(idx) == 0) next
    pb_meta$cluster[i] <- obj_meta$.pb_cluster[idx[1]]
    pb_meta$sample[i]  <- obj_meta$.pb_sample[idx[1]]
    pb_meta$group[i]   <- obj_meta$.pb_group[idx[1]]
    pb_meta$n_cells[i] <- length(idx)
  }

  keep_cols <- !is.na(pb_meta$cluster) & pb_meta$n_cells >= min_cells_per_pb
  if (sum(keep_cols) == 0) return(list())
  mat     <- mat[, keep_cols, drop = FALSE]
  pb_meta <- pb_meta[keep_cols, , drop = FALSE]

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    # Return matrices + metadata rather than failing — caller can still use
    # rowMeans / cor.
    out <- split(seq_len(nrow(pb_meta)), pb_meta$cluster)
    return(lapply(out, function(idx) {
      list(counts = mat[, idx, drop = FALSE],
           coldata = pb_meta[idx, , drop = FALSE])
    }))
  }

  cluster_ids <- sort(unique(pb_meta$cluster))
  lapply(stats::setNames(cluster_ids, cluster_ids), function(ck) {
    idx <- which(pb_meta$cluster == ck)
    SummarizedExperiment::SummarizedExperiment(
      assays  = list(counts = mat[, idx, drop = FALSE]),
      colData = S4Vectors::DataFrame(pb_meta[idx, , drop = FALSE]))
  })
}
