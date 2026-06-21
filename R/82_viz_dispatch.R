# R/82_viz_dispatch.R
# Master visualization dispatcher + shared viz_* helpers used by both Figure 1
# (full atlas) and Figure 2 (eye sub-atlas). Per-figure entry points live in:
#   R/83_viz_full.R     — run_visualizations_full + F1-only viz blocks
#   R/84_viz_eye.R      — run_visualizations_eye
#   R/85_viz_myeloid.R  — run_visualizations_myeloid (Figure 3)
#   R/86_viz_bcell.R    — run_visualizations_bcell   (Figure 4)
#   R/88_viz_tcell.R    — run_visualizations_tcell   (Figure 5)
suppressPackageStartupMessages({
  library(Seurat)
  library(scplotter)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(forcats)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Helper: pick merged clusters when available ---
.cluster_col <- function(obj) {
  if ("merged.celltype.cluster" %in% colnames(obj[[]]))
    "merged.celltype.cluster" else "knn.leiden.cluster"
}

# --- Helper: safe plot wrapper ---
safe_plot <- function(expr, file, width = 10, height = 7) {
  tryCatch({
    p <- eval(expr, envir = parent.frame())
    if (!is.null(p)) {
      ensure_dir(dirname(file))
      # ComplexHeatmap objects (HeatmapList, Heatmap) need pdf/draw/dev.off
      if (inherits(p, c("HeatmapList", "Heatmap", "AdditiveUnit"))) {
        ext <- tolower(tools::file_ext(file))
        if (ext == "pdf") {
          grDevices::pdf(file, width = width, height = height)
        } else {
          grDevices::pdf(file, width = width, height = height, units = "in")
        }
        tryCatch(ComplexHeatmap::draw(p), finally = grDevices::dev.off())
      } else {
        ggsave(file, p, width = width, height = height)
      }
      log_message("  Saved: ", basename(file))
    }
  }, error = function(e) {
    log_message("  WARN plot failed (", basename(file), "): ", conditionMessage(e))
  })
}

# ============================================================================
# MASTER ENTRY POINT — dispatches to per-figure file
# ============================================================================
# run_visualizations is the public entry point called from run_pipeline.R. It
# routes to the per-figure entry points defined in 83_viz_full.R,
# 84_viz_eye.R, 85_viz_myeloid.R, 86_viz_bcell.R, 88_viz_tcell.R.
run_visualizations <- function(cfg, target = c("all", "eye", "myeloid", "bcell", "tcell",
                                                "cross_compartment")) {
  target <- match.arg(target)
  switch(target,
    all     = run_visualizations_full(cfg),
    eye     = run_visualizations_eye(cfg),
    myeloid = run_visualizations_myeloid(cfg),
    bcell   = run_visualizations_bcell(cfg),
    tcell   = run_visualizations_tcell(cfg),
    cross_compartment = run_visualizations_cross_compartment(cfg))
}

# ============================================================================
# 1. QC SUMMARY
# ============================================================================
viz_qc_summary <- function(obj, cfg, paths = cfg$paths) {
  log_message("Visualizing: QC summary")
  vdir <- viz_subdir(paths, "qc")

  # Cell counts per sample, colored by tissue
  safe_plot(quote(
    CellStatPlot(obj, ident = "orig.ident", group_by = "Tissue_1",
                 plot_type = "bar", frac = "none",
                 ylab = "Number of Cells",
                 x_text_angle = 60)
  ), file.path(vdir, "qc_cell_counts_by_sample.pdf"), width = 12, height = 6)

  # QC metrics violin by sample, split by tissue
  for (feat in c("nFeature_RNA", "nCount_RNA", "mito.perc")) {
    if (!feat %in% colnames(obj[[]])) next
    safe_plot(bquote(
      FeatureStatPlot(obj, features = .(feat), ident = "orig.ident",
                      split_by = "Tissue_1", plot_type = "violin",
                      x_text_angle = 60)
    ), file.path(vdir, paste0("qc_violin_", feat, ".pdf")), width = 14, height = 5)
  }

  # QC metrics violin by tissue (aggregated)
  qc_feats <- intersect(c("nFeature_RNA", "nCount_RNA", "mito.perc"), colnames(obj[[]]))
  if (length(qc_feats) > 0) {
    safe_plot(bquote(
      FeatureStatPlot(obj, features = .(qc_feats), ident = "Tissue_1",
                      plot_type = "violin")
    ), file.path(vdir, "qc_violin_by_tissue.pdf"), width = 10, height = 5)
  }

  # Doublet rate and cell filtering from QC summary table
  qc_summary_path <- file.path(paths$qc_dir, "qc_summary.csv")
  if (file.exists(qc_summary_path)) {
    qc_summary <- read.csv(qc_summary_path, check.names = FALSE)

    safe_plot(bquote({
      ggplot(.(qc_summary), aes(x = reorder(sample, -pct_doublets), y = pct_doublets,
                                 fill = Tissue_1)) +
        geom_col(alpha = 0.8) +
        theme_minimal() +
        labs(title = "Doublet Detection Rate per Sample",
             x = "Sample", y = "% Doublets") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
    }), file.path(vdir, "qc_doublet_rate_per_sample.pdf"), width = 12, height = 5)

    safe_plot(bquote({
      plot_df <- .(qc_summary) %>%
        tidyr::pivot_longer(cols = c(cells_pre_filter, cells_post_doublet),
                            names_to = "stage", values_to = "n_cells") %>%
        mutate(stage = factor(stage,
                              levels = c("cells_pre_filter", "cells_post_doublet"),
                              labels = c("Pre-filter", "Post-QC")))
      ggplot(plot_df, aes(x = reorder(sample, -n_cells), y = n_cells, fill = stage)) +
        geom_col(position = "dodge", alpha = 0.8) +
        theme_minimal() +
        labs(title = "Cell Counts: Pre- vs Post-QC Filtering",
             x = "Sample", y = "Number of Cells", fill = "Stage") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
    }), file.path(vdir, "qc_cells_pre_post_filter.pdf"), width = 14, height = 5)
  }
}

# ============================================================================
# 2. INTEGRATION & CLUSTERING
# ============================================================================
viz_integration <- function(obj, cfg, paths = cfg$paths) {
  log_message("Visualizing: Integration & clustering")
  vdir <- viz_subdir(paths, "integration")
  ccol <- .cluster_col(obj)
  umap_sz <- cfg$visualization$umap_size %||% 8
  stripped_bg <- cfg$visualization$stripped_bg %||% "black"

  # UMAP by cluster (dual: labeled + stripped)
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = ccol,
                reduction = "UMAP",
                label = TRUE,
                label_insitu = TRUE,
                show_stat = FALSE,
                raster = FALSE)
  ), file.path(vdir, "umap_cluster"), width = umap_sz, height = umap_sz,
  stripped_bg = stripped_bg)

  # UMAP by sample
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = "orig.ident", reduction = "UMAP",
                show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_sample"), width = umap_sz, height = umap_sz,
  stripped_bg = stripped_bg)

  # UMAP by tissue
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = "Tissue_1", reduction = "UMAP",
                show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_tissue"), width = umap_sz, height = umap_sz,
  stripped_bg = stripped_bg)

  # UMAP by etiology
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = "Etiology", reduction = "UMAP",
                show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_etiology"), width = umap_sz, height = umap_sz,
  stripped_bg = stripped_bg)

  # UMAP by phenotype
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = "Phenotype_2", reduction = "UMAP",
                show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_phenotype2"), width = umap_sz, height = umap_sz,
  stripped_bg = stripped_bg)

  # UMAP split by tissue (wider for facets, not square)
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = ccol, split_by = "Tissue_1",
                reduction = "UMAP", label = TRUE, label_insitu = TRUE,
                show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_cluster_split_tissue"), width = 16, height = 7,
  stripped_bg = stripped_bg)

  # UMAP split by phenotype
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = ccol, split_by = "Phenotype_2",
                reduction = "UMAP", label = TRUE,
                show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_cluster_split_phenotype"), width = 20, height = 7,
  stripped_bg = stripped_bg)

  # UMAP density by tissue
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = "Tissue_1", reduction = "UMAP",
                add_density = TRUE, show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_density_tissue"), width = umap_sz, height = umap_sz,
  stripped_bg = stripped_bg)

  # UMAP with cluster pie stats by tissue
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = ccol, reduction = "UMAP",
                stat_by = "Tissue_1", stat_plot_type = "pie",
                label = TRUE, label_insitu = TRUE,
                show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_cluster_pie_tissue"), width = 10, height = 10,
  stripped_bg = stripped_bg)

  # Cluster composition: fraction by tissue
  safe_plot(quote(
    CellStatPlot(obj, ident = ccol, group_by = "Tissue_1",
                 frac = "group", plot_type = "bar")
  ), file.path(vdir, "cluster_composition_tissue.pdf"), width = 10, height = 6)

  # Cluster composition: fraction by phenotype
  safe_plot(quote(
    CellStatPlot(obj,
                 ident = ccol,
                 group_by = "Phenotype_2",
                 frac = "ident", 
                 facet_by = "Tissue_1",
                 plot_type = "bar",
                 palette = "Set1",
                 position = "stack")
  ), file.path(vdir, "cluster_composition_phenotype.pdf"), width = 10, height = 6)

  # Sankey: cluster -> tissue -> phenotype
  safe_plot(quote(
    CellStatPlot(obj, plot_type = "sankey",
                 group_by = c(ccol, "Tissue_1", "Phenotype_2"))
  ), file.path(vdir, "sankey_cluster_tissue_phenotype.pdf"), width = 12, height = 10)

  # Cluster proportion heatmap
  safe_plot(quote(
    CellStatPlot(obj,
                 ident = ccol,
                 group_by = "Tissue_1",
                 plot_type = "heatmap",
                 frac = "group",
                 palette = "viridis",
                 columns_split_by = "Subject",
                 swap = TRUE) 
  ), file.path(vdir, "cluster_heatmap_tissue.pdf"), width = 8, height = 8)

  # PCA elbow plot (variance from fastMNN)
  if ("fastMNN" %in% names(obj@reductions)) {
    emb <- obj@reductions$fastMNN@cell.embeddings
    var_per_dim <- apply(emb, 2, var)
    pct_var <- var_per_dim / sum(var_per_dim) * 100

    safe_plot(bquote({
      df <- data.frame(PC = seq_along(.(pct_var)), pct_var = .(pct_var))
      ggplot(df, aes(x = PC, y = pct_var)) +
        geom_point(size = 2) +
        geom_line() +
        theme_minimal() +
        labs(title = "fastMNN: Variance Explained per Dimension",
             x = "Dimension", y = "% Variance") +
        geom_vline(xintercept = 20, linetype = "dashed", color = "red", alpha = 0.5)
    }), file.path(vdir, "pca_elbow.pdf"), width = 8, height = 5)
  }
}

# ============================================================================
# 3. CELL TYPE ANNOTATION
# ============================================================================
viz_celltypes <- function(obj, cfg, paths = cfg$paths) {
  log_message("Visualizing: Cell type annotation")
  vdir <- viz_subdir(paths, "celltypes")

  if (!"celltype_broad" %in% colnames(obj[[]])) {
    log_message("  No celltype_broad column found. Skipping.")
    return(invisible(NULL))
  }

  umap_sz <- cfg$visualization$umap_size %||% 8
  stripped_bg <- cfg$visualization$stripped_bg %||% "black"

  # UMAP by celltype_broad
  dual_save_plot(quote(
    CellDimPlot(obj, group_by = "celltype_broad", reduction = "UMAP",
                label = TRUE, label_insitu = TRUE,
                show_stat = FALSE, raster = FALSE)
  ), file.path(vdir, "umap_celltype_broad"), width = umap_sz, height = umap_sz,
  stripped_bg = stripped_bg)

  # UMAP by celltype (fine)
  if ("celltype" %in% colnames(obj[[]])) {
    dual_save_plot(quote(
      CellDimPlot(obj, group_by = "celltype", reduction = "UMAP",
                  label = TRUE, label_insitu = TRUE,
                  show_stat = FALSE, raster = FALSE)
    ), file.path(vdir, "umap_celltype_fine"), width = umap_sz, height = umap_sz,
    stripped_bg = stripped_bg)
  }

  # Canonical marker dot plot
  CANONICAL_MARKERS <- list(
    "T cell"              = c("CD3D", "CD3E", "CD3G", "TRAC"),
    "CD4 T cell"          = c("CD4", "IL7R", "CCR7", "LEF1"),
    "CD8 T cell"          = c("CD8A", "CD8B", "GZMK", "GZMB"),
    "NK"                  = c("NKG7", "GNLY", "KLRD1", "NCAM1"),
    "B cell"              = c("MS4A1", "CD19", "CD79A", "CD79B"),
    "Plasma"              = c("MZB1", "JCHAIN", "XBP1", "SDC1"),
    "Monocyte/Macrophage" = c("CD14", "LYZ", "S100A8", "S100A9", "CD68", "FCGR3A"),
    "DC"                  = c("FCER1A", "CLEC10A", "CD1C", "LILRA4", "IRF7"),
    "Platelet"            = c("PPBP", "PF4", "GP9"),
    "RBC"                 = c("HBA1", "HBA2", "HBB")
  )

  # Filter to only genes present in the object
  available_genes <- rownames(obj)
  markers_filtered <- lapply(CANONICAL_MARKERS, function(g) intersect(g, available_genes))
  markers_filtered <- markers_filtered[sapply(markers_filtered, length) > 0]

  if (length(markers_filtered) > 0) {
    safe_plot(bquote(
      FeatureStatPlot(obj, features = .(markers_filtered),
                      ident = "celltype_broad", 
                      plot_type = "dot", 
                      palette = "viridis",
                      row_name_annotation = FALSE)
    ), file.path(vdir, "canonical_markers_celltype_dotplot.pdf"), width = 10, height = 16)
  }
  
  if (length(markers_filtered) > 0) {
    safe_plot(bquote(
      FeatureStatPlot(obj, 
                      features = .(markers_filtered),
                      ident = "merged.celltype.cluster", 
                      plot_type = "dot", 
                      palette = "viridis",
                      row_name_annotation = FALSE)
    ), file.path(vdir, "canonical_markers_cluster_dotplot.pdf"), width = 10, height = 16)
  }

  # Annotation concordance heatmap
  concordance_path <- file.path(paths$results_tables, "cluster_celltype_mapping.csv")
  if (file.exists(concordance_path)) {
    conc_df <- read.csv(concordance_path, check.names = FALSE)
    mode_cols <- grep("^mode_", colnames(conc_df), value = TRUE)
    if (length(mode_cols) >= 2) {
      ref_matrix <- conc_df[, c("cluster", mode_cols, "celltype_broad")]
      safe_plot(bquote({
        plot_df <- .(ref_matrix) %>%
          tidyr::pivot_longer(cols = -cluster, names_to = "reference", values_to = "label") %>%
          mutate(reference = sub("^mode_", "", reference))
        ggplot(plot_df, aes(x = reference, y = as.factor(cluster), fill = label)) +
          geom_tile(color = "white") +
          theme_minimal() +
          labs(title = "Annotation Concordance: References vs Consensus",
               x = "Reference", y = "Cluster", fill = "Label") +
          theme(axis.text.y = element_text(size = 7),
                axis.text.x = element_text(angle = 45, hjust = 1))
      }), file.path(vdir, "annotation_concordance_heatmap.pdf"), width = 12, height = 10)
    }
  }

  # Canonical markers on UMAP (feature plots)
  key_markers <- c("CD3D", "CD8A", "CD4", "FOXP3",
                   "MS4A1", "MZB1", #B/Plasma
                   "CD68", "FCGR1A",  #Monocytes
                   
                   "NKG7", "KLRD1", #KN
                   "LILRA4", "IL3RA") #pDC)
  key_markers <- intersect(key_markers, rownames(obj))
  if (length(key_markers) > 0) {
    safe_plot(bquote(
      FeatureStatPlot(obj, 
                      raster = FALSE,
                      features = .(key_markers), 
                      plot_type = "dim",
                      reduction = "UMAP", 
                      ncol = 4, 
                      palette = "viridis", 
                      bg_cutoff = -Inf, 
                      hex = TRUE) &
        theme(axis.title = element_blank(), 
              axis.text = element_blank(),
              axis.ticks = element_blank())
    ), file.path(vdir, "canonical_markers_featureplot.pdf"), width = 16, height = 10)
  }

}

# ============================================================================
# 4. CLUSTER MARKERS
# ============================================================================
viz_markers <- function(obj, cfg, paths = cfg$paths) {
  log_message("Visualizing: Cluster markers")
  vdir <- viz_subdir(paths, "markers")
  ccol <- .cluster_col(obj)

  markers_path <- file.path(paths$results_tables, "FindAllMarkers_objclusters.csv")
  if (!file.exists(markers_path)) {
    log_message("  No FindAllMarkers results found. Skipping.")
    return(invisible(NULL))
  }

  allmarkers <- read.csv(markers_path, check.names = FALSE)

  # Select top markers per cluster
  top5 <- allmarkers %>%
    dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 5) %>%
    ungroup()
  top5_genes <- unique(top5$gene)

  top8 <- allmarkers %>%
    filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 8) %>%
    ungroup()
  
  top8_genes <- unique(top8$gene)
  top8_genes <- top8_genes[!grepl("LINC|ENSG", top8_genes)]

  # Dot plot: top 5 markers per cluster
  if (length(top5_genes) > 0) {
    safe_plot(bquote(
      
      FeatureStatPlot(obj, 
                      features = .(top8_genes),
                      ident = .(ccol),
                      group_by = .(ccol),
                      plot_type = "dot", 
                      palette = "viridis")
  
    ), file.path(vdir, "markers_dotplot_top5.pdf"), width = 8, height = 20)
  }

  # Heatmap: top 5 markers per cluster
  if (length(top5_genes) > 0) {
    safe_plot(bquote({
      sub_obj <- ScaleData(obj, features = .(top5_genes), verbose = FALSE)
      DoHeatmap(sub_obj, features = .(top5_genes), group.by = .(ccol)) +
        ggtitle("Top 5 Markers per Cluster")
    }), file.path(vdir, "markers_heatmap_top5.pdf"), width = 14, height = 12)
  }

  # Volcano per cluster (faceted)
  safe_plot(bquote({
    df <- .(allmarkers) %>%
      filter(!is.na(p_val_adj)) %>%
      mutate(sig = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
    ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
      geom_point(size = 0.8, alpha = 0.5) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60"), guide = "none") +
      facet_wrap(~cluster, ncol = 4, scales = "free_y") +
      theme_minimal() +
      labs(title = "Volcano Plots per Cluster", x = "avg log2FC", y = "-log10(adj. P)")
  }), file.path(vdir, "markers_volcano_per_cluster.pdf"), width = 20, height = 16)

  # Heatmap of log2FC across clusters
  safe_plot(bquote({
    mat <- .(top8) %>%
      dplyr::select(gene, cluster, avg_log2FC) %>%
      tidyr::pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)
    mat_m <- as.matrix(mat[, -1])
    rownames(mat_m) <- mat$gene
    df_long <- .(top8) %>% dplyr::select(gene, cluster, avg_log2FC)
    ggplot(df_long, aes(x = as.factor(cluster), y = gene, fill = avg_log2FC)) +
      geom_tile(color = "white") +
      scale_fill_viridis() +
      theme_minimal() +
      labs(title = "Marker log2FC Heatmap", x = "Cluster", y = "", fill = "avg log2FC") +
      theme(axis.text.y = element_text(size = 6))
  }), file.path(vdir, "markers_heatmap_log2fc.pdf"), width = 14, height = 14)

  # Dot plot: log2FC + significance across clusters
  safe_plot(bquote({
    df <- .(top8) %>%
      mutate(neg_log_p = -log10(p_val_adj))
    ggplot(df, aes(x = as.factor(cluster), y = gene,
                   size = neg_log_p, color = avg_log2FC)) +
      geom_point() +
      scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027") +
      scale_size_continuous(range = c(1, 6)) +
      theme_minimal() +
      labs(title = "Marker Significance & Effect Size", x = "Cluster", y = "",
           size = "-log10(padj)", color = "avg log2FC") +
      theme(axis.text.y = element_text(size = 6))
  }), file.path(vdir, "markers_dot_log2fc.pdf"), width = 14, height = 14)
}

# ============================================================================
# 5. DIFFERENTIAL GENE EXPRESSION
# ============================================================================
viz_dge <- function(obj, cfg, paths = cfg$paths, target = "all") {
  log_message("Visualizing: Differential gene expression")
  vdir <- viz_subdir(paths, "dge")

  contrasts <- cfg$dge$contrasts
  if (is.null(contrasts)) return(invisible(NULL))
  # Keep contrasts tagged for this target. target == "any" is a wildcard that
  # fires for every target, mirroring run_dge() in R/30_markers.R. Without the
  # "any" branch, "any"-tagged contrasts (e.g. Autoimmune_vs_Viral) were
  # silently dropped from the eye / compartment plot suites even though
  # run_dge had written their DGE tables (regression from commit 23c0029).
  contrasts <- Filter(function(c) {
    ct <- c$target %||% "all"
    ct == target || ct == "any"
  }, contrasts)
  if (length(contrasts) == 0) {
    log_message("  No DGE contrasts for target=", target, ".")
    return(invisible(NULL))
  }

  # Build "immune-related" gene set once: GO:0002376 (immune system process)
  # + all descendants, mapped to HGNC symbols via org.Hs.eg.db. The GOALL
  # keytype expands the parent term to its full subtree, so this picks up
  # cytokine signaling, leukocyte activation, antigen processing, IFN
  # response, complement, etc. Returns NULL on failure so callers can fall
  # back to no immune filtering.
  immune_genes <- tryCatch({
    if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      sym <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys     = "GO:0002376",
        columns  = "SYMBOL",
        keytype  = "GOALL"
      )
      sort(unique(stats::na.omit(sym$SYMBOL)))
    } else NULL
  }, error = function(e) {
    log_message("  WARN: failed to build immune gene set: ",
                conditionMessage(e))
    NULL
  })
  if (!is.null(immune_genes)) {
    log_message("  Immune gene set: ", length(immune_genes),
                " HGNC symbols (GO:0002376 + descendants).")
  }

  for (contrast in contrasts) {
    cname <- contrast$name

    # Try pseudobulk first, then wilcoxon
    for (method in c("pseudobulk", "wilcox")) {
      dge_path <- file.path(paths$results_tables,
                            paste0("DGE_", method, "_", cname, ".csv"))
      if (!file.exists(dge_path)) next

      tryCatch({
        dge_df <- read.csv(dge_path, check.names = FALSE)
        if (nrow(dge_df) == 0) next

        # Standardize column names: DESeq2 uses log2FoldChange/padj,
        # Seurat wilcox uses avg_log2FC/p_val_adj
        fc_col <- if ("log2FoldChange" %in% colnames(dge_df)) "log2FoldChange" else "avg_log2FC"
        pv_col <- if ("padj" %in% colnames(dge_df)) "padj" else "p_val_adj"
        dge_df$fc  <- dge_df[[fc_col]]
        dge_df$pv  <- dge_df[[pv_col]]

        prefix <- paste0("dge_", method, "_", cname)

        # Volcano: global (direction-aware, top labels per direction)
        global_df <- dge_df %>% filter(cluster == "global")
        if (nrow(global_df) > 0) {
          fc_thr <- 0.5; pv_thr <- 0.05; n_label <- 12
          g1 <- contrast$group1; g2 <- contrast$group2
          prepared <- global_df %>%
            filter(!is.na(pv), !is.na(fc)) %>%
            mutate(
              neg_log10_p = -log10(pmax(pv, 1e-300)),
              direction = case_when(
                pv < pv_thr & fc >  fc_thr ~ "up_g1",
                pv < pv_thr & fc < -fc_thr ~ "up_g2",
                TRUE                       ~ "ns"
              ),
              score = neg_log10_p * abs(fc)
            )
          top_up <- prepared %>%
            filter(direction == "up_g1") %>%
            slice_max(score, n = n_label, with_ties = FALSE) %>%
            pull(gene)
          top_dn <- prepared %>%
            filter(direction == "up_g2") %>%
            slice_max(score, n = n_label, with_ties = FALSE) %>%
            pull(gene)
          prepared$label <- ifelse(prepared$gene %in% c(top_up, top_dn),
                                   prepared$gene, NA_character_)
          n_up <- sum(prepared$direction == "up_g1")
          n_dn <- sum(prepared$direction == "up_g2")

          safe_plot(bquote({
            df <- .(prepared)
            df$direction <- factor(df$direction,
                                   levels = c("up_g2", "ns", "up_g1"))
            pal <- c(up_g2 = "#1F78B4", ns = "grey70", up_g1 = "#D7301F")
            lbl <- c(up_g2 = paste("↑ in", .(g2)),
                     ns    = "ns",
                     up_g1 = paste("↑ in", .(g1)))
            ggplot(df, aes(x = fc, y = neg_log10_p, color = direction)) +
              geom_vline(xintercept = c(-.(fc_thr), .(fc_thr)),
                         linetype = "dashed", color = "grey60") +
              geom_hline(yintercept = -log10(.(pv_thr)),
                         linetype = "dashed", color = "grey60") +
              geom_point(size = 1.0, alpha = 0.55) +
              scale_color_manual(values = pal, labels = lbl, name = NULL,
                                 drop = FALSE) +
              ggrepel::geom_text_repel(
                aes(label = label), size = 3, max.overlaps = 30,
                na.rm = TRUE, segment.color = "grey50",
                segment.size = 0.3, min.segment.length = 0.1,
                box.padding = 0.4
              ) +
              annotate("text", x = Inf, y = Inf,
                       hjust = 1.05, vjust = 1.4,
                       label = paste0(.(n_up), " ↑ in ", .(g1)),
                       color = "#D7301F", size = 3.4, fontface = "bold") +
              annotate("text", x = -Inf, y = Inf,
                       hjust = -0.05, vjust = 1.4,
                       label = paste0(.(n_dn), " ↑ in ", .(g2)),
                       color = "#1F78B4", size = 3.4, fontface = "bold") +
              theme_minimal() +
              theme(legend.position = "top") +
              labs(
                title = paste("Volcano (global):", .(cname), "-", .(method)),
                subtitle = paste0("|log2FC| > ", .(fc_thr),
                                  " and padj < ", .(pv_thr),
                                  "; top ", .(n_label),
                                  " labelled per direction by ",
                                  "-log10(padj) × |log2FC|"),
                x = paste0("log2 Fold Change   ← ↑ in ", .(g2),
                           "   |   ↑ in ", .(g1), " →"),
                y = "-log10(adj. P)"
              )
          }), file.path(vdir, paste0(prefix, "_volcano_global.pdf")),
          width = 6, height = 6)
        }

        # Per-cluster volcano (faceted)
        cluster_df <- dge_df %>% filter(cluster != "global")
        if (nrow(cluster_df) > 0) {
          safe_plot(bquote({
            df <- .(cluster_df) %>%
              filter(!is.na(pv)) %>%
              mutate(sig = pv < 0.05 & abs(fc) > 0.5)
            ggplot(df, aes(x = fc, y = -log10(pv), color = sig)) +
              geom_point(size = 0.6, alpha = 0.4) +
              scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60"), guide = "none") +
              facet_wrap(~cluster, ncol = 4, scales = "free_y") +
              theme_minimal() +
              labs(title = paste("Volcano per Cluster:", .(cname)),
                   x = "log2 Fold Change", y = "-log10(adj. P)")
          }), file.path(vdir, paste0(prefix, "_volcano_per_cluster.pdf")),
          width = 20, height = 16)

          # Heatmap of log2FC across clusters (top 10 per cluster)
          safe_plot(bquote({
            top_per <- .(cluster_df) %>%
              filter(!is.na(pv), pv < 0.05) %>%
              group_by(cluster) %>%
              slice_max(abs(fc), n = 10) %>%
              ungroup()
            ggplot(top_per, aes(x = as.factor(cluster), y = gene, fill = fc)) +
              geom_tile(color = "white") +
              scale_fill_viridis() +
              theme_minimal() +
              labs(title = paste("DEG log2FC Heatmap:", .(cname)),
                   x = "Cluster", y = "", fill = "log2FC") +
              theme(axis.text.y = element_text(size = 5))
          }), file.path(vdir, paste0(prefix, "_heatmap_log2fc.pdf")),
          width = 14, height = 12)

          # Compartment contribution decomposition (REPLACES the old consensus
          # bubble). The bubble colored a gene x cluster grid by within-cluster
          # log2FC, so the paracrine ISG block lit up in every compartment and
          # the panel read as "the antiviral signature is everywhere", arguing
          # against the localization claim in the text. Instead, for each
          # leading gene we decompose its bulk signal into the fraction of raw
          # counts contributed by each lineage, separately for the viral-up and
          # autoimmune-up programs. This answers where each program actually
          # lives. Only emitted for the pseudobulk method (the file the
          # manuscript Fig 2G is built from).
          if (method == "pseudobulk") tryCatch({
            g1 <- contrast$group1; g2 <- contrast$group2   # NIU, Viral
            ccol    <- .cluster_col(obj)
            # Counts must come from the gene-level RNA assay, not SCT/integrated.
            assay   <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
            present <- rownames(obj[[assay]])

            # Curated functional submodules per direction. Kept only if present
            # in the object AND significant in the matching direction of the
            # GLOBAL contrast (fc/pv added to dge_df upstream).
            mods_viral <- list(
              ISG               = c("IFI27","IFI6","IFI44L","IFIT2","IFIT3","ISG15","MX1","RSAD2"),
              `CXCR3 chemokine` = c("CXCL9","CXCL10","CXCL11"),
              `Mono chemokine`  = c("CCL2","CCL7","CCL8","CCL18"),
              Complement        = c("C1QC","CFB","SERPING1"),
              Checkpoint        = c("LAG3","IL10","CD274"))
            mods_auto <- list(
              `DC / cDC2`      = c("FCER1A","CD1C","CD1E","CD1A","CLEC10A"),
              `Classical mono` = c("FCN1","S100A12","FOLR3","CD300H","ADGRE1","CTSG","VCAN","SELL"),
              `Th17 / IL-23`   = c("RORC","IL23A","CCR6","IL22"),
              `T costim`       = c("TNFRSF8","CD40LG"),
              `B lineage`      = c("CD79B","MS4A1","TCL1A"))

            glob  <- dge_df %>% filter(cluster == "global", !is.na(pv), !is.na(fc))
            up_g1 <- glob %>% filter(fc >  0.5, pv < 0.05) %>% pull(gene)  # NIU up
            up_g2 <- glob %>% filter(fc < -0.5, pv < 0.05) %>% pull(gene)  # Viral up

            # Collapse fine cell-type labels to broad lineages. Extend the
            # regexes if your label scheme differs.
            lineage_of <- function(x) dplyr::case_when(
              grepl("Mono|Macro|DC|cDC|pDC|Langerhans|Myeloid", x) ~ "Myeloid",
              grepl("^CD8", x)                                     ~ "CD8 T",
              grepl("^CD4|Treg", x)                                ~ "CD4 T",
              grepl("gdT|MAIT|dnT", x)                             ~ "Other T",
              grepl("NK", x)                                       ~ "NK",
              grepl("Plasma|^B$|^B[ /_-]", x)                      ~ "B/Plasma",
              TRUE                                                 ~ "Mixed")
            lin_levels <- c("Myeloid","CD8 T","CD4 T","Other T","NK","B/Plasma","Mixed")
            lin_pal <- c(Myeloid = "#E4572E", "CD8 T" = "#2E86AB", "CD4 T" = "#5BC0BE",
                         "Other T" = "#9BC53D", NK = "#7B2CBF", "B/Plasma" = "#F2A65A",
                         Mixed = "#B0B0B0")

            # Per-cluster summed counts via a sparse indicator matrix. This
            # keeps exact cluster labels (no Seurat name sanitization) and stays
            # sparse until the final small genes x clusters result.
            grp  <- factor(as.character(obj[[ccol]][, 1]))
            cnts <- tryCatch(
              SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
              error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = "counts"))
            ind  <- Matrix::fac2sparse(grp)                # levels x cells

            build_long <- function(mods, sig_genes) {
              pres <- lapply(mods, function(v) intersect(intersect(v, present), sig_genes))
              pres <- pres[lengths(pres) > 0]
              genes <- unlist(pres, use.names = FALSE)
              if (length(genes) < 2) return(NULL)
              submod_of <- stats::setNames(rep(names(pres), lengths(pres)), genes)
              agg <- as.matrix(cnts[genes, , drop = FALSE] %*% Matrix::t(ind))
              colnames(agg) <- levels(grp)
              frac <- sweep(agg, 1, pmax(rowSums(agg), 1), "/")
              tibble::as_tibble(frac, rownames = "gene") %>%
                tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "frac") %>%
                mutate(lineage   = factor(lineage_of(cluster), lin_levels),
                       submodule = factor(submod_of[gene], names(pres))) %>%
                group_by(gene, submodule, lineage) %>%
                summarise(frac = sum(frac), .groups = "drop") %>%
                mutate(gene = factor(gene, rev(genes)))
            }

            viral_long <- build_long(mods_viral, up_g2)    # up in g2 (Viral)
            auto_long  <- build_long(mods_auto,  up_g1)    # up in g1 (NIU)

            if (!is.null(viral_long) || !is.null(auto_long)) {
              n_rows <- length(unique(c(as.character(viral_long$gene),
                                        as.character(auto_long$gene))))
              safe_plot(bquote({
                vl <- .(viral_long); al <- .(auto_long); pal <- .(lin_pal)
                mk <- function(d, ttl) {
                  if (is.null(d)) return(patchwork::plot_spacer())
                  ggplot(d, aes(frac, gene, fill = lineage)) +
                    geom_col(width = 0.8, color = "white", linewidth = 0.2) +
                    geom_vline(xintercept = 0.5, linetype = "dashed",
                               color = "grey50") +
                    scale_fill_manual(values = pal, name = "Lineage", drop = FALSE) +
                    facet_grid(submodule ~ ., scales = "free_y",
                               space = "free_y", switch = "y") +
                    theme_minimal(base_size = 10) +
                    theme(strip.text.y.left = element_text(angle = 0, face = "bold"),
                          panel.grid.minor = element_blank()) +
                    labs(title = ttl, x = "Fraction of total expression", y = NULL)
                }
                (mk(vl, paste("Up in", .(g2), "- antiviral effector program")) /
                 mk(al, paste("Up in", .(g1), "- autoimmune program"))) +
                  patchwork::plot_layout(guides = "collect") +
                  patchwork::plot_annotation(
                    title = paste("Compartment contribution to the bulk signal:",
                                  .(cname)),
                    subtitle = paste0("Fraction of each gene's total counts ",
                                      "contributed by each lineage"))
              }), file.path(vdir, paste0(prefix, "_consensus_bubble.pdf")),
              width = 9, height = max(7, 0.26 * n_rows + 3))

              if (!is.null(viral_long))
                write.csv(viral_long, file.path(paths$results_tables,
                  paste0("dge_", method, "_", cname, "_contribution_viral.csv")),
                  row.names = FALSE)
              if (!is.null(auto_long))
                write.csv(auto_long, file.path(paths$results_tables,
                  paste0("dge_", method, "_", cname, "_contribution_autoimmune.csv")),
                  row.names = FALSE)
              log_message("  Saved: ", prefix,
                          "_consensus_bubble.pdf (compartment contribution) + CSVs.")
            } else {
              log_message("  Contribution panel: no qualifying genes for ",
                          cname, "/", method, "; skipping.")
            }
          }, error = function(e) {
            log_message("  WARN: contribution panel (", cname, "/", method,
                        ") failed: ", conditionMessage(e))
          })
        }

        # MA plot (baseMean vs log2FC) for pseudobulk global
        if (method == "pseudobulk" && "baseMean" %in% colnames(dge_df)) {
          global_ma <- dge_df %>% filter(cluster == "global", !is.na(pv))
          if (nrow(global_ma) > 0) {
            safe_plot(bquote({
              df <- .(global_ma) %>%
                mutate(sig = pv < 0.05 & abs(fc) > 0.5)
              ggplot(df, aes(x = log10(baseMean + 1), y = fc, color = sig)) +
                geom_point(size = 0.8, alpha = 0.5) +
                scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60"),
                                   name = "Significant") +
                geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
                theme_minimal() +
                labs(title = paste("MA Plot:", .(cname)),
                     x = "log10(Mean Expression)", y = "log2 Fold Change")
            }), file.path(vdir, paste0(prefix, "_ma_plot.pdf")), width = 8, height = 6)
          }
        }

        # Top DEGs violin per celltype (for primary method only)
        if (method == "pseudobulk" && "celltype_broad" %in% colnames(obj[[]])) {
          top_genes <- dge_df %>%
            filter(cluster == "global", !is.na(pv)) %>%
            arrange(pv) %>%
            head(8) %>%
            pull(gene)

          top_genes <- intersect(top_genes, rownames(obj))
          if (length(top_genes) > 0) {
            safe_plot(bquote(
              FeatureStatPlot(obj, 
                              raster = FALSE,
                              features = .(top_genes),
                              ident = "celltype_broad",
                              group_by = .(contrast$group_col),
                              plot_type = "violin")
            ), file.path(vdir, paste0(prefix, "_top_degs_violin_celltype.pdf")),
            width = 16, height = 12)

            # Feature dim split by comparison variable
            safe_plot(bquote(
              FeatureStatPlot(obj, 
                              raster = FALSE,
                              features = .(top_genes[1:min(4, length(top_genes))]),
                              plot_type = "dim",
                              split_by = .(contrast$group_col))
            ), file.path(vdir, paste0(prefix, "_top_degs_dimplot.pdf")),
            width = 14, height = 10)
          }
        }
      }, error = function(e) {
        log_message("  WARN DGE viz failed (", method, " ", cname, "): ", conditionMessage(e))
      })
    }
  }

  # DGE summary barplot (DEG counts per cluster)
  summary_path <- file.path(paths$results_tables, "dge_summary_stats.csv")
  if (file.exists(summary_path)) {
    dge_stats <- read.csv(summary_path, check.names = FALSE)
    for (cname_s in unique(dge_stats$contrast)) {
      cstats <- dge_stats %>%
        filter(contrast == cname_s, cluster != "global")
      if (nrow(cstats) > 0) {
        safe_plot(bquote({
          plot_df <- .(cstats) %>%
            tidyr::pivot_longer(cols = c(n_up, n_down), names_to = "direction",
                                values_to = "count") %>%
            mutate(direction = ifelse(direction == "n_up", "Up", "Down"),
                   count = ifelse(direction == "Down", -count, count))
          ggplot(plot_df, aes(x = as.factor(cluster), y = count, fill = direction)) +
            geom_col() +
            scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
            coord_flip() +
            theme_minimal() +
            labs(title = paste("DEG Counts per Cluster:", .(cname_s)),
                 x = "Cluster", y = "Number of DEGs (Up / Down)", fill = "Direction") +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
        }), file.path(vdir, paste0("dge_summary_barplot_", cname_s, ".pdf")),
        width = 10, height = 8)
      }
    }
  }

  # UpSet plot of shared DEGs across clusters
  for (contrast in contrasts) {
    cname_up <- contrast$name
    pb_path <- file.path(paths$results_tables, paste0("DGE_pseudobulk_", cname_up, ".csv"))
    if (!file.exists(pb_path)) next

    pb_df <- read.csv(pb_path, check.names = FALSE)
    cluster_sig <- pb_df %>%
      filter(cluster != "global", !is.na(padj), padj < 0.05) %>%
      distinct(cluster, gene)

    clusters_with_degs <- unique(cluster_sig$cluster)
    if (length(clusters_with_degs) < 2) next

    tryCatch({
      library(UpSetR)
      upset_list <- split(cluster_sig$gene, cluster_sig$cluster)
      upset_list <- upset_list[sapply(upset_list, length) > 0]

      if (length(upset_list) >= 2) {
        ensure_dir(vdir)  # raw pdf() needs the bucket to exist
        pdf(file.path(vdir, paste0("dge_upset_", cname_up, ".pdf")),
            width = 12, height = 7)
        print(upset(fromList(upset_list), order.by = "freq",
                    nsets = min(20, length(upset_list)),
                    text.scale = 1.3, mainbar.y.label = "Shared DEGs",
                    sets.x.label = "DEGs per Cluster"))
        dev.off()
        log_message("  Saved: dge_upset_", cname_up, ".pdf")
      }
    }, error = function(e) {
      log_message("  WARN UpSet plot failed: ", conditionMessage(e))
    })
  }
}

# ============================================================================
# 7. MILOR DIFFERENTIAL ABUNDANCE
# ============================================================================
viz_milo <- function(obj, cfg, paths = cfg$paths, target = "all") {
  log_message("Visualizing: MiloR differential abundance")
  vdir <- viz_subdir(paths, "milo")

  suppressPackageStartupMessages({
    library(miloR)
    library(SingleCellExperiment)
    library(ggbeeswarm)
    library(ggrepel)
    library(viridis)
    library(RColorBrewer)
    library(ggnewscale)
  })

  milo_path <- file.path(paths$results_objects, "MiloObject.rds")
  # Annotation columns come from 40_milo.R::annotateNhoods. The eye object
  # and the three compartment sub-atlases skip merge_clusters, so the
  # per-cluster grouping is the fine-grained knn.leiden.cluster (re-clustered
  # on the eye-only / compartment-only graph).
  is_compartment <- target %in% c("myeloid", "bcell", "tcell")
  groupings <- if (target == "eye" || is_compartment) {
    c("knn.leiden.cluster")
  } else {
    c("merged.celltype.cluster")
  }
  contrast_names <- if (target == "eye") {
    c("etiology_eye", "granulom_eye")
  } else if (is_compartment) {
    # 40_milo.R appends "_eye" to every contrast in the eye/compartment loop,
    # so files on disk are e.g. Milo_etiology_myeloid_eye_DA.csv. Match that.
    paste0("etiology_", target, "_eye")
  } else {
    c("tissue")
  }

  # --- helpers ----------------------------------------------------------------
  stars_label <- function(q) {
    ifelse(is.na(q), "",
      ifelse(q < 0.001, "***",
      ifelse(q < 0.01,  "**",
      ifelse(q < 0.05,  "*",
      ifelse(q < 0.1,   ".", "ns")))))
  }

  paired_palette <- function(n) {
    base <- RColorBrewer::brewer.pal(12, "Paired")
    if (n <= 12) base[seq_len(n)] else colorRampPalette(base)(n)
  }

  drop_edge_layers <- function(p) {
    keep <- vapply(p$layers, function(l) {
      !any(grepl("edge", class(l$geom), ignore.case = TRUE))
    }, logical(1))
    p$layers <- p$layers[keep]
    p
  }

  # ---------------------------------------------------------------------------

  for (contrast_name in contrast_names) {
    res_path <- file.path(paths$results_tables,
                          paste0("Milo_", contrast_name, "_DA.csv"))
    if (!file.exists(res_path)) next
    res <- read.csv(res_path, check.names = FALSE)
    if (!all(c("logFC", "SpatialFDR", "Nhood") %in% colnames(res))) next
    ensure_dir(vdir)  # raw ggsave() calls below need the bucket to exist

    # -------- Beeswarm + Boxplot for each grouping ---------------------------
    for (gcol in groupings) {
      if (!gcol %in% colnames(res)) {
        log_message("  Skipping ", gcol, " for ", contrast_name, " (column missing).")
        next
      }
      purity_col <- paste0(gcol, "_fraction")

      df <- res %>% filter(!is.na(.data[[gcol]]))
      if (purity_col %in% colnames(df)) {
        df <- df %>% filter(.data[[purity_col]] >= 0.7)
      }
      if (nrow(df) < 5) {
        log_message("  Too few nhoods after purity filter for ", gcol,
                    " / ", contrast_name, "; skipping.")
        next
      }

      # Wilcoxon signed-rank vs 0 per group, BH-adjusted across groups
      stats_df <- df %>%
        group_by(.data[[gcol]]) %>%
        summarise(
          n            = dplyr::n(),
          n_sig        = sum(SpatialFDR < 0.1, na.rm = TRUE),
          n_up_sig     = sum(SpatialFDR < 0.1 & logFC > 0, na.rm = TRUE),
          n_down_sig   = sum(SpatialFDR < 0.1 & logFC < 0, na.rm = TRUE),
          median_logFC = median(logFC, na.rm = TRUE),
          wilcox_p     = if (dplyr::n() >= 3)
            tryCatch(wilcox.test(logFC, mu = 0, exact = FALSE)$p.value,
                     error = function(e) NA_real_)
            else NA_real_,
          .groups = "drop"
        ) %>%
        mutate(
          wilcox_q  = p.adjust(wilcox_p, method = "BH"),
          direction = ifelse(is.na(wilcox_q) | wilcox_q >= 0.1, "",
                      ifelse(median_logFC > 0, "\u2191", "\u2193")),
          sig_label = trimws(paste(direction, stars_label(wilcox_q)))
        )

      # Order groups by median logFC (most negative bottom, most positive top)
      group_order <- stats_df %>% arrange(median_logFC) %>% pull(.data[[gcol]])
      df[[gcol]]       <- factor(df[[gcol]],       levels = group_order)
      stats_df[[gcol]] <- factor(stats_df[[gcol]], levels = group_order)

      # Annotation x-position (just right of max logFC)
      xrange   <- range(df$logFC, na.rm = TRUE)
      x_ann    <- xrange[2] + diff(xrange) * 0.08
      stats_df$x_ann <- x_ann

      n_groups <- length(group_order)
      plt_h    <- max(4, 0.35 * n_groups + 2.5)

      # --- Boxplot (fill = group, Paired palette) ---
      tryCatch({
        cols <- paired_palette(n_groups)
        names(cols) <- group_order

        p_bx <- ggplot(df, aes(x = logFC, y = .data[[gcol]],
                               fill = .data[[gcol]])) +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
          geom_boxplot(outlier.shape = NA, alpha = 0.85, color = "grey20") +
          scale_fill_manual(values = cols, guide = "none") +
          geom_text(data = stats_df,
                    aes(x = x_ann, y = .data[[gcol]], label = sig_label),
                    inherit.aes = FALSE, hjust = 0, size = 3.3) +
          coord_cartesian(clip = "off") +
          theme_minimal() +
          theme(plot.margin = margin(5, 55, 5, 5)) +
          labs(title = paste("MiloR DA:", contrast_name, "by", gcol),
               subtitle = "Wilcoxon signed-rank vs 0, BH-adjusted (\u2191/\u2193 sig@q<0.1)",
               x = "logFC", y = NULL)

        ggsave(file.path(vdir, paste0("milo_", contrast_name,
                                      "_boxplot_", gcol, ".pdf")),
               p_bx, width = 10, height = plt_h)
        log_message("  Saved: milo_", contrast_name, "_boxplot_", gcol, ".pdf")
      }, error = function(e) {
        log_message("  WARN: boxplot (", gcol, ") failed: ",
                    conditionMessage(e))
      })

      # --- Summary boxplot + FDR strip (one ggplot, dual fill via ggnewscale) ---
      tryCatch({
        cols <- paired_palette(n_groups)
        names(cols) <- group_order

        # Right-side strip geometry (just outside the boxplot data range)
        strip_w <- diff(xrange) * 0.05
        strip_x <- xrange[2] + diff(xrange) * 0.10
        star_x  <- strip_x + strip_w * 0.7 + diff(xrange) * 0.02

        sdf <- stats_df %>%
          mutate(
            neg_log10_q = -log10(pmax(wilcox_q, 1e-300)),
            strip_x     = strip_x,
            star_x      = star_x
          )

        p_sum <- ggplot(df, aes(x = logFC, y = .data[[gcol]])) +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
          geom_boxplot(aes(fill = .data[[gcol]]),
                       outlier.shape = NA, alpha = 0.9, color = "grey20") +
          scale_fill_manual(values = cols, guide = "none") +
          ggnewscale::new_scale_fill() +
          geom_tile(data = sdf,
                    aes(x = strip_x, y = .data[[gcol]], fill = neg_log10_q),
                    width = strip_w, height = 0.75,
                    color = "grey30", inherit.aes = FALSE) +
          viridis::scale_fill_viridis(
            option = "rocket", direction = -1,
            name = expression(-log[10]*"(q)"),
            limits = c(0, NA), na.value = "grey90"
          ) +
          geom_text(data = sdf,
                    aes(x = star_x, y = .data[[gcol]], label = sig_label),
                    inherit.aes = FALSE, hjust = 0, size = 3.3) +
          coord_cartesian(clip = "off") +
          theme_minimal() +
          theme(plot.margin = margin(5, 90, 5, 5)) +
          labs(
            title = paste("MiloR DA:", contrast_name, "by", gcol),
            subtitle = paste0(
              "Box: per-nhood logFC; tile: cluster-level Wilcoxon q (BH); ",
              "stars: q<0.001/0.01/0.05/0.1"
            ),
            x = "logFC", y = NULL
          )

        ggsave(file.path(vdir, paste0("milo_", contrast_name,
                                      "_summary_", gcol, ".pdf")),
               p_sum, width = 11, height = plt_h)
        log_message("  Saved: milo_", contrast_name, "_summary_", gcol, ".pdf")
      }, error = function(e) {
        log_message("  WARN: summary plot (", gcol, ") failed: ",
                    conditionMessage(e))
      })

      # Per-group stats table
      write.csv(stats_df %>% dplyr::select(-x_ann),
                file.path(paths$results_tables,
                          paste0("milo_", contrast_name, "_",
                                 gcol, "_wilcox_stats.csv")),
                row.names = FALSE)

      # --- Normalized stacked composition bar (milo-ordered + Fisher/GLMM q) ---
      # Maps the milo contrast_name back to the cell-level contrast metadata
      # (column + reference/alt groups). The reference/alt order matches
      # 40_milo.R::eye_contrasts so the GLMM beta sign tracks milo's logFC.
      contrast_meta <- list(
        tissue                = list(col = "Tissue_1",    groups = c("Eye",  "Blood")),
        etiology_eye          = list(col = "Phenotype_2", groups = c("NIU",  "Viral")),
        granulom_eye          = list(col = "Phenotype_2", groups = c("Gran", "Non_Gran")),
        etiology_myeloid_eye  = list(col = "Phenotype_2", groups = c("NIU",  "Viral")),
        etiology_bcell_eye    = list(col = "Phenotype_2", groups = c("NIU",  "Viral")),
        etiology_tcell_eye    = list(col = "Phenotype_2", groups = c("NIU",  "Viral"))
      )
      cm <- contrast_meta[[contrast_name]]
      if (!is.null(cm) && cm$col %in% colnames(obj[[]]) &&
          gcol %in% colnames(obj[[]]) &&
          "Subject" %in% colnames(obj[[]])) {
        tryCatch({
          md <- obj[[]]
          keep <- md[[cm$col]] %in% cm$groups & !is.na(md[[gcol]])
          md_k <- md[keep, , drop = FALSE]
          md_k[[cm$col]] <- factor(as.character(md_k[[cm$col]]),
                                   levels = cm$groups)
          # Restrict to clusters present in the milo group_order so the bar
          # plot lines up vertically with the milo boxplot.
          md_k <- md_k[as.character(md_k[[gcol]]) %in%
                         as.character(group_order), , drop = FALSE]
          md_k[[gcol]] <- factor(as.character(md_k[[gcol]]),
                                 levels = group_order)

          if (nrow(md_k) > 0) {
            N_cond <- as.numeric(table(md_k[[cm$col]]))
            names(N_cond) <- levels(md_k[[cm$col]])

            cnt <- as.data.frame(table(md_k[[gcol]], md_k[[cm$col]]))
            colnames(cnt) <- c("cluster", "group", "n")
            cnt$prop_in_cond <- cnt$n / N_cond[as.character(cnt$group)]
            cnt <- cnt %>%
              dplyr::group_by(cluster) %>%
              dplyr::mutate(
                denom = sum(prop_in_cond, na.rm = TRUE),
                frac  = ifelse(denom > 0, prop_in_cond / denom, NA_real_)
              ) %>%
              dplyr::ungroup()
            cnt$cluster <- factor(as.character(cnt$cluster),
                                  levels = group_order)
            cnt$group   <- factor(as.character(cnt$group), levels = cm$groups)

            fg <- run_fisher_glmm_per_cluster(
              meta         = md_k,
              contrast_col = cm$col,
              groups       = cm$groups,
              cluster_col  = gcol,
              sample_col   = "Subject"
            )
            fg$cluster <- factor(as.character(fg$cluster),
                                 levels = group_order)
            fg$q_show  <- ifelse(is.na(fg$q_glmm), fg$q_fisher, fg$q_glmm)
            fg$lbl     <- trimws(paste0(
              stars_label(fg$q_show),
              ifelse(is.na(fg$q_show), "",
                     sprintf(" q=%.2g", fg$q_show))
            ))

            write.csv(fg,
                      file.path(paths$results_tables,
                                paste0("milo_", contrast_name, "_", gcol,
                                       "_fisher_glmm.csv")),
                      row.names = FALSE)

            # Palette: use the project NIU/Viral hex codes when this is an
            # etiology contrast so the bar matches viz_compartment_composition
            # downstream. Falls back to Set1 for tissue / granulom / unknown
            # contrasts whose levels aren't in ETIOLOGY_GROUP_COLORS.
            grp_pal <- if (all(cm$groups %in% names(ETIOLOGY_GROUP_COLORS))) {
              ETIOLOGY_GROUP_COLORS[cm$groups]
            } else {
              setNames(
                RColorBrewer::brewer.pal(max(3, length(cm$groups)),
                                         "Set1")[seq_along(cm$groups)],
                cm$groups
              )
            }

            p_sb <- ggplot(cnt, aes(x = frac, y = cluster, fill = group)) +
              geom_col(position = "stack", color = "grey25", width = 0.78) +
              geom_vline(xintercept = 0.5,
                         linetype = "dashed", color = "grey30") +
              scale_fill_manual(values = grp_pal, name = NULL) +
              scale_x_continuous(
                expand = c(0, 0), limits = c(0, 1.20),
                breaks = c(0, 0.25, 0.5, 0.75, 1),
                labels = scales::percent_format(accuracy = 1)
              ) +
              geom_text(data = fg,
                        aes(x = 1.02, y = cluster, label = lbl),
                        inherit.aes = FALSE, hjust = 0, size = 3.0) +
              coord_cartesian(clip = "off") +
              theme_minimal() +
              theme(panel.grid.major.y = element_blank(),
                    plot.margin = margin(5, 95, 5, 5)) +
              labs(
                title = paste0("Normalized composition: ", contrast_name,
                               " by ", gcol),
                subtitle = paste0(
                  "Per-condition counts normalized by total cells in each ",
                  cm$col, " group, then rescaled to fill per cluster; ",
                  "50% line = baseline-no-DA reference. q from binomial GLMM",
                  " (cbind(n_in,n_out) ~ ", cm$col,
                  " + (1|Subject), LRT, BH); falls back to Fisher q if GLMM NA."
                ),
                x = "Phenotype-normalised share (per cluster)", y = NULL
              )

            ggsave(file.path(vdir, paste0("milo_", contrast_name,
                                          "_composition_", gcol, ".pdf")),
                   p_sb, width = 11, height = plt_h)
            log_message("  Saved: milo_", contrast_name,
                        "_composition_", gcol, ".pdf")
          } else {
            log_message("  Composition: no overlapping cells for ",
                        contrast_name, " / ", gcol, "; skipping.")
          }
        }, error = function(e) {
          log_message("  WARN: composition plot (", gcol, ") failed: ",
                      conditionMessage(e))
        })
      }
    }

    # -------- Volcano (viridis) ---------------------------------------------
    safe_plot(bquote({
      df <- .(res) %>% filter(!is.na(SpatialFDR), !is.na(logFC))
      ggplot(df, aes(x = logFC, y = -log10(SpatialFDR))) +
        geom_point(aes(color = -log10(SpatialFDR)), size = 1.5, alpha = 0.75) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = "-log10(SpatialFDR)") +
        geom_hline(yintercept = -log10(0.1),
                   linetype = "dashed", color = "red") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        theme_minimal() +
        labs(title = paste("MiloR DA:", .(contrast_name)),
             x = "Log Fold Change", y = "-log10(SpatialFDR)")
    }), file.path(vdir, paste0("milo_", contrast_name, "_volcano.pdf")),
    width = 8, height = 6)

  }

  # -------- Nhood size distribution (unchanged) -----------------------------
  if (file.exists(milo_path)) {
    tryCatch({
      milo <- readRDS(milo_path)
      nhood_sizes <- colSums(nhoods(milo))

      safe_plot(bquote({
        ggplot(data.frame(size = .(nhood_sizes)), aes(x = size)) +
          geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
          theme_minimal() +
          labs(title = "MiloR Neighborhood Size Distribution",
               x = "Neighborhood Size (cells)", y = "Count") +
          geom_vline(xintercept = median(.(nhood_sizes)),
                     linetype = "dashed", color = "red")
      }), file.path(vdir, "milo_nhood_size_hist.pdf"), width = 7, height = 5)

      rm(milo)
    }, error = function(e) {
      log_message("  WARN: nhood size hist failed: ", conditionMessage(e))
    })
  }
}

# ============================================================================
# 12. COMPOSITION TESTING
# ============================================================================
viz_composition <- function(cfg, paths = cfg$paths, target = "all") {
  log_message("Visualizing: Composition testing")
  vdir <- viz_subdir(paths, "composition")

  comp_names <- if (target == "eye") {
    c("Autoimmune_vs_Viral", "Gran_vs_NonGran")
  } else {
    c("Eye_vs_Blood")
  }

  for (cname in comp_names) {
    res_path <- file.path(paths$results_tables, paste0("composition_test_", cname, ".csv"))
    if (!file.exists(res_path)) next

    res <- read.csv(res_path, check.names = FALSE, row.names = 1)

    # Significance barplot
    if ("FDR" %in% colnames(res)) {
      safe_plot(bquote({
        df <- .(res) %>%
          tibble::rownames_to_column("celltype") %>%
          mutate(sig = FDR < 0.05,
                 neg_log_fdr = -log10(FDR))
        ggplot(df, aes(x = reorder(celltype, neg_log_fdr), y = neg_log_fdr, fill = sig)) +
          geom_col(alpha = 0.8) +
          scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "grey60"),
                            name = "FDR < 0.05") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
          coord_flip() +
          theme_minimal() +
          labs(title = paste("Composition Test:", .(cname)),
               x = "Cell Type", y = "-log10(FDR)")
      }), file.path(vdir, paste0("composition_test_barplot_", cname, ".pdf")),
      width = 8, height = 6)
    }
  }
}


# ============================================================================
# 13. BCR LINEAGE
# ============================================================================
viz_bcr_lineage <- function(cfg, paths = cfg$paths) {
  log_message("Visualizing: BCR lineage trees")
  # Trees are already generated in 52_bcr_lineage.R
  tree_dir <- viz_subdir(paths, "bcr_lineage")
  if (dir.exists(tree_dir)) {
    n_trees <- length(list.files(tree_dir, pattern = "\\.pdf$"))
    log_message("  BCR lineage trees directory: ", n_trees, " tree PDFs found.")
  } else {
    log_message("  No BCR lineage tree directory found.")
  }
}
