# R/89_viz_tcr_advanced.R
# Visualization dispatcher for the advanced TCR analysis modules (67-71).
# All panels write to outputs/viz/eye/tcell/09_tcr_motif/tcr_advanced/ and reuse the
# shared helpers from R/01_setup_utils.R, R/81_viz_compartment_helpers.R,
# and R/88_viz_tcell.R.
#
# Panels:
#   viz_gliph_tcrdist_correlation()  — UMAP_TCR + overlap heatmap + violin
#   viz_gliph_network_metrics()      — per-cluster network + metric heatmap
#   viz_vdjdb_overlay()              — per-species UMAP overlays + GLIPH crosswalk
#   viz_tcr_genex_signatures()       — UCell heatmap, volcano, coupling diagram
#   viz_meta_motifs()                — VGKWY-style assembly figure + UMAP overlay

`%||%` <- function(x, y) if (is.null(x)) y else x

.tcra_paths <- function(cfg) {
  paths_tcell <- get_target_paths(cfg, "tcell")
  base <- file.path(viz_subdir(paths_tcell, "tcr_motif"), "tcr_advanced")
  ensure_dir(base)
  list(
    base       = base,
    tables_rep = file.path(cfg$paths$results_tables, "repertoire"),
    tables_t   = paths_tcell$results_tables,
    objs       = cfg$paths$results_objects,
    obj_t      = paths_tcell$results_objects
  )
}

# ---------------------------------------------------------------------------
# Module 1: GLIPH x tcrdist correlation
# ---------------------------------------------------------------------------
viz_gliph_tcrdist_correlation <- function(cfg) {
  p <- .tcra_paths(cfg)
  joint_rds <- file.path(p$objs, "GliphTcrdistJoint.rds")
  conc_csv  <- file.path(p$tables_rep, "gliph_tcrdist_concordance.csv")
  ov_csv    <- file.path(p$tables_rep, "gliph_clustcr_overlap.csv")
  tcr_rds   <- file.path(p$objs, "ImmLynxTcrdistResults.rds")
  if (!file.exists(joint_rds) || !file.exists(conc_csv)) {
    log_message("  viz_gliph_tcrdist_correlation: required tables missing; skipping.")
    return(invisible(NULL))
  }

  joint <- readRDS(joint_rds)
  conc  <- utils::read.csv(conc_csv, stringsAsFactors = FALSE)

  # --- UMAP_TCR colored by primary GLIPH cluster (top 12 by cell count) ----
  if (file.exists(tcr_rds)) {
    tcr <- readRDS(tcr_rds)
    umap <- tcr$umap
    if (!is.null(umap)) {
      bcs <- rownames(umap)
      df <- data.frame(barcode = bcs,
                       UMAP1   = umap[, 1],
                       UMAP2   = umap[, 2],
                       stringsAsFactors = FALSE) |>
        dplyr::left_join(
          dplyr::select(joint, barcode, gliph_primary_cluster), by = "barcode")
      top_clusters <- names(sort(table(df$gliph_primary_cluster),
                                 decreasing = TRUE))[1:12]
      df$cluster_lbl <- ifelse(df$gliph_primary_cluster %in% top_clusters,
                               df$gliph_primary_cluster, "other/none")
      df$cluster_lbl[is.na(df$gliph_primary_cluster)] <- "none"
      p1 <- ggplot2::ggplot(df,
                            ggplot2::aes(UMAP1, UMAP2, color = cluster_lbl)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.7) +
        ggplot2::scale_color_viridis_d(option = "turbo", na.value = "grey85") +
        ggplot2::labs(title  = "tcrdist UMAP colored by GLIPH cluster (top 12)",
                      color  = "GLIPH cluster") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(legend.position = "right")
      save_pdf_png(p1, file.path(p$base, "tcrdist_umap_by_gliph"),
                   w = 8, h = 7)
    }
  }

  # --- GLIPH x clusTCR log2OR tile heatmap --------------------------------
  # With ~500 GLIPH clusters x ~180 clusTCR clusters every tick label
  # collapses into noise. Restrict to FDR-significant pairs (with a
  # sensible upper cap) and drop axis ticks when too many clusters remain.
  if (file.exists(ov_csv)) {
    ov <- utils::read.csv(ov_csv, stringsAsFactors = FALSE)
    if (nrow(ov) > 0) {
      ov_sig <- ov[!is.na(ov$fdr) & ov$fdr < 0.1, , drop = FALSE]
      if (nrow(ov_sig) > 0) {
        ov_plot <- ov_sig[order(-abs(ov_sig$log2OR)), , drop = FALSE]
        ov_plot <- utils::head(ov_plot, 600L)
        # Order rows by cluster's max |log2OR| so the strongest pairs sit
        # near the top, and columns the same way.
        g_ord <- names(sort(tapply(abs(ov_plot$log2OR),
                                   ov_plot$gliph_cluster_id, max),
                            decreasing = TRUE))
        c_ord <- names(sort(tapply(abs(ov_plot$log2OR),
                                   ov_plot$clustcr_cluster, max),
                            decreasing = TRUE))
        ov_plot$gliph_cluster_id <- factor(ov_plot$gliph_cluster_id,
                                           levels = g_ord)
        ov_plot$clustcr_cluster  <- factor(ov_plot$clustcr_cluster,
                                           levels = c_ord)
        n_g <- nlevels(ov_plot$gliph_cluster_id)
        n_c <- nlevels(ov_plot$clustcr_cluster)
        drop_x <- n_c > 40L
        drop_y <- n_g > 40L
        p2 <- ggplot2::ggplot(ov_plot,
                              ggplot2::aes(x = clustcr_cluster,
                                           y = gliph_cluster_id,
                                           fill = log2OR)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white",
                                        high = "#B40426", midpoint = 0) +
          ggplot2::labs(title = paste0("GLIPH x clusTCR co-occurrence ",
                                       "(top ", nrow(ov_plot),
                                       " FDR<0.1 pairs by |log2 OR|)"),
                        subtitle = paste0(n_g, " GLIPH x ", n_c,
                                          " clusTCR clusters; axis labels ",
                                          "suppressed above 40 categories"),
                        x = "clusTCR cluster", y = "GLIPH cluster") +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::theme(
            plot.subtitle    = ggplot2::element_text(size = 9),
            axis.text.x      = if (drop_x) ggplot2::element_blank() else
                                 ggplot2::element_text(angle = 45,
                                                       hjust = 1, size = 7),
            axis.text.y      = if (drop_y) ggplot2::element_blank() else
                                 ggplot2::element_text(size = 7),
            axis.ticks.x     = if (drop_x) ggplot2::element_blank() else
                                 ggplot2::element_line(),
            axis.ticks.y     = if (drop_y) ggplot2::element_blank() else
                                 ggplot2::element_line(),
            panel.grid.minor = ggplot2::element_blank())
        # When labels are suppressed, fix the canvas so the file stays
        # readable; when labels are kept (small category counts), grow
        # the canvas modestly to accommodate them.
        w_out <- if (drop_x) 10 else min(20, max(9, n_c * 0.18 + 4))
        h_out <- if (drop_y) 10 else min(20, max(7, n_g * 0.18 + 4))
        save_pdf_png(p2, file.path(p$base, "gliph_clustcr_log2OR_heatmap"),
                     w = w_out, h = h_out)
      } else {
        log_message("  gliph_clustcr_log2OR_heatmap: no FDR<0.1 pairs; ",
                    "skipping plot.")
      }
    }
  }

  # --- Within vs between tcrdist scatter per cluster ----------------------
  d <- conc[!is.na(conc$mean_within_tcrdist), , drop = FALSE]
  if (nrow(d) > 0) {
    d$sig <- ifelse(!is.na(d$fdr) & d$fdr < 0.1, "FDR<0.1", "n.s.")
    p3 <- ggplot2::ggplot(d, ggplot2::aes(mean_between_tcrdist,
                                          mean_within_tcrdist,
                                          color = sig, size = n_cells)) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::scale_color_manual(values = c("FDR<0.1" = "#B40426",
                                             "n.s."   = "grey60")) +
      ggplot2::labs(title = "Within- vs between-cluster tcrdist per GLIPH group",
                    x = "Mean tcrdist, background (Subject-stratified)",
                    y = "Mean tcrdist, within cluster",
                    color = "perm test", size = "n cells") +
      ggplot2::theme_bw(base_size = 11)
    save_pdf_png(p3, file.path(p$base, "gliph_within_vs_between_tcrdist"),
                 w = 7, h = 6)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Module 2: GLIPH network metrics
# ---------------------------------------------------------------------------
viz_gliph_network_metrics <- function(cfg) {
  p <- .tcra_paths(cfg)
  metrics_csv <- file.path(p$tables_rep, "gliph_network_metrics.csv")
  nodes_csv   <- file.path(p$tables_rep, "gliph_network_node_attrs.csv")
  if (!file.exists(metrics_csv) || !file.exists(nodes_csv)) {
    log_message("  viz_gliph_network_metrics: tables missing; skipping.")
    return(invisible(NULL))
  }

  m <- utils::read.csv(metrics_csv, stringsAsFactors = FALSE)
  n <- utils::read.csv(nodes_csv,   stringsAsFactors = FALSE)

  # --- Metric heatmap: cluster x metric, faceted by edge_type --------------
  metric_cols <- c("density", "transitivity", "modularity",
                   "assortativity_phenotype")
  long <- m |>
    dplyr::select(cluster_id, edge_type, dplyr::all_of(metric_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(metric_cols),
                        names_to = "metric", values_to = "value")
  if (nrow(long) > 0) {
    long$cluster_id <- factor(long$cluster_id,
                              levels = unique(long$cluster_id))
    p1 <- ggplot2::ggplot(long,
                          ggplot2::aes(metric, cluster_id, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(option = "viridis", na.value = "grey92") +
      ggplot2::facet_wrap(~ edge_type, nrow = 1) +
      ggplot2::labs(title = "GLIPH network metrics per cluster",
                    x = NULL, y = "GLIPH cluster") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         hjust = 1))
    save_pdf_png(p1, file.path(p$base, "gliph_network_metric_heatmap"),
                 w = 9, h = 7)
  }

  # --- Centrality vs log10_Pgen per edge_type -------------------------------
  if ("log10_Pgen" %in% colnames(n)) {
    nd <- n[is.finite(n$log10_Pgen) & is.finite(n$eigen_centrality), ,
            drop = FALSE]
    if (nrow(nd) > 0) {
      p2 <- ggplot2::ggplot(nd,
                            ggplot2::aes(log10_Pgen, eigen_centrality,
                                         color = dominant_phenotype)) +
        ggplot2::geom_point(alpha = 0.6, size = 1.4) +
        ggplot2::scale_color_manual(values = ETIOLOGY_GROUP_COLORS,
                                    na.value = "grey60") +
        ggplot2::facet_wrap(~ edge_type, nrow = 1) +
        ggplot2::labs(title = "Hub centrality vs generation probability",
                      x = "log10(P_gen)", y = "Eigenvector centrality",
                      color = "Phenotype") +
        ggplot2::theme_bw(base_size = 11)
      save_pdf_png(p2, file.path(p$base, "gliph_hub_eigen_vs_pgen"),
                   w = 9, h = 5)
    }
  }

  # --- Per-cluster ggraph layouts (top N clusters by n_nodes) --------------
  top_n <- 8L
  top_clusters <- m |>
    dplyr::filter(edge_type == "binary") |>
    dplyr::arrange(dplyr::desc(n_nodes)) |>
    utils::head(top_n) |>
    dplyr::pull(cluster_id)

  for (cl in top_clusters) {
    nodes_cl <- n[n$cluster_id == cl & n$edge_type == "tcrdist", , drop = FALSE]
    if (nrow(nodes_cl) < 3) next
    # Re-derive a simple edge list from membership for layout: clique between
    # nodes that appear in the same tcrdist row (cheap; layout-only)
    el <- t(utils::combn(nodes_cl$CDR3b, 2))
    g  <- tryCatch(igraph::graph_from_edgelist(el, directed = FALSE),
                   error = function(e) NULL)
    if (is.null(g)) next
    tg <- tidygraph::as_tbl_graph(g) |>
      tidygraph::activate(nodes) |>
      dplyr::mutate(degree = igraph::degree(g),
                    phen   = nodes_cl$dominant_phenotype[
                      match(name, nodes_cl$CDR3b)])

    pg <- ggraph::ggraph(tg, layout = "fr") +
      ggraph::geom_edge_link(alpha = 0.15) +
      ggraph::geom_node_point(ggplot2::aes(size = degree, color = phen)) +
      ggplot2::scale_color_manual(values = ETIOLOGY_GROUP_COLORS,
                                  na.value = "grey60") +
      ggplot2::labs(title = paste0("GLIPH cluster ", cl,
                                   " (n=", nrow(nodes_cl), " CDR3)"),
                    color = "Phenotype") +
      ggplot2::theme_void(base_size = 10) +
      ggplot2::theme(legend.position = "right")
    save_pdf_png(pg,
                 file.path(p$base, paste0("gliph_network_", cl)),
                 w = 6, h = 5)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Module 3: VDJdb overlay
# ---------------------------------------------------------------------------
viz_vdjdb_overlay <- function(cfg) {
  p <- .tcra_paths(cfg)
  ann_csv <- file.path(p$tables_rep, "vdjdb_annotations.csv")
  enr_csv <- file.path(p$tables_rep, "vdjdb_substate_enrichment.csv")
  ov_csv  <- file.path(p$tables_rep, "vdjdb_gliph_overlap.csv")
  if (!file.exists(ann_csv)) {
    log_message("  viz_vdjdb_overlay: vdjdb_annotations.csv missing; skipping.")
    return(invisible(NULL))
  }
  ann <- utils::read.csv(ann_csv, stringsAsFactors = FALSE)

  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  if (!file.exists(tcell_rds)) {
    log_message("  T cell object missing; skipping UMAP overlays.")
  } else {
    obj <- readRDS(tcell_rds)
    top_species <- names(sort(table(ann$antigen_species),
                              decreasing = TRUE))[1:5]
    top_species <- top_species[!is.na(top_species)]
    enr <- if (file.exists(enr_csv))
             utils::read.csv(enr_csv, stringsAsFactors = FALSE) else NULL
    umap_sz <- cfg$visualization$umap_size %||% 8
    for (sp in top_species) {
      hit_bcs <- ann$barcode[ann$antigen_species == sp]
      mask <- colnames(obj) %in% hit_bcs
      if (sum(mask) == 0) next
      obj_sp <- .stamp_highlight_col(obj, "vdjdb_hit", mask,
                                     label_yes = paste0("VDJdb: ", sp))
      # Pull species color from the shared PATHOGEN_COLORS palette so
      # every per-pathogen UMAP, fraction bar, and motif logo uses the
      # same hex for the same species.
      sp_col <- if (sp %in% names(PATHOGEN_COLORS)) PATHOGEN_COLORS[[sp]] else
                  unname(PATHOGEN_COLORS["Other"])
      pal <- c(Other = "grey85")
      pal[paste0("VDJdb: ", sp)] <- sp_col
      subtitle <- ""
      if (!is.null(enr)) {
        # Prefer CD8_only rows when present (the hypothesis-aligned
        # contrast); fall back to all_T. Older enrichment tables lack
        # the `subset` column, so handle both shapes.
        if ("subset" %in% colnames(enr)) {
          sub_row <- enr[enr$antigen_species == sp &
                         enr$subset == "CD8_only", , drop = FALSE]
          if (nrow(sub_row) == 0)
            sub_row <- enr[enr$antigen_species == sp &
                           enr$subset == "all_T", , drop = FALSE]
        } else {
          sub_row <- enr[enr$antigen_species == sp, , drop = FALSE]
        }
        if (nrow(sub_row) > 0) {
          best <- sub_row[which.min(sub_row$fdr)[1], ]
          subtitle <- sprintf("n=%d  top substate=%s  FDR=%.2g",
                              sum(mask), best$substate_key, best$fdr)
        } else {
          subtitle <- sprintf("n=%d cells", sum(mask))
        }
      } else {
        subtitle <- sprintf("n=%d cells", sum(mask))
      }
      base <- file.path(p$base,
                        paste0("umap_vdjdb_", make.names(sp)))
      title_txt <- sprintf("VDJdb antigen species: %s", sp)
      dual_save_plot(bquote(
        scplotter::CellDimPlot(.(obj_sp),
                               group_by = "vdjdb_hit",
                               reduction = "UMAP",
                               highlight = .(sprintf('vdjdb_hit == "VDJdb: %s"', sp)),
                               palcolor = .(pal),
                               bg_color = "grey92",
                               pt_alpha = 0.5,
                               highlight_size = 2.0,
                               highlight_stroke = 0.5,
                               raster = FALSE,
                               show_stat = FALSE) +
          ggplot2::ggtitle(.(title_txt), subtitle = .(subtitle))
      ), base, width = umap_sz, height = umap_sz)
    }
    rm(obj); invisible(gc())
  }

  # --- GLIPH x VDJdb antigen-species fraction heatmap ---------------------
  # Top-N filter: ~580 GLIPH clusters overwhelms a single bar plot, so keep
  # the most informative ones (highest n_annotated, breaking ties by
  # fraction). This focuses the figure on clusters with both broad VDJdb
  # coverage and strong species enrichment.
  if (file.exists(ov_csv)) {
    ov <- utils::read.csv(ov_csv, stringsAsFactors = FALSE)
    if (nrow(ov) > 0) {
      ov$sig <- ifelse(!is.na(ov$fdr) & ov$fdr < 0.1, "*", "")
      top_n <- min(40L, nrow(ov))
      ov_top <- ov[order(-ov$n_annotated, -ov$fraction_annotated), ,
                   drop = FALSE]
      ov_top <- utils::head(ov_top, top_n)
      # Pull palette from the shared PATHOGEN_COLORS map; species not in
      # the palette fall through to "Other".
      sp_levels <- sort(unique(ov_top$dominant_antigen_species))
      sp_pal <- vapply(sp_levels, function(s) {
        if (s %in% names(PATHOGEN_COLORS)) PATHOGEN_COLORS[[s]] else
          unname(PATHOGEN_COLORS["Other"])
      }, character(1))
      names(sp_pal) <- sp_levels
      p1 <- ggplot2::ggplot(ov_top,
                            ggplot2::aes(reorder(cluster_id,
                                                 fraction_annotated),
                                         fraction_annotated,
                                         fill = dominant_antigen_species)) +
        ggplot2::geom_col() +
        ggplot2::geom_text(ggplot2::aes(label = sig),
                           hjust = -0.2, size = 3.5) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = sp_pal, na.value = "grey85") +
        ggplot2::labs(title  = paste0("GLIPH cluster annotation by VDJdb ",
                                      "(top ", top_n,
                                      " by n_annotated; * FDR<0.1)"),
                      x = NULL, y = "Fraction of cluster CDR3 annotated",
                      fill = "Dominant antigen species") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                       axis.text.x = ggplot2::element_text(size = 9),
                       legend.position = "right")
      save_pdf_png(p1, file.path(p$base, "gliph_vdjdb_fraction_bar"),
                   w = 9, h = max(6, top_n * 0.22 + 2))
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Module 3b: Antigen-class substate composition bars
# Stacked bars with x = antigen class (VDJdb species n>=10, plus HLA-B27
# pathogenic motif) and fill = T cell substate, position="fill" so the
# y-axis reads as relative cluster composition. A cell that is e.g. both
# CMV+ via VDJdb and HLA-B27 pathogenic by TRAV21+CDR3b motif contributes
# to *both* bars; the n shown above each bar is the cell count for that
# antigen class. Two panels written:
#   antigen_substate_composition_overall          (all phenotypes pooled)
#   antigen_substate_composition_by_phenotype     (facet_wrap NIU vs Viral)
# Side-effect CSV antigen_substate_composition_counts.csv carries the
# (antigen_class, Phenotype_2, substate, n) triples used to build the bars.
# ---------------------------------------------------------------------------
viz_antigen_substate_composition <- function(cfg) {
  p <- .tcra_paths(cfg)
  ann_csv <- file.path(p$tables_rep, "vdjdb_annotations.csv")
  if (!file.exists(ann_csv)) {
    log_message("  viz_antigen_substate_composition: vdjdb_annotations.csv ",
                "missing; skipping.")
    return(invisible(NULL))
  }
  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  if (!file.exists(tcell_rds)) {
    log_message("  viz_antigen_substate_composition: T cell object missing; ",
                "skipping.")
    return(invisible(NULL))
  }

  obj <- readRDS(tcell_rds)
  meta <- obj@meta.data
  meta$barcode  <- colnames(obj)
  meta$substate <- substate_labels(cfg, "tcell", meta$knn.leiden.cluster)
  rm(obj); invisible(gc())

  # VDJdb species cells: one (barcode, species) row per distinct match,
  # restricted to species with >= 10 unique T cell barcodes so the panel
  # does not pick up singleton hits.
  # Note (2026-05-22): HLA-B27 pathogenic motif removed from this panel —
  # the first set of figures focuses on viral specificities only. The
  # HLA-B27 arm is carried by its own panels (Fig 5C per-subject pathogenic
  # fraction, Fig 5H per-subject GLIPH cluster recovery).
  vdjdb <- utils::read.csv(ann_csv, stringsAsFactors = FALSE)
  vdjdb <- dplyr::distinct(vdjdb, .data$barcode, .data$antigen_species)
  vdjdb <- vdjdb[vdjdb$barcode %in% meta$barcode, , drop = FALSE]
  sp_n    <- table(vdjdb$antigen_species)
  keep_sp <- names(sp_n)[sp_n >= 10]
  vdjdb   <- vdjdb[vdjdb$antigen_species %in% keep_sp, , drop = FALSE]
  long  <- vdjdb |>
    dplyr::left_join(meta[, c("barcode", "substate", "Phenotype_2")],
                     by = "barcode") |>
    dplyr::rename(antigen_class = "antigen_species")

  if (nrow(long) == 0) {
    log_message("  viz_antigen_substate_composition: no annotated cells; ",
                "skipping.")
    return(invisible(NULL))
  }

  cls_levels <- names(sort(table(long$antigen_class), decreasing = TRUE))
  long$antigen_class <- factor(long$antigen_class, levels = cls_levels)
  sub_levels <- sort(unique(long$substate))
  long$substate <- factor(long$substate, levels = sub_levels)

  # Match the project's Paired-palette convention (mirrors paired_palette()
  # in R/82_viz_dispatch.R:1268). brewer.pal requires n>=3, so we always
  # pull all 12 Paired slots, take the first n we need, and ramp only when
  # the substate count exceeds 12.
  paired_base <- RColorBrewer::brewer.pal(12, "Paired")
  pal_subs <- if (length(sub_levels) <= 12) {
    paired_base[seq_along(sub_levels)]
  } else {
    grDevices::colorRampPalette(paired_base)(length(sub_levels))
  }
  names(pal_subs) <- sub_levels

  n_overall <- dplyr::count(long, .data$antigen_class, name = "n")
  n_pheno   <- dplyr::count(
    dplyr::filter(long, .data$Phenotype_2 %in% c("NIU", "Viral")),
    .data$Phenotype_2, .data$antigen_class, name = "n")

  out_counts <- dplyr::bind_rows(
    dplyr::count(long, .data$antigen_class, .data$substate, name = "n") |>
      dplyr::mutate(Phenotype_2 = "ALL"),
    dplyr::count(long, .data$Phenotype_2, .data$antigen_class,
                 .data$substate, name = "n")
  )
  utils::write.csv(out_counts,
                   file.path(p$tables_rep,
                             "antigen_substate_composition_counts.csv"),
                   row.names = FALSE)

  pA <- ggplot2::ggplot(long,
                        ggplot2::aes(x = .data$antigen_class,
                                     fill = .data$substate)) +
    ggplot2::geom_bar(position = "fill", width = 0.78,
                      color = "white", linewidth = 0.2) +
    ggplot2::geom_text(data = n_overall,
                       ggplot2::aes(x = .data$antigen_class, y = 1.02,
                                    label = paste0("n=", .data$n)),
                       inherit.aes = FALSE, vjust = 0, size = 3) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(add = c(0, 0.08))) +
    ggplot2::scale_fill_manual(values = pal_subs,
                               name = "T cell substate") +
    ggplot2::labs(title = "T cell substate composition by antigen class",
                  subtitle = paste0("VDJdb species (n>=10) + HLA-B27 ",
                                    "pathogenic motif. Cells appear in ",
                                    "every class they match."),
                  x = NULL, y = "Relative frequency") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title  = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 30,
                                                       hjust = 1))
  save_pdf_png(pA,
               file.path(p$base, "antigen_substate_composition_overall"),
               w = max(6, length(cls_levels) * 0.6 + 3), h = 5.5)

  long_pheno <- dplyr::filter(long,
                              .data$Phenotype_2 %in% c("NIU", "Viral"))
  if (nrow(long_pheno) > 0) {
    # Display labels: NIU -> "Autoimmune" (matches the published Fig 2/3/4
    # phenotype labels and the screenshot reference for this panel).
    long_pheno$pheno_lbl <- factor(
      dplyr::recode(as.character(long_pheno$Phenotype_2),
                    NIU = "Autoimmune", Viral = "Viral"),
      levels = c("Autoimmune", "Viral"))

    # Panel B palette: RColorBrewer "Paired" (project convention). Ramp
    # when substate count exceeds 12.
    paired_base <- RColorBrewer::brewer.pal(12, "Paired")
    pal_subs_pB <- if (length(sub_levels) <= 12) {
      paired_base[seq_along(sub_levels)]
    } else {
      grDevices::colorRampPalette(paired_base)(length(sub_levels))
    }
    names(pal_subs_pB) <- sub_levels

    # Pad missing (Phenotype, antigen_class) combinations with a grey "No
    # data" sentinel row so the bar renders as a solid grey strip when an
    # antigen class has zero cells in one phenotype (e.g., HTLV-1 in
    # Autoimmune).
    full_grid <- expand.grid(pheno_lbl     = c("Autoimmune", "Viral"),
                             antigen_class = cls_levels,
                             stringsAsFactors = FALSE)
    have <- dplyr::distinct(long_pheno, .data$pheno_lbl, .data$antigen_class)
    have$present <- TRUE
    missing_grid <- dplyr::anti_join(full_grid, have,
                                     by = c("pheno_lbl", "antigen_class"))
    if (nrow(missing_grid) > 0L) {
      missing_grid$substate <- "No data"
      missing_grid$barcode  <- NA_character_
      missing_grid$Phenotype_2 <- ifelse(missing_grid$pheno_lbl ==
                                          "Autoimmune", "NIU", "Viral")
      missing_grid <- missing_grid[, c("barcode", "substate", "Phenotype_2",
                                       "antigen_class", "pheno_lbl")]
      long_pheno_aug <- dplyr::bind_rows(
        long_pheno[, c("barcode", "substate", "Phenotype_2",
                       "antigen_class", "pheno_lbl")],
        missing_grid)
    } else {
      long_pheno_aug <- long_pheno
    }
    long_pheno_aug$substate <- factor(long_pheno_aug$substate,
                                      levels = c(sub_levels, "No data"))
    # Reverse antigen factor so the top of the y-axis is the highest-count
    # species (matches the screenshot's CMV-on-top ordering).
    long_pheno_aug$antigen_class <- factor(long_pheno_aug$antigen_class,
                                           levels = rev(cls_levels))
    long_pheno_aug$pheno_lbl <- factor(long_pheno_aug$pheno_lbl,
                                       levels = c("Autoimmune", "Viral"))
    pal_subs_pB_full <- c(pal_subs_pB, `No data` = "grey75")

    # n_pheno labels keyed by the display label too.
    n_pheno_lbl <- n_pheno |>
      dplyr::mutate(pheno_lbl = factor(
        dplyr::recode(as.character(.data$Phenotype_2),
                      NIU = "Autoimmune", Viral = "Viral"),
        levels = c("Autoimmune", "Viral")))
    n_pheno_lbl$antigen_class <- factor(n_pheno_lbl$antigen_class,
                                        levels = rev(cls_levels))

    pB <- ggplot2::ggplot(long_pheno_aug,
                          ggplot2::aes(y = .data$antigen_class,
                                       fill = .data$substate)) +
      ggplot2::geom_bar(position = "fill", width = 0.78,
                        color = "white", linewidth = 0.2) +
      ggplot2::geom_text(data = n_pheno_lbl,
                         ggplot2::aes(y = .data$antigen_class, x = 1.02,
                                      label = paste0("n=", .data$n)),
                         inherit.aes = FALSE, hjust = 0, size = 3) +
      ggh4x::facet_grid2(rows = ggplot2::vars(.data$pheno_lbl),
                         scales = "free_y", space = "free_y",
                         strip = ggh4x::strip_themed(
                           background_y = ggh4x::elem_list_rect(fill = NA),
                           text_y = ggh4x::elem_list_text(angle = -90,
                                                          face = "bold"))) +
      ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::scale_x_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = ggplot2::expansion(add = c(0, 0.12))) +
      ggplot2::scale_fill_manual(values = pal_subs_pB_full,
                                 name = "T cell substate",
                                 drop = FALSE) +
      ggplot2::labs(x = "Relative frequency", y = NULL) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(strip.background.y = ggplot2::element_rect(fill = NA,
                                                                color = NA),
                     strip.text.y.right = ggplot2::element_text(angle = -90,
                                                                face = "bold"),
                     axis.text.y       = ggplot2::element_text(size = 10),
                     legend.position   = "right")
    save_pdf_png(pB,
                 file.path(p$base,
                           "antigen_substate_composition_by_phenotype"),
                 w = 8, h = max(4, length(cls_levels) * 0.45 + 2))
  }
  log_message("  viz_antigen_substate_composition: wrote overall + ",
              "by-phenotype panels (n classes = ", length(cls_levels), ").")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Module 3c: Single stacked bar — substate composition of HLA-B27 pathogenic
# motif T cells (TRAV21 + [YF]S[TS] CDR3-beta), eye tissue only. Companion
# to viz_antigen_substate_composition but restricted to one antigen class
# (HLA-B27 pathogenic motif) so there is exactly one x-category. Substates
# stack along the y-axis (position="fill"), n is annotated above the bar.
# Side-effect CSV hla_b27_substate_composition_counts.csv carries
# (substate, n) for the bar.
# ---------------------------------------------------------------------------
viz_hla_b27_substate_composition <- function(cfg) {
  p <- .tcra_paths(cfg)
  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  if (!file.exists(tcell_rds)) {
    log_message("  viz_hla_b27_substate_composition: T cell object missing; ",
                "skipping.")
    return(invisible(NULL))
  }
  obj <- readRDS(tcell_rds)
  meta <- obj@meta.data
  meta$barcode  <- colnames(obj)
  meta$substate <- substate_labels(cfg, "tcell", meta$knn.leiden.cluster)
  meta$pathogenic <- flag_hla_b27_pathogenic_tcr(meta$CTgene, meta$CTaa)
  rm(obj); invisible(gc())

  if ("Tissue_1" %in% colnames(meta)) {
    meta <- meta[meta$Tissue_1 == "Eye", , drop = FALSE]
  }
  meta <- meta[meta$pathogenic, , drop = FALSE]
  if (!nrow(meta)) {
    log_message("  viz_hla_b27_substate_composition: no eye cells carry ",
                "the HLA-B27 pathogenic motif; skipping.")
    return(invisible(NULL))
  }

  sub_levels <- sort(unique(meta$substate))
  meta$substate <- factor(meta$substate, levels = sub_levels)

  paired_base <- RColorBrewer::brewer.pal(12, "Paired")
  pal_subs <- if (length(sub_levels) <= 12) {
    paired_base[seq_along(sub_levels)]
  } else {
    grDevices::colorRampPalette(paired_base)(length(sub_levels))
  }
  names(pal_subs) <- sub_levels

  counts <- dplyr::count(meta, .data$substate, name = "n")
  utils::write.csv(counts,
                   file.path(p$tables_rep,
                             "hla_b27_substate_composition_counts.csv"),
                   row.names = FALSE)

  n_total <- nrow(meta)
  plot_df <- data.frame(x = "HLA-B27", meta[, "substate", drop = FALSE])

  pH <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = .data$x, fill = .data$substate)) +
    ggplot2::geom_bar(position = "fill", width = 0.55,
                      color = "white", linewidth = 0.2) +
    ggplot2::annotate("text", x = 1, y = 1.02,
                      label = paste0("n=", n_total),
                      vjust = 0, size = 3.2) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(add = c(0, 0.08))) +
    ggplot2::scale_fill_manual(values = pal_subs,
                               name = "T cell substate") +
    ggplot2::labs(
      title    = paste0("T cell substate composition of HLA-B27 ",
                        "pathogenic TCR cells"),
      subtitle = paste0("Intraocular T cells carrying TRAV21 + [YF]S[TS] ",
                        "CDR3-beta motif."),
      x = NULL, y = "Relative frequency") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title  = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 11))
  save_pdf_png(pH,
               file.path(p$base, "hla_b27_substate_composition"),
               w = 4.2, h = 5.5)
  log_message("  viz_hla_b27_substate_composition: wrote single-bar panel ",
              "(n cells = ", n_total, ", n substates = ",
              length(sub_levels), ").")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Module 4: TCR-genex signatures
# ---------------------------------------------------------------------------
viz_tcr_genex_signatures <- function(cfg) {
  p <- .tcra_paths(cfg)
  sig_csv  <- file.path(p$tables_t, "tcr_signature_scores_by_clone_group.csv")
  cpl_csv  <- file.path(p$tables_t, "tcr_myeloid_coupling.csv")

  if (file.exists(sig_csv)) {
    sig <- utils::read.csv(sig_csv, stringsAsFactors = FALSE)
    if (nrow(sig) > 0) {
      sig$signed_log10 <- -log10(pmax(sig$p, 1e-300)) * sign(sig$delta)
      p1 <- ggplot2::ggplot(sig,
                            ggplot2::aes(panel, substate_key,
                                         fill = signed_log10)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white",
                                      high = "#B40426", midpoint = 0) +
        ggplot2::facet_wrap(~ clone_group, ncol = 1) +
        ggplot2::labs(title = "UCell signatures: in_group vs out_group per substate",
                      fill = "signed -log10 p", x = "Signature panel",
                      y = "T cell substate") +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                           hjust = 1))
      save_pdf_png(p1, file.path(p$base, "tcr_signature_heatmap"),
                   w = 8, h = 9)
    }
  }

  # Volcanoes per clone group, reading the pseudobulk DE tables.
  for (g in c("expanded", "public", "gliph")) {
    de_csv <- file.path(p$tables_t,
                        paste0("tcr_pseudobulk_DE_", g, ".csv"))
    if (!file.exists(de_csv)) next
    de <- utils::read.csv(de_csv, stringsAsFactors = FALSE)
    if (nrow(de) == 0) next
    de$sig <- ifelse(!is.na(de$padj) & de$padj < 0.05 & abs(de$log2FoldChange) > 0.5,
                     "sig", "ns")
    pv <- ggplot2::ggplot(de,
                          ggplot2::aes(log2FoldChange, -log10(pmax(padj, 1e-300)),
                                       color = sig)) +
      ggplot2::geom_point(alpha = 0.5, size = 0.7) +
      ggplot2::scale_color_manual(values = c(sig = "#B40426", ns = "grey75")) +
      ggplot2::labs(title = paste0("Pseudobulk DE: ", g, " in_group vs out_group"),
                    x = "log2 fold change", y = "-log10 padj") +
      ggplot2::theme_bw(base_size = 10)
    save_pdf_png(pv,
                 file.path(p$base, paste0("tcr_volcano_", g)),
                 w = 6.5, h = 5)
  }

  # Coupling diagram (bipartite ggplot via geom_segment).
  if (file.exists(cpl_csv)) {
    cpl <- utils::read.csv(cpl_csv, stringsAsFactors = FALSE)
    if (nrow(cpl) > 0) {
      tcell_states   <- sort(unique(cpl$tcell_substate))
      myeloid_states <- sort(unique(cpl$myeloid_substate))
      yT <- setNames(seq_along(tcell_states),   tcell_states)
      yM <- setNames(seq_along(myeloid_states), myeloid_states)
      cpl$x_left  <- 0
      cpl$x_right <- 1
      cpl$y_left  <- yT[cpl$tcell_substate]
      cpl$y_right <- yM[cpl$myeloid_substate]
      pc <- ggplot2::ggplot(cpl) +
        ggplot2::geom_segment(ggplot2::aes(x = x_left, xend = x_right,
                                           y = y_left, yend = y_right,
                                           linewidth = coupling_strength,
                                           alpha = coupling_strength),
                              color = "#7A1A30") +
        ggplot2::geom_text(data = data.frame(state = names(yT), y = yT),
                           ggplot2::aes(x = -0.05, y = y, label = state),
                           hjust = 1, size = 3) +
        ggplot2::geom_text(data = data.frame(state = names(yM), y = yM),
                           ggplot2::aes(x = 1.05, y = y, label = state),
                           hjust = 0, size = 3) +
        ggplot2::scale_linewidth(range = c(0.3, 2.5), guide = "none") +
        ggplot2::scale_alpha(range = c(0.3, 0.9), guide = "none") +
        ggplot2::scale_x_continuous(limits = c(-0.6, 1.6)) +
        ggplot2::labs(title = "T cell <-> myeloid coupling (TCR-driven substates)",
                      x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = 11)
      save_pdf_png(pc, file.path(p$base, "tcr_myeloid_coupling_diagram"),
                   w = 8, h = 7)
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Module 5: meta-motif assembly figure (VGKWY-style)
# ---------------------------------------------------------------------------
viz_meta_motifs <- function(cfg) {
  p <- .tcra_paths(cfg)
  asm_rds <- file.path(p$objs, "GliphMotifAssemblies.rds")
  sum_csv <- file.path(p$tables_rep, "meta_motif_assemblies.csv")
  if (!file.exists(asm_rds) || !file.exists(sum_csv)) {
    log_message("  viz_meta_motifs: assembly outputs missing; skipping.")
    return(invisible(NULL))
  }
  asm_data <- readRDS(asm_rds)
  summary_df <- utils::read.csv(sum_csv, stringsAsFactors = FALSE)
  if (length(asm_data$assemblies) == 0) return(invisible(NULL))

  # Build the per-meta-motif tile figure: rows = constituent motifs aligned
  # at their offsets, columns = position. Consensus shown above as text.
  for (a in asm_data$assemblies) {
    consts  <- a$constituent_motifs
    offsets <- a$position_offsets
    if (length(consts) < 2) next
    max_end <- max(offsets + nchar(consts))
    tiles <- do.call(rbind, lapply(seq_along(consts), function(i) {
      cc <- strsplit(consts[i], "")[[1]]
      data.frame(motif = consts[i],
                 row = i,
                 col = seq_along(cc) + offsets[i],
                 aa  = cc,
                 stringsAsFactors = FALSE)
    }))
    cons_chars <- strsplit(a$consensus, "")[[1]]
    cons_df <- data.frame(col = seq_along(cons_chars), aa = cons_chars)
    pal <- c(A = "#7E91C3", C = "#FFD040", D = "#D03B3B", E = "#D03B3B",
             F = "#7E5BBA", G = "#9CCB9C", H = "#7CC8C5", I = "#9CCB9C",
             K = "#7E91C3", L = "#9CCB9C", M = "#9CCB9C", N = "#FFB575",
             P = "#FFD040", Q = "#FFB575", R = "#7E91C3", S = "#FFB575",
             T = "#FFB575", V = "#9CCB9C", W = "#7E5BBA", Y = "#7E5BBA")
    p1 <- ggplot2::ggplot(tiles,
                          ggplot2::aes(col, factor(row), fill = aa, label = aa)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(size = 3) +
      ggplot2::geom_label(data = cons_df,
                          ggplot2::aes(col, length(consts) + 1, label = aa),
                          inherit.aes = FALSE,
                          fill = "#C02C2C", color = "white",
                          label.padding = ggplot2::unit(0.1, "lines"),
                          size = 3.2) +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::scale_y_discrete(breaks = seq_along(consts),
                                labels = consts) +
      ggplot2::scale_x_continuous(breaks = seq_len(max_end)) +
      ggplot2::labs(title = paste0(a$meta_motif_id, " consensus: ", a$consensus),
                    x = "Position", y = NULL) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    save_pdf_png(p1,
                 file.path(p$base,
                           paste0("meta_motif_", a$meta_motif_id)),
                 w = 7, h = 1.5 + 0.5 * length(consts))
  }

  # Companion UMAP overlays per meta-motif on the T cell compartment object
  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  if (file.exists(tcell_rds)) {
    obj <- readRDS(tcell_rds)
    trb_aa <- .cell_trb_cdr3(obj)
    umap_sz <- cfg$visualization$umap_size %||% 8
    for (a in asm_data$assemblies) {
      hits <- vapply(trb_aa, function(s)
        any(vapply(a$constituent_motifs,
                   function(m) grepl(m, s, fixed = TRUE),
                   logical(1))),
        logical(1))
      hits[is.na(hits)] <- FALSE
      if (sum(hits) == 0) next
      obj_h <- .stamp_highlight_col(obj, "meta_motif_hit", hits,
                                    label_yes = a$meta_motif_id)
      pal <- c(Other = "grey85")
      pal[a$meta_motif_id] <- "#7A1A30"
      title <- sprintf("%s (n=%d cells, consensus %s)",
                       a$meta_motif_id, sum(hits), a$consensus)
      base <- file.path(p$base, paste0("umap_meta_", a$meta_motif_id))
      dual_save_plot(bquote(
        scplotter::CellDimPlot(.(obj_h),
                               group_by = "meta_motif_hit",
                               reduction = "UMAP",
                               highlight = .(sprintf('meta_motif_hit == "%s"',
                                                     a$meta_motif_id)),
                               palcolor = .(pal),
                               bg_color = "grey92",
                               pt_alpha = 0.5,
                               highlight_size = 2.2,
                               highlight_stroke = 0.5,
                               raster = FALSE,
                               show_stat = FALSE) +
          ggplot2::ggtitle(.(title))
      ), base, width = umap_sz, height = umap_sz)
    }
    rm(obj); invisible(gc())
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Module 5b: Directional GLIPH motifs
#   Two outputs that work even when the OLC meta-motif assembler in R/71
#   returns empty (which it does here because eye-cohort clones rarely
#   meet the shared-expanded/shared-public weight floor):
#     - gliph_motif_logos_by_direction: per-cluster CDR3 sequence logos
#       for the top-N Viral- and NIU-enriched GLIPH clusters, length-matched
#       within each cluster. Colored by the shared PATHOGEN_COLORS palette
#       where the cluster's dominant VDJdb annotation matches a known
#       pathogen, otherwise grey.
#     - gliph_motif_overlap_network: nodes = motifs, edges = pairwise
#       overlaps (shared CDR3 count), node fill = NIU/Viral/none direction.
#       The R/71 OLC weights are too strict to assemble meta-motifs, but
#       the overlap edges themselves still reveal structure.
# ---------------------------------------------------------------------------
viz_gliph_motif_directional <- function(cfg, top_n = 12L) {
  p <- .tcra_paths(cfg)
  enr_csv  <- file.path(p$tables_rep, "gliph_enrichment_viral_vs_niu.csv")
  cls_csv  <- file.path(p$tables_rep, "gliph_clusters.csv")
  edge_csv <- file.path(p$tables_rep, "motif_overlap_edges.csv")
  ov_csv   <- file.path(p$tables_rep, "vdjdb_gliph_overlap.csv")
  if (!file.exists(enr_csv) || !file.exists(cls_csv)) {
    log_message("  viz_gliph_motif_directional: required tables missing; ",
                "skipping.")
    return(invisible(NULL))
  }
  enr <- utils::read.csv(enr_csv, stringsAsFactors = FALSE)
  cls <- utils::read.csv(cls_csv, stringsAsFactors = FALSE)
  enr$direction <- ifelse(grepl("^Viral", enr$direction, ignore.case = TRUE),
                          "Viral_enriched",
                   ifelse(grepl("^NIU",   enr$direction, ignore.case = TRUE),
                          "NIU_enriched", enr$direction))

  # ---- 5b.i: per-direction CDR3 sequence logos --------------------------
  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    log_message("  ggseqlogo not installed; skipping motif logos.")
  } else {
    top_viral <- enr[enr$direction == "Viral_enriched" & !is.na(enr$FDR), ,
                     drop = FALSE]
    top_viral <- utils::head(top_viral[order(top_viral$FDR), , drop = FALSE],
                             top_n)
    top_niu   <- enr[enr$direction == "NIU_enriched"   & !is.na(enr$FDR), ,
                     drop = FALSE]
    top_niu   <- utils::head(top_niu[order(top_niu$FDR), , drop = FALSE],
                             top_n)
    picks <- rbind(top_viral, top_niu)

    seq_lists <- list()
    for (i in seq_len(nrow(picks))) {
      cid <- picks$cluster_id[i]
      s   <- cls$CDR3b[cls$cluster_id == cid]
      s   <- s[!is.na(s) & nchar(s) > 0]
      if (length(s) < 3) next
      # Length-match: take the modal length so the logo is a clean PWM.
      L <- as.integer(names(sort(table(nchar(s)), decreasing = TRUE))[1])
      s <- s[nchar(s) == L]
      if (length(s) < 3) next
      side <- if (picks$direction[i] == "Viral_enriched") "Viral" else "NIU"
      lab  <- sprintf("%s | %s (n=%d, FDR=%.2g)",
                      side, picks$motif[i], length(s), picks$FDR[i])
      seq_lists[[lab]] <- s
    }
    if (length(seq_lists) > 0) {
      lp <- ggseqlogo::ggseqlogo(seq_lists, ncol = 3) +
        ggplot2::labs(
          title    = paste0("GLIPH cluster CDR3 logos by NIU vs Viral ",
                            "enrichment (top ", top_n, " per direction)"),
          subtitle = "Length-matched CDR3-beta within each convergence group; FDR from R/64") +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(strip.text   = ggplot2::element_text(size = 7,
                                                            face = "bold"),
                       plot.title    = ggplot2::element_text(face = "bold"),
                       plot.subtitle = ggplot2::element_text(size = 9))
      save_pdf_png(lp, file.path(p$base,
                                 "gliph_motif_logos_by_direction"),
                   w = 13, h = max(7, ceiling(length(seq_lists) / 3) * 1.6))
    } else {
      log_message("  No GLIPH clusters had >=3 same-length CDR3s for logo.")
    }
  }

  # ---- 5b.ii: motif overlap network -------------------------------------
  # Use ANY overlap (shared_cdr3 > 0) rather than R/71's strict
  # shared_expanded + 0.5*shared_public >= 1 floor, since the eye cohort
  # rarely meets that threshold but raw motif overlap still tells the
  # structural story.
  if (file.exists(edge_csv) &&
      requireNamespace("igraph", quietly = TRUE) &&
      requireNamespace("ggraph", quietly = TRUE)) {
    edges <- utils::read.csv(edge_csv, stringsAsFactors = FALSE)
    edges <- edges[!is.na(edges$shared_cdr3) & edges$shared_cdr3 > 0, ,
                   drop = FALSE]
    # Motif -> dominant direction lookup. A motif belongs to >=1 GLIPH
    # cluster; take the direction of the most significant cluster (lowest
    # FDR) using that motif.
    motif_dir <- enr[!is.na(enr$motif) & nchar(enr$motif) > 0, , drop = FALSE]
    motif_dir <- motif_dir[order(motif_dir$FDR), , drop = FALSE]
    motif_dir <- motif_dir[!duplicated(motif_dir$motif), , drop = FALSE]
    # Strip GLIPH wildcards so motif_dir keys match the edges table's
    # cleaned motif strings (R/71 uses .gmo_clean_motif).
    clean <- function(s) gsub("%", "", gsub("\\{[^}]+\\}", "", toupper(s)))
    motif_dir$motif_clean <- clean(motif_dir$motif)
    dir_lookup <- setNames(motif_dir$direction, motif_dir$motif_clean)

    nodes <- unique(c(edges$motif_a, edges$motif_b))
    if (length(nodes) >= 2 && nrow(edges) > 0) {
      node_df <- data.frame(
        name = nodes,
        direction = unname(dir_lookup[nodes]),
        stringsAsFactors = FALSE)
      node_df$direction[is.na(node_df$direction)] <- "Unannotated"
      g <- igraph::graph_from_data_frame(
        d = edges[, c("motif_a", "motif_b", "edge_weight", "shared_cdr3",
                      "overlap_len")],
        directed = TRUE, vertices = node_df)
      # Keep only the largest connected components (top 6 by node count)
      # to keep the plot legible.
      comps <- igraph::components(g, mode = "weak")
      keep_comps <- order(comps$csize, decreasing = TRUE)[
        seq_len(min(6, comps$no))]
      keep_v <- which(comps$membership %in% keep_comps)
      g_sub <- igraph::induced_subgraph(g, vids = keep_v)
      dir_colors <- c(Viral_enriched = unname(ETIOLOGY_GROUP_COLORS["Viral"]),
                      NIU_enriched   = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
                      Unannotated    = "grey75")
      np <- ggraph::ggraph(g_sub, layout = "fr") +
        ggraph::geom_edge_link(ggplot2::aes(width = .data$shared_cdr3,
                                            alpha = .data$shared_cdr3),
                               color = "grey45", arrow = NULL,
                               end_cap = ggraph::circle(3, "mm")) +
        ggraph::geom_node_point(ggplot2::aes(color = .data$direction),
                                size = 6, alpha = 0.92) +
        ggraph::geom_node_text(ggplot2::aes(label = .data$name),
                               size = 3, repel = TRUE,
                               max.overlaps = 50) +
        ggplot2::scale_color_manual(values = dir_colors,
                                    name   = "Direction") +
        ggraph::scale_edge_width(range = c(0.3, 2.5),
                                 name = "Shared CDR3") +
        ggraph::scale_edge_alpha(range = c(0.4, 0.9), guide = "none") +
        ggplot2::labs(
          title    = "GLIPH motif overlap network",
          subtitle = paste0(igraph::vcount(g_sub), " motifs in top-",
                            length(keep_comps),
                            " components; edges = pairwise suffix/prefix ",
                            "overlap weighted by shared CDR3 count")) +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                       plot.subtitle = ggplot2::element_text(size = 9),
                       legend.position = "right")
      save_pdf_png(np, file.path(p$base,
                                 "gliph_motif_overlap_network"),
                   w = 12, h = 10)
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Module 6: Clone-group <-> myeloid->T cell LR coupling
#   For each clone group (expanded / public / gliph), join the in-group vs
#   out-group pseudobulk DE results to the receptors of the condition-specific
#   LR pairs identified in G4b (R/47 liana output). Answers: "do the cells
#   carrying expanded / public / GLIPH-enriched clones differentially express
#   the receptors of the NIU- and Viral-driving LR pairs?" — i.e., are these
#   clone groups participating in the cross-compartment signaling axis the
#   myeloid->T cell L:R analysis surfaced.
#
#   Default scope: T cell DA substates 0, 1, 4, 5 (the substates with
#   differential abundance between conditions). Genes split by direction
#   (Viral-unique LR receptors / NIU-unique LR receptors).
# ---------------------------------------------------------------------------
.cgl_split_subunits <- function(s) {
  s <- as.character(s)
  s <- s[!is.na(s) & nchar(s) > 0]
  unique(unlist(strsplit(s, "[_&+]"), use.names = FALSE))
}

viz_clone_group_lr_coupling <- function(cfg,
                                        tcell_keep = c("0","1","4","5"),
                                        clone_groups = c("expanded","public","gliph"),
                                        padj_cut    = 0.05) {
  p <- .tcra_paths(cfg)
  paths_t  <- get_target_paths(cfg, "tcell")
  cc_paths <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/")
  liana_csv <- file.path(cc_paths$tables, "liana_myeloid_to_tcell_combined.csv")
  if (!file.exists(liana_csv)) {
    log_message("  viz_clone_group_lr_coupling: LIANA table missing; skipping.")
    return(invisible(NULL))
  }
  liana <- utils::read.csv(liana_csv, stringsAsFactors = FALSE)
  # Restrict to myeloid -> tcell rows targeting the DA substates.
  liana <- liana[grepl("^myeloid_", liana$source) &
                 grepl("^tcell_",   liana$target), , drop = FALSE]
  liana$target_id <- sub("^tcell_", "", liana$target)
  liana <- liana[liana$target_id %in% tcell_keep, , drop = FALSE]
  if (!nrow(liana)) {
    log_message("  viz_clone_group_lr_coupling: no LIANA rows after target filter.")
    return(invisible(NULL))
  }
  if (!"unique_to" %in% colnames(liana)) {
    log_message("  viz_clone_group_lr_coupling: unique_to column missing; skipping.")
    return(invisible(NULL))
  }

  receptors <- list(
    Viral = .cgl_split_subunits(
      liana$receptor_complex[liana$unique_to == "Viral"]),
    NIU   = .cgl_split_subunits(
      liana$receptor_complex[liana$unique_to == "NIU"]))
  all_receptors <- unique(c(receptors$Viral, receptors$NIU))
  if (!length(all_receptors)) {
    log_message("  viz_clone_group_lr_coupling: no condition-unique receptors.")
    return(invisible(NULL))
  }
  log_message("  Condition-unique LR receptors: Viral=", length(receptors$Viral),
              ", NIU=", length(receptors$NIU))

  # Pull pseudobulk DE for each clone group; restrict to those receptors + DA substates.
  de_rows <- list()
  for (cg in clone_groups) {
    de_csv <- file.path(paths_t$results_tables,
                        paste0("tcr_pseudobulk_DE_", cg, ".csv"))
    if (!file.exists(de_csv)) {
      log_message("  Missing ", basename(de_csv), "; clone group ",
                  cg, " omitted.")
      next
    }
    d <- utils::read.csv(de_csv, stringsAsFactors = FALSE)
    if (!all(c("gene","log2FoldChange","padj","cluster") %in% colnames(d)))
      next
    d$target_id <- sub("^tcell_", "", d$cluster)
    d <- d[d$target_id %in% tcell_keep & d$gene %in% all_receptors, ,
           drop = FALSE]
    if (!nrow(d)) next
    d$clone_group <- cg
    de_rows[[cg]] <- d[, c("gene","log2FoldChange","padj","target_id","clone_group")]
  }
  if (length(de_rows) == 0L) {
    log_message("  viz_clone_group_lr_coupling: no clone-group DE rows ",
                "matched receptors; skipping.")
    return(invisible(NULL))
  }
  de <- do.call(rbind, de_rows)
  # Filter to receptors that move meaningfully in at least one cell of the
  # grid: padj < cutoff OR |log2FoldChange| >= 1. Without this every receptor
  # ever annotated by LIANA shows up as a row, including ~150 quiet ones that
  # add no signal and flatten the color scale.
  keep_genes <- de |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
      min_padj    = min(padj, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::filter(min_padj < padj_cut | max_abs_lfc >= 1) |>
    dplyr::pull(gene)
  if (!length(keep_genes)) {
    log_message("  viz_clone_group_lr_coupling: no receptors passed ",
                "padj<", padj_cut, " or |LFC|>=1; skipping plot.")
    return(invisible(NULL))
  }
  log_message("  Receptors with non-trivial effect: ", length(keep_genes),
              " (of ", length(all_receptors), ")")
  de <- de[de$gene %in% keep_genes, , drop = FALSE]
  # Annotate receptor direction (Viral / NIU / Both for receptors that appear
  # in both directions).
  de$lr_direction <- ifelse(
    de$gene %in% intersect(receptors$Viral, receptors$NIU), "Both",
    ifelse(de$gene %in% receptors$Viral, "Viral",
           ifelse(de$gene %in% receptors$NIU, "NIU", "Other")))
  de$lr_direction <- factor(de$lr_direction,
                            levels = c("Viral","NIU","Both","Other"))

  # Substate labels (e.g. "0: Naive/CM CD4") for x-axis readability.
  de$substate_display <- vapply(de$target_id,
    function(id) get_substate_display(cfg, "tcell", id), character(1))

  # Order genes by direction then by max |LFC|.
  ord_gene <- de |>
    dplyr::group_by(gene, lr_direction) |>
    dplyr::summarise(max_abs = max(abs(log2FoldChange), na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(lr_direction, dplyr::desc(max_abs))
  de$gene <- factor(de$gene, levels = rev(ord_gene$gene))
  de$sig  <- ifelse(!is.na(de$padj) & de$padj < padj_cut, "*", "")

  lim <- max(abs(de$log2FoldChange), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  hp <- ggplot2::ggplot(de,
                        ggplot2::aes(.data$substate_display, .data$gene,
                                     fill = .data$log2FoldChange)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig),
                       color = "black", size = 3.2, vjust = 0.7) +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white",
                                  high = "#B40426", midpoint = 0,
                                  limits = c(-lim, lim),
                                  name = "log2FC\n(in-group vs out)") +
    ggh4x::facet_grid2(rows = ggplot2::vars(.data$lr_direction),
                       cols = ggplot2::vars(.data$clone_group),
                       scales = "free_y", space = "free_y",
                       strip = ggh4x::strip_themed(
                         text_y = ggh4x::elem_list_text(face = "bold",
                                                        size = 9))) +
    ggplot2::labs(
      title    = paste0("Clone-group expression of myeloid->T cell LR ",
                        "receptors (G4b condition-unique pairs)"),
      subtitle = paste0("Rows = receptor genes split by direction of the LR ",
                        "pair (Viral-unique / NIU-unique / both). Cols = ",
                        "clone group. T cell targets restricted to DA ",
                        "clusters (", paste(tcell_keep, collapse = ", "),
                        "). * padj<", padj_cut, "."),
      x = "T cell target substate", y = "Receptor gene") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title       = ggplot2::element_text(face = "bold"),
                   plot.subtitle    = ggplot2::element_text(size = 9),
                   strip.text.x     = ggplot2::element_text(face = "bold",
                                                            size = 9),
                   strip.text.y     = ggplot2::element_text(face = "bold",
                                                            size = 9,
                                                            angle = 0,
                                                            hjust = 0),
                   axis.text.x      = ggplot2::element_text(angle = 35,
                                                            hjust = 1,
                                                            size = 8),
                   axis.text.y      = ggplot2::element_text(size = 8),
                   panel.spacing    = grid::unit(0.3, "lines"))
  n_genes <- length(unique(de$gene))
  save_pdf_png(hp, file.path(p$base,
                             "clone_group_lr_receptor_heatmap"),
               w = 11, h = max(7, n_genes * 0.25 + 3))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Module 7: Antigen-response signature focus on DA T cell substates
#   Companion to viz_tcr_genex_signatures(): same data
#   (tcr_signature_scores_by_clone_group.csv) but restricted to the DA
#   substates so the slot of biological interest reads cleanly. Dot size =
#   -log10(p), color = delta (median in - median out).
# ---------------------------------------------------------------------------
viz_clone_group_antigen_signatures <- function(cfg,
                                               tcell_keep = c("0","1","4","5")) {
  p <- .tcra_paths(cfg)
  sig_csv <- file.path(p$tables_t, "tcr_signature_scores_by_clone_group.csv")
  if (!file.exists(sig_csv)) return(invisible(NULL))
  sig <- utils::read.csv(sig_csv, stringsAsFactors = FALSE)
  if (!nrow(sig)) return(invisible(NULL))
  sig$target_id <- sub("^tcell_", "", sig$substate_key)
  sig <- sig[sig$target_id %in% tcell_keep, , drop = FALSE]
  if (!nrow(sig)) return(invisible(NULL))
  sig$substate_display <- vapply(sig$target_id,
    function(id) get_substate_display(cfg, "tcell", id), character(1))
  sig$neglog10_p <- -log10(pmax(sig$p, 1e-10))

  lim <- max(abs(sig$delta), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  dp <- ggplot2::ggplot(sig,
                        ggplot2::aes(.data$substate_display, .data$panel,
                                     size  = .data$neglog10_p,
                                     color = .data$delta)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(low = "#3B4CC0", mid = "white",
                                   high = "#B40426", midpoint = 0,
                                   limits = c(-lim, lim),
                                   name = "delta\n(in - out)") +
    ggplot2::scale_size_continuous(name = "-log10 p", range = c(1.2, 6)) +
    ggplot2::facet_wrap(~ clone_group, nrow = 1) +
    ggplot2::labs(
      title    = "Clone-group UCell signatures in DA T cell substates",
      subtitle = paste0("Antigen-response panels (effector/exhaustion/IFN-I/",
                        "TCR-activation/Treg) per clone group. T cell ",
                        "targets restricted to DA clusters (",
                        paste(tcell_keep, collapse = ", "), ")."),
      x = "T cell target substate", y = "Signature panel") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   strip.text    = ggplot2::element_text(face = "bold"),
                   axis.text.x   = ggplot2::element_text(angle = 35,
                                                         hjust = 1,
                                                         size = 8))
  save_pdf_png(dp, file.path(p$base,
                             "clone_group_antigen_signature_dotplot"),
               w = 11, h = 5)
  invisible(TRUE)
}

# ===========================================================================
# Figure 5 panels (Acts 1-3). Each function below is keyed to a single panel
# letter and is self-contained — call signature is uniformly viz_fig5_<X>(cfg).
# Outputs land in outputs/viz/eye/tcell/09_tcr_motif/tcr_advanced/fig5_<x>_*.{pdf,png}.
# ===========================================================================

# Panel A: VDJdb antigen-class overlay on intraocular T cell UMAP, split by
# phenotype. Uses scplotter::CellDimPlot with `highlight` so the antigen-
# annotated cells render on top of a greyed background — matches the
# established legacy `viz_vdjdb_overlay()` style and is rendered with a
# single CellDimPlot call (no manual ggplot reconstruction). The two
# phenotype facets are produced via `split_by = "Phenotype_2"`. Cell-level
# antigen class is stamped onto the object via .stamp_highlight_col so
# CellDimPlot can highlight by the resulting factor.
viz_fig5_a_vdjdb_umap_composite <- function(cfg) {
  p <- .tcra_paths(cfg)
  ann_csv <- file.path(p$tables_rep, "vdjdb_annotations.csv")
  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  if (!file.exists(ann_csv) || !file.exists(tcell_rds)) {
    log_message("  viz_fig5_a: required inputs missing; skipping.")
    return(invisible(NULL))
  }
  if (!requireNamespace("scplotter", quietly = TRUE)) {
    log_message("  viz_fig5_a: scplotter not installed; skipping.")
    return(invisible(NULL))
  }
  obj <- readRDS(tcell_rds)
  ann <- utils::read.csv(ann_csv, stringsAsFactors = FALSE)
  top_species <- cfg$tcr_advanced$fig5$panel_a_top_species %||%
                   c("CMV", "EBV", "VZV", "InfluenzaA", "SARS-CoV-2", "HTLV-1")

  # Reduction name (case-sensitive); accept either casing.
  rd_name <- intersect(c("UMAP", "umap"), names(obj@reductions))[1]
  if (is.na(rd_name)) {
    log_message("  viz_fig5_a: UMAP reduction missing; skipping.")
    return(invisible(NULL))
  }

  # Stamp a per-cell antigen-class factor onto the object. Cells with multiple
  # species hits get the highest-priority species so the highlight is
  # unambiguous in the legend.
  ann_priority <- ann[ann$antigen_species %in% top_species, , drop = FALSE]
  ann_priority$prio <- match(ann_priority$antigen_species, top_species)
  ann_priority <- ann_priority[order(ann_priority$barcode,
                                     ann_priority$prio), ]
  ann_priority <- ann_priority[!duplicated(ann_priority$barcode), ]
  bc <- colnames(obj)
  antigen <- ann_priority$antigen_species[match(bc, ann_priority$barcode)]
  antigen[is.na(antigen)] <- "Unannotated"
  obj$fig5_antigen <- factor(antigen, levels = c(top_species, "Unannotated"))
  # Restrict to NIU and Viral cells (drop Healthy/other for this panel).
  ph <- as.character(obj$Phenotype_2)
  keep_cells <- which(ph %in% c("NIU", "Viral"))
  if (length(keep_cells) == 0L) {
    log_message("  viz_fig5_a: no NIU/Viral cells; skipping.")
    return(invisible(NULL))
  }
  obj_sub <- obj[, keep_cells]
  obj_sub$Phenotype_2 <- factor(as.character(obj_sub$Phenotype_2),
                                levels = c("NIU", "Viral"))
  rm(obj); invisible(gc())

  pal <- vapply(top_species, function(s) {
    if (s %in% names(PATHOGEN_COLORS)) PATHOGEN_COLORS[[s]] else
      unname(PATHOGEN_COLORS["Other"])
  }, character(1))
  names(pal) <- top_species
  pal["Unannotated"] <- "grey88"

  # Highlight expression: any cell whose antigen != Unannotated.
  hl_expr <- 'fig5_antigen != "Unannotated"'

  # Use save_pdf_png (single labeled PDF + PNG) instead of dual_save_plot —
  # the latter writes a `_stripped.pdf` variant with black background that
  # is meant for slide overlays, not the manuscript panel. Other Fig 5
  # panels are all PDF+PNG, so this keeps the figure pack uniform.
  pA <- scplotter::CellDimPlot(obj_sub,
                               group_by    = "fig5_antigen",
                               split_by    = "Phenotype_2",
                               reduction   = rd_name,
                               highlight   = hl_expr,
                               palcolor    = pal,
                               bg_color    = "grey92",
                               pt_alpha    = 0.45,
                               highlight_size   = 1.6,
                               highlight_stroke = 0.3,
                               raster      = FALSE,
                               show_stat   = FALSE) +
    ggplot2::ggtitle(
      "Fig 5A: intraocular T cells by VDJdb antigen species",
      subtitle = paste0(
        "Highlighted: cells matching one of ", length(top_species),
        " curated pathogens (CellDimPlot highlight on VDJdb 2026-05-16). ",
        "Split: Phenotype_2.\n",
        "Caveat: VDJdb is biased toward well-studied specificities ",
        "(CMV n=27,545; InfluenzaA n=6,780; EBV n=6,488; SARS-CoV-2 ",
        "n=5,523 TRB records) and barely covers the clinical etiologies ",
        "(VZV n=31, a single IE62/ALWALPHAA epitope; HSV-2 n=17; HSV-1 ",
        "n=0). Absence of VZV/HSV annotation reflects this detection ",
        "limit, not biological absence.")) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8))
  save_pdf_png(pA, file.path(p$base, "fig5_a_vdjdb_umap_composite"),
               w = 13, h = 7)
  invisible(NULL)
}

# Panel C: HLA-B27 pathogenic motif recovery per subject.
# Per-subject fraction of intraocular T cells carrying TRAV21 + [YF]S[TS]
# CDR3-beta motif, grouped by HLA-B27 status x phenotype. Wilcoxon between
# HLA-B27+ NIU and other groups. Dots = subjects.
viz_fig5_c_hla_b27_per_subject <- function(cfg) {
  p <- .tcra_paths(cfg)
  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  if (!file.exists(tcell_rds)) {
    log_message("  viz_fig5_c: T cell object missing; skipping.")
    return(invisible(NULL))
  }
  obj <- readRDS(tcell_rds)
  meta <- obj@meta.data
  meta$barcode <- colnames(obj)
  meta$pathogenic <- flag_hla_b27_pathogenic_tcr(meta$CTgene, meta$CTaa)
  # Cohort group: HLA-B27+ NIU | HLA-B27- NIU | Viral | Healthy.
  # Etiology HLA_B27 is the carrier metadata flag for NIU.
  is_niu     <- meta$Phenotype_2 == "NIU"
  is_viral   <- meta$Phenotype_2 == "Viral"
  is_healthy <- meta$Phenotype_2 %in% c("Healthy", "HC")
  is_b27carrier <- (meta$Etiology %in% c("HLA_B27", "B27", "HLA-B27")) |
                    isTRUE(meta$HLA_B27_status %in%
                           c("Positive", "+", "TRUE", "POS"))
  meta$cohort_group <- NA_character_
  meta$cohort_group[is_niu &  is_b27carrier] <- "HLA-B27+ NIU"
  meta$cohort_group[is_niu & !is_b27carrier] <- "HLA-B27- NIU"
  meta$cohort_group[is_viral]                <- "Viral"
  meta$cohort_group[is_healthy]              <- "Healthy"
  rm(obj); invisible(gc())

  # Eye-only T cells (Figure 5 is intraocular).
  if ("Tissue_1" %in% colnames(meta)) {
    meta <- meta[meta$Tissue_1 == "Eye", , drop = FALSE]
  }
  meta <- meta[!is.na(meta$cohort_group), , drop = FALSE]
  if (!nrow(meta)) {
    log_message("  viz_fig5_c: no cells survived cohort filter; skipping.")
    return(invisible(NULL))
  }

  per_subj <- meta |>
    dplyr::group_by(.data$Subject, .data$cohort_group) |>
    dplyr::summarise(
      n_cells       = dplyr::n(),
      n_pathogenic  = sum(.data$pathogenic, na.rm = TRUE),
      frac          = .data$n_pathogenic / .data$n_cells,
      .groups = "drop") |>
    dplyr::filter(.data$n_cells >= 30L) # baseline T-cell denominator

  min_subj <- cfg$tcr_advanced$fig5$panel_c_min_subjects_per_group %||% 2L
  group_n <- per_subj |>
    dplyr::count(.data$cohort_group, name = "n_subj")
  keep_groups <- group_n$cohort_group[group_n$n_subj >= min_subj]
  per_subj <- per_subj[per_subj$cohort_group %in% keep_groups, , drop = FALSE]

  if (!nrow(per_subj)) {
    log_message("  viz_fig5_c: no group passed min_subjects floor; skipping.")
    return(invisible(NULL))
  }

  # Pairwise Wilcoxon: each group vs HLA-B27+ NIU (when present).
  comparisons <- list()
  if ("HLA-B27+ NIU" %in% keep_groups) {
    ref <- per_subj$frac[per_subj$cohort_group == "HLA-B27+ NIU"]
    for (g in setdiff(keep_groups, "HLA-B27+ NIU")) {
      alt <- per_subj$frac[per_subj$cohort_group == g]
      if (length(alt) < 2L) next
      pv <- tryCatch(stats::wilcox.test(ref, alt, exact = FALSE)$p.value,
                     error = function(e) NA_real_)
      comparisons[[g]] <- data.frame(group = g, p = pv)
    }
  }
  comp_df <- if (length(comparisons))
               dplyr::bind_rows(comparisons) else data.frame()
  if (nrow(comp_df) > 0L) {
    comp_df$label <- sprintf("p=%.3g", comp_df$p)
  }

  utils::write.csv(per_subj,
                   file.path(p$tables_rep,
                             "fig5_c_hla_b27_per_subject.csv"),
                   row.names = FALSE)

  group_order <- intersect(c("HLA-B27+ NIU", "HLA-B27- NIU", "Viral",
                              "Healthy"), keep_groups)
  per_subj$cohort_group <- factor(per_subj$cohort_group,
                                  levels = group_order)
  pal <- c("HLA-B27+ NIU" = unname(ETIOLOGY_GROUP_COLORS["NIU"]),
           "HLA-B27- NIU" = "#F8B6B8",
           "Viral"        = unname(ETIOLOGY_GROUP_COLORS["Viral"]),
           "Healthy"      = unname(ETIOLOGY_GROUP_COLORS["Healthy"]))
  y_top <- max(per_subj$frac, na.rm = TRUE)
  pC <- ggplot2::ggplot(per_subj,
                        ggplot2::aes(.data$cohort_group, .data$frac,
                                     color = .data$cohort_group)) +
    ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA,
                          color = "grey40") +
    ggplot2::geom_jitter(width = 0.18, height = 0, size = 2.4, alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                                expand = ggplot2::expansion(mult = c(0.05, 0.20))) +
    ggplot2::labs(
      title    = paste0("Fig 5C: intraocular T cells carrying HLA-B27 ",
                        "pathogenic TCR (TRAV21 + [YF]S[TS] CDR3-beta)"),
      subtitle = paste0("One dot = one subject. Subjects with <30 eye T cells ",
                        "excluded. Groups with < ", min_subj,
                        " subjects suppressed."),
      x = NULL,
      y = "Fraction of eye T cells with pathogenic TCR") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   axis.text.x   = ggplot2::element_text(angle = 15,
                                                         hjust = 1),
                   aspect.ratio  = 1)
  if (nrow(comp_df) > 0L) {
    comp_df$cohort_group <- factor(comp_df$group, levels = group_order)
    pC <- pC + ggplot2::geom_text(
      data = comp_df, ggplot2::aes(x = .data$cohort_group,
                                    y = y_top * 1.10,
                                    label = .data$label),
      inherit.aes = FALSE, size = 3.2)
  }
  # 1:1 aspect requested — make the saved canvas square so the panel area
  # matches the theme's aspect.ratio = 1.
  save_pdf_png(pC, file.path(p$base, "fig5_c_hla_b27_per_subject"),
               w = 6, h = 6)
  invisible(NULL)
}

# Panel F composite: tcrdist UMAP colored by top NIU-enriched GLIPH clusters,
# with VDJdb-annotated cells overlaid. Two layers on the same canvas:
#   layer 1 (color)  : top NIU-enriched GLIPH cluster IDs
#   layer 2 (shape)  : open circles for VDJdb-annotated cells (any species)
# Background = greyed tcrdist UMAP.
.fig5_load_candidates <- function(cfg) {
  csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                   "repertoire", "novel_tcr_candidates_ranked.csv")
  if (!file.exists(csv)) {
    log_message("  fig5 panel: novel_tcr_candidates_ranked.csv missing; ",
                "run cfg$steps$novel_tcr_discovery first.")
    return(NULL)
  }
  cand <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (!nrow(cand)) {
    log_message("  fig5 panel: candidate table empty; nothing to render.")
    return(NULL)
  }
  cand
}

# Derive a per-clone (CDR3b) ranking from the motif-level candidate table.
# Each CDR3 inherits its parent motif's composite_rank (best/lowest when
# the CDR3 belongs to multiple motif clusters). Per-clone n_cells and
# n_subjects are pulled from the Pgen table so we can break ties and
# annotate the row labels. Returns a data.frame of top-N CDR3 clones with:
#   CDR3b, parent_cluster_id, motif, composite_rank, composite_score,
#   is_b27_known, dominant_tcell_cluster, n_cells, n_subjects
.fig5_top_clones <- function(cfg, top_n) {
  cand <- .fig5_load_candidates(cfg)
  if (is.null(cand)) return(NULL)
  exclude_b27 <- isTRUE(cfg$tcr_advanced$fig5$novel_tcr_discovery$exclude_b27_in_viz)
  if (exclude_b27 && any(cand$is_b27_known, na.rm = TRUE))
    cand <- cand[!cand$is_b27_known, , drop = FALSE]
  if (!nrow(cand)) return(NULL)

  cls_csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                       "repertoire", "gliph_clusters.csv")
  if (!file.exists(cls_csv)) {
    log_message("  .fig5_top_clones: gliph_clusters.csv missing.")
    return(NULL)
  }
  cls <- utils::read.csv(cls_csv, stringsAsFactors = FALSE)
  cls <- cls[!is.na(cls$CDR3b) & nchar(cls$CDR3b) > 0L, , drop = FALSE]

  pgen_csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                        "repertoire", "olga_pgen_per_clone.csv")
  pgen <- if (file.exists(pgen_csv))
            utils::read.csv(pgen_csv, stringsAsFactors = FALSE)
          else NULL
  per_cdr3 <- if (!is.null(pgen) && nrow(pgen))
    pgen |>
      dplyr::group_by(.data$CDR3b) |>
      dplyr::summarise(
        n_cells    = dplyr::n(),
        n_subjects = dplyr::n_distinct(.data$Subject),
        median_log10_Pgen = stats::median(.data$log10_Pgen[
          is.finite(.data$log10_Pgen)]),
        .groups = "drop") else
    data.frame(CDR3b = character(0), n_cells = integer(0),
               n_subjects = integer(0), median_log10_Pgen = numeric(0),
               stringsAsFactors = FALSE)

  cand_keep <- cand[, c("cluster_id", "motif", "composite_rank",
                        "composite_score", "is_b27_known",
                        "dominant_tcell_cluster")]
  per_clone <- cls[, c("cluster_id", "CDR3b")] |>
    dplyr::distinct() |>
    dplyr::inner_join(cand_keep, by = "cluster_id") |>
    dplyr::group_by(.data$CDR3b) |>
    dplyr::slice_min(.data$composite_rank, n = 1L,
                     with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::rename(parent_cluster_id = "cluster_id") |>
    dplyr::left_join(per_cdr3, by = "CDR3b")
  per_clone$n_cells[is.na(per_clone$n_cells)] <- 0L
  per_clone$n_subjects[is.na(per_clone$n_subjects)] <- 0L

  # Rank clones by parent motif composite_rank ascending, then by clone
  # cell count descending (tiebreaker prefers more-expanded clones).
  per_clone <- per_clone[order(per_clone$composite_rank,
                               -per_clone$n_cells), , drop = FALSE]
  per_clone$clone_rank <- seq_len(nrow(per_clone))
  utils::head(per_clone, top_n)
}

# Panel F (novel discovery landscape) — was: tcrdist NIU-enriched cluster overlay.
# Plots every GLIPH motif candidate as a point in (NIU vs Viral FDR, antigen
# selection) space, colored by composite candidate rank. Labels every motif
# above a -log10(FDR) threshold (config knob) instead of the top-N by
# composite rank, so the points that visibly stand out get named. Dashed
# reference lines mark where the known HLA-B27 pathogenic motif sits so
# reviewers can see novel candidates relative to that benchmark.
viz_fig5_f_tcrdist_niu_vs_vdjdb <- function(cfg, top_n = 10L) {
  p <- .tcra_paths(cfg)
  cand <- .fig5_load_candidates(cfg)
  if (is.null(cand)) return(invisible(NULL))
  exclude_b27 <- isTRUE(cfg$tcr_advanced$fig5$novel_tcr_discovery$exclude_b27_in_viz)
  label_thr <- cfg$tcr_advanced$fig5$novel_tcr_discovery$panel_f_label_neglog_fdr_threshold %||% 5

  # Locate the B27 motif reference (used for dashed lines + as a benchmark
  # point even when excluded from the candidate ranking).
  b27_row <- cand[isTRUE(any(cand$is_b27_known)) & cand$is_b27_known, ,
                  drop = FALSE]
  if (nrow(b27_row)) b27_row <- b27_row[1, , drop = FALSE]

  d <- cand
  if (exclude_b27 && any(d$is_b27_known, na.rm = TRUE))
    d <- d[!d$is_b27_known, , drop = FALSE]
  if (!nrow(d)) {
    log_message("  viz_fig5_f: no novel candidates after B27 exclusion; skipping.")
    return(invisible(NULL))
  }
  # Fallback to the raw p column if a CSV from before the FDR change is
  # still on disk; otherwise prefer the BH-adjusted column.
  if (!"niu_vs_viral_fisher_FDR" %in% colnames(d)) {
    d$niu_vs_viral_fisher_FDR <- stats::p.adjust(d$niu_vs_viral_fisher_p,
                                                 method = "BH")
  }
  d$neglog_fdr <- -log10(pmax(d$niu_vs_viral_fisher_FDR, 1e-300))
  d$neg_pgen   <- -d$median_log10_Pgen
  d$shape      <- ifelse(d$is_b27_known, "HLA-B27 known", "Novel candidate")

  label_df <- d[!is.na(d$neglog_fdr) & d$neglog_fdr >= label_thr, ,
                drop = FALSE]
  log_message(sprintf("  viz_fig5_f: labeling %d motifs above -log10(FDR) = %g.",
                      nrow(label_df), label_thr))

  pF <- ggplot2::ggplot(d, ggplot2::aes(.data$neglog_fdr, .data$neg_pgen)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$composite_rank,
                                     size  = .data$n_cells,
                                     shape = .data$shape),
                        alpha = 0.85, stroke = 0.4) +
    ggplot2::scale_color_viridis_c(direction = -1,
                                   name = "Composite rank") +
    ggplot2::scale_shape_manual(values = c(`Novel candidate` = 16,
                                           `HLA-B27 known`   = 17),
                                name = NULL) +
    ggplot2::scale_size_continuous(range = c(2.5, 10),
                                   name = "Cells") +
    ggrepel::geom_text_repel(data = label_df,
                             ggplot2::aes(label = .data$motif),
                             size = 4, max.overlaps = Inf,
                             min.segment.length = 0.1,
                             box.padding = 0.55,
                             point.padding = 0.3,
                             force = 3,
                             show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = label_thr,
                        linetype = "dotted", color = "grey50",
                        linewidth = 0.3) +
    ggplot2::labs(
      title    = paste0("Fig 5F: novel autoimmune-pathogenic TCR motif ",
                        "discovery landscape"),
      subtitle = sprintf("GLIPH motif clusters. Up = lower Pgen (antigen-selected). Right = NIU vs Viral, BH-adjusted. Labeled: -log10(FDR) >= %g.",
                         label_thr),
      x = bquote("-log"[10]*"(FDR, NIU vs Viral)"),
      y = bquote("-median log"[10]*"(P"[gen]*")")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   panel.grid.minor = ggplot2::element_blank())

  if (nrow(b27_row)) {
    b27_fdr <- if ("niu_vs_viral_fisher_FDR" %in% colnames(b27_row))
                 b27_row$niu_vs_viral_fisher_FDR
               else stats::p.adjust(b27_row$niu_vs_viral_fisher_p,
                                    method = "BH")
    pF <- pF +
      ggplot2::geom_vline(xintercept = -log10(pmax(b27_fdr, 1e-300)),
                          linetype = "dashed", color = "#9D0208",
                          linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = -b27_row$median_log10_Pgen,
                          linetype = "dashed", color = "#9D0208",
                          linewidth = 0.4) +
      ggplot2::annotate("text", x = -log10(pmax(b27_fdr, 1e-300)),
                        y = max(d$neg_pgen, na.rm = TRUE),
                        label = "HLA-B27 motif benchmark",
                        color = "#9D0208", size = 3, hjust = -0.05,
                        vjust = 1.1)
  }
  save_pdf_png(pF, file.path(p$base, "fig5_f_novel_landscape"),
               w = 7, h = 7)
  invisible(NULL)
}

# Panel G (Pgen x clone size, colored by composite candidate rank) — was:
# NIU-enriched GLIPH motif logos. Reimagines the legacy
# .tcell_pgen_vs_clone_size_gliph() (R/75) with the candidate-ranking lens.
# Each point is a CDR3 from the per-clone Pgen table; inherits the parent
# GLIPH cluster's composite_rank from the ranked CSV. B27-known motif
# members get a red outline ring so they sit visibly alongside the novel
# discoveries.
viz_fig5_g_gliph_logos_with_b27 <- function(cfg) {
  p <- .tcra_paths(cfg)
  cand <- .fig5_load_candidates(cfg)
  if (is.null(cand)) return(invisible(NULL))
  exclude_b27 <- isTRUE(cfg$tcr_advanced$fig5$novel_tcr_discovery$exclude_b27_in_viz)
  top_n <- cfg$tcr_advanced$fig5$novel_tcr_discovery$top_n_for_viz %||% 10L

  pgen_csv <- file.path(get_target_paths(cfg, "all")$results_tables,
                        "repertoire", "olga_pgen_per_clone.csv")
  cls_csv  <- file.path(p$tables_rep, "gliph_clusters.csv")
  if (!file.exists(pgen_csv) || !file.exists(cls_csv)) {
    log_message("  viz_fig5_g: pgen or gliph_clusters table missing; skipping.")
    return(invisible(NULL))
  }
  pgen <- utils::read.csv(pgen_csv, stringsAsFactors = FALSE)
  cls  <- utils::read.csv(cls_csv,  stringsAsFactors = FALSE)
  pgen <- pgen[is.finite(pgen$log10_Pgen), , drop = FALSE]
  if (!nrow(pgen)) {
    log_message("  viz_fig5_g: empty Pgen table; skipping.")
    return(invisible(NULL))
  }

  # Per-CDR3 summary (one point per CDR3).
  per_cdr3 <- pgen |>
    dplyr::group_by(.data$CDR3b) |>
    dplyr::summarise(
      n_cells    = dplyr::n(),
      log10_Pgen = stats::median(.data$log10_Pgen),
      n_subjects = dplyr::n_distinct(.data$Subject),
      .groups = "drop"
    ) |>
    dplyr::mutate(sharing = ifelse(.data$n_subjects > 1, "public", "private"))

  # Map CDR3 -> best (lowest) composite_rank across its parent GLIPH clusters.
  # B27-known motifs are flagged separately so we can ring them.
  if (exclude_b27) {
    cand_active <- cand[!cand$is_b27_known, , drop = FALSE]
  } else {
    cand_active <- cand
  }
  cls_active <- cls[cls$cluster_id %in% cand_active$cluster_id, , drop = FALSE]
  cls_active <- dplyr::left_join(cls_active,
                                 dplyr::select(cand_active, "cluster_id",
                                               "composite_rank"),
                                 by = "cluster_id")
  cdr3_rank <- cls_active |>
    dplyr::group_by(.data$CDR3b) |>
    dplyr::summarise(composite_rank = min(.data$composite_rank, na.rm = TRUE),
                     .groups = "drop")
  per_cdr3 <- dplyr::left_join(per_cdr3, cdr3_rank, by = "CDR3b")

  b27_cluster_ids <- cand$cluster_id[isTRUE(any(cand$is_b27_known)) &
                                     cand$is_b27_known]
  if (length(b27_cluster_ids)) {
    b27_cdr3 <- unique(cls$CDR3b[cls$cluster_id %in% b27_cluster_ids])
  } else {
    b27_cdr3 <- character(0)
  }
  per_cdr3$is_b27_known <- per_cdr3$CDR3b %in% b27_cdr3

  # Background = CDR3s without any candidate parent (or in B27 motif when excluded).
  bg <- per_cdr3[is.na(per_cdr3$composite_rank), , drop = FALSE]
  fg <- per_cdr3[!is.na(per_cdr3$composite_rank), , drop = FALSE]

  # Top-N motif CDR3 labels — for each of the top-N motifs, pick its most
  # expanded member CDR3.
  top_clusters <- utils::head(cand_active[order(cand_active$composite_rank), ,
                                          drop = FALSE], top_n)
  cls_top <- cls[cls$cluster_id %in% top_clusters$cluster_id, , drop = FALSE]
  per_cdr3$top_motif <- per_cdr3$CDR3b %in% cls_top$CDR3b
  label_df <- per_cdr3 |>
    dplyr::filter(.data$top_motif, !is.na(.data$composite_rank)) |>
    dplyr::group_by(.data$composite_rank) |>
    dplyr::slice_max(.data$n_cells, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  size_vals <- c(public = 4.8, private = 1.4)

  pG <- ggplot2::ggplot() +
    ggplot2::geom_point(data = bg,
                        ggplot2::aes(.data$log10_Pgen, .data$n_cells,
                                     size = .data$sharing),
                        color = "grey80", alpha = 0.35) +
    ggplot2::geom_point(data = fg,
                        ggplot2::aes(.data$log10_Pgen, .data$n_cells,
                                     color = .data$composite_rank,
                                     size = .data$sharing),
                        alpha = 0.9) +
    ggplot2::geom_point(data = fg[fg$is_b27_known, , drop = FALSE],
                        ggplot2::aes(.data$log10_Pgen, .data$n_cells),
                        shape = 1, color = "#9D0208",
                        size = 5, stroke = 0.7) +
    ggrepel::geom_text_repel(data = label_df,
                             ggplot2::aes(.data$log10_Pgen, .data$n_cells,
                                          label = .data$CDR3b),
                             size = 2.4, max.overlaps = Inf,
                             min.segment.length = 0.1,
                             box.padding = 0.35,
                             force = 2,
                             show.legend = FALSE) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_viridis_c(direction = -1,
                                   name = "Composite rank") +
    ggplot2::scale_size_manual(values = size_vals, name = "Sharing") +
    ggplot2::labs(
      title    = "Fig 5G: generation probability vs clonal expansion, ranked candidates",
      subtitle = "Each point = CDR3-beta. Color = composite candidate rank of parent motif. Red ring = HLA-B27 known motif member.",
      x = bquote("log"[10]*"(P"[gen]*")"),
      y = "Cells per CDR3 (log10)") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   aspect.ratio  = 1)
  save_pdf_png(pG, file.path(p$base, "fig5_g_pgen_vs_clone_ranked"),
               w = 8, h = 7)
  invisible(NULL)
}

# Panel H (top-N clonal candidate occupancy vs NIU baseline). Two separate
# alluvial plots — by T cell cluster (H1) and by autoimmune indication
# (H2) — comparing the top-N CDR3 clones (ranked by parent GLIPH motif
# composite_rank, ties broken by clone size) against all NIU eye T cells.
viz_fig5_h_top_niu_motif_overlap <- function(cfg) {
  p <- .tcra_paths(cfg)
  top_n <- cfg$tcr_advanced$fig5$novel_tcr_discovery$panel_h_top_n %||% 100L
  clones <- .fig5_top_clones(cfg, top_n)
  if (is.null(clones) || !nrow(clones)) {
    log_message("  viz_fig5_h: no clonal candidates; skipping.")
    return(invisible(NULL))
  }
  cand_cdr3 <- unique(clones$CDR3b)
  log_message(sprintf("  viz_fig5_h: %d top clonal candidates (top-N motifs sweep).",
                      length(cand_cdr3)))

  # Load tcell metadata once. Restrict to NIU eye cells (baseline universe).
  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  if (!file.exists(tcell_rds)) {
    log_message("  viz_fig5_h: tcell object missing; skipping.")
    return(invisible(NULL))
  }
  meta <- local({
    obj <- readRDS(tcell_rds)
    md  <- obj@meta.data
    trb <- stringr::str_split(as.character(md$CTaa), "_", simplify = TRUE)
    md$CDR3b <- if (ncol(trb) >= 2) trb[, 2] else NA_character_
    rm(obj); invisible(gc())
    md
  })
  clu_col <- if ("knn.leiden.cluster" %in% colnames(meta))
               "knn.leiden.cluster" else NULL
  if (is.null(clu_col)) {
    log_message("  viz_fig5_h: knn.leiden.cluster missing; skipping.")
    return(invisible(NULL))
  }
  meta[[clu_col]] <- as.character(meta[[clu_col]])
  meta <- meta[meta$Phenotype_2 == "NIU", , drop = FALSE]
  if ("Tissue_1" %in% colnames(meta))
    meta <- meta[meta$Tissue_1 == "Eye", , drop = FALSE]
  if (!nrow(meta)) {
    log_message("  viz_fig5_h: no NIU eye cells available; skipping.")
    return(invisible(NULL))
  }
  is_cand <- !is.na(meta$CDR3b) & meta$CDR3b != "" &
             meta$CDR3b %in% cand_cdr3

  # ---- Helper: build a long frame with two groups (Candidates, Baseline)
  # and Fisher tests of (in-category × in-candidates) vs (not, not). The
  # baseline universe is all NIU eye cells; the candidate subset is the
  # is_cand mask.
  build_comp <- function(category_vec) {
    keep <- !is.na(category_vec) & category_vec != ""
    cat_v <- category_vec[keep]
    cand_v <- is_cand[keep]
    cats <- sort(unique(cat_v))
    rows <- lapply(cats, function(ct) {
      a <- sum(cand_v &  cat_v == ct)   # candidate, in category
      b <- sum(cand_v &  cat_v != ct)   # candidate, not
      c <- sum(!cand_v & cat_v == ct)   # baseline-only, in
      d2 <- sum(!cand_v & cat_v != ct)  # baseline-only, not
      ft <- tryCatch(stats::fisher.test(matrix(c(a, c, b, d2), nrow = 2)),
                     error = function(e) NULL)
      data.frame(
        category   = ct,
        n_cand     = a,
        n_baseline = c,
        OR         = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        p          = if (!is.null(ft)) ft$p.value           else NA_real_,
        stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    out$FDR <- stats::p.adjust(out$p, method = "BH")
    out
  }

  # Common cluster + etiology level orderings.
  clu_levels <- sort(unique(as.character(meta[[clu_col]])))
  eti_levels <- sort(unique(meta$Etiology[!is.na(meta$Etiology) &
                                           meta$Etiology != ""]))

  long_frame <- function(cat_vec, levels) {
    # Build a 2-bar (Candidates vs Baseline) stacked-bar data.frame.
    keep <- !is.na(cat_vec) & cat_vec != ""
    df_cand <- data.frame(group = "Candidates",
                          category = factor(cat_vec[keep & is_cand],
                                            levels = levels))
    df_base <- data.frame(group = "Baseline (all NIU eye T)",
                          category = factor(cat_vec[keep],
                                            levels = levels))
    rbind(
      .group_to_frac(df_cand, levels),
      .group_to_frac(df_base, levels)
    )
  }

  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    log_message("  viz_fig5_h: ggalluvial not installed; skipping.")
    return(invisible(NULL))
  }

  cat_pal_n <- function(n) {
    base <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
              "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f")
    if (n <= length(base)) base[seq_len(n)]
    else grDevices::colorRampPalette(base)(n)
  }

  # Build a lode-form data frame: one row per (category, axis), where axis 1
  # is "Overall" (baseline within-group fraction) and axis 2 is "Enriched
  # clones" (candidate within-group fraction). Each category becomes a
  # single alluvium that connects its proportion under each condition.
  lode_form <- function(cat_vec, levels, axis_labels) {
    keep <- !is.na(cat_vec) & cat_vec != ""
    cat_v  <- cat_vec[keep]
    cand_v <- is_cand[keep]
    base_frac <- vapply(levels, function(ct) mean(cat_v == ct),
                        numeric(1))
    cand_frac <- vapply(levels, function(ct) {
      denom <- sum(cand_v)
      if (!denom) return(0)
      sum(cand_v & cat_v == ct) / denom
    }, numeric(1))
    rbind(
      data.frame(alluvium = levels,
                 x        = axis_labels[1],
                 stratum  = levels,
                 y        = base_frac,
                 stringsAsFactors = FALSE),
      data.frame(alluvium = levels,
                 x        = axis_labels[2],
                 stratum  = levels,
                 y        = cand_frac,
                 stringsAsFactors = FALSE)
    )
  }

  axis_labels_h <- c("Overall (all NIU eye T)", "Enriched clones (top-N)")

  # ===================================================================
  # Plot 1: by T cell cluster (H1) — overall vs enriched alluvial
  # ===================================================================
  cluster_fdr  <- build_comp(meta[[clu_col]])
  utils::write.csv(cluster_fdr,
                   file.path(p$tables_rep,
                             "fig5_h_cluster_enrichment_vs_baseline.csv"),
                   row.names = FALSE)
  cluster_lode <- lode_form(meta[[clu_col]], clu_levels, axis_labels_h)
  cluster_lode$x       <- factor(cluster_lode$x, levels = axis_labels_h)
  cluster_lode$stratum <- factor(cluster_lode$stratum, levels = clu_levels)
  cluster_lode$alluvium <- factor(cluster_lode$alluvium, levels = clu_levels)

  pal_clu <- setNames(cat_pal_n(length(clu_levels)), clu_levels)
  sub_lbl <- if (exists("substate_labels"))
               substate_labels(cfg, "tcell", clu_levels)
             else paste0(clu_levels, ": cluster_", clu_levels)
  names(sub_lbl) <- clu_levels
  enr_sig <- cluster_fdr[!is.na(cluster_fdr$FDR) &
                          cluster_fdr$FDR < 0.05, , drop = FALSE]
  enr_sig$dir <- ifelse(enr_sig$OR > 1, "up", "down")
  sig_str <- if (nrow(enr_sig))
    paste(sprintf("%s %s (OR=%.2f, FDR=%.1g)",
                  enr_sig$category, enr_sig$dir, enr_sig$OR, enr_sig$FDR),
          collapse = "; ")
  else "no sig clusters"

  pH1 <- ggplot2::ggplot(cluster_lode,
                         ggplot2::aes(x = .data$x, y = .data$y,
                                      stratum  = .data$stratum,
                                      alluvium = .data$alluvium,
                                      fill     = .data$stratum)) +
    ggalluvial::geom_flow(alpha = 0.85, knot.pos = 0.4,
                          color = "white", linewidth = 0.2) +
    ggalluvial::geom_stratum(color = "grey25", linewidth = 0.4,
                             width = 0.32) +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(stratum)),
                       size = 3.4, color = "white", fontface = "bold") +
    ggplot2::scale_fill_manual(values = pal_clu, labels = sub_lbl,
                               name = "T cell cluster") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                name = "Within-column fraction",
                                expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::labs(
      title    = sprintf("Fig 5H1: cluster occupancy — overall vs enriched (top %d clones)",
                         length(cand_cdr3)),
      subtitle = paste0("Left = all NIU eye T cells. Right = cells in top-",
                        length(cand_cdr3),
                        " clonal candidates (CDR3 ranked by parent motif). Cand vs baseline FDR (BH): ",
                        sig_str, ".")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   axis.title.x  = ggplot2::element_blank(),
                   axis.text.x   = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  save_pdf_png(pH1, file.path(p$base, "fig5_h1_cluster_occupancy"),
               w = 8, h = 7)

  # ===================================================================
  # Plot 2: by autoimmune indication / etiology (H2) — overall vs enriched
  # ===================================================================
  eti_fdr <- build_comp(meta$Etiology)
  utils::write.csv(eti_fdr,
                   file.path(p$tables_rep,
                             "fig5_h_etiology_enrichment_vs_baseline.csv"),
                   row.names = FALSE)
  eti_lode <- lode_form(meta$Etiology, eti_levels, axis_labels_h)
  eti_lode$x        <- factor(eti_lode$x, levels = axis_labels_h)
  eti_lode$stratum  <- factor(eti_lode$stratum, levels = eti_levels)
  eti_lode$alluvium <- factor(eti_lode$alluvium, levels = eti_levels)

  pal_eti <- setNames(cat_pal_n(length(eti_levels)), eti_levels)
  sig_e <- eti_fdr[!is.na(eti_fdr$FDR) & eti_fdr$FDR < 0.05, , drop = FALSE]
  sig_e$dir <- ifelse(sig_e$OR > 1, "up", "down")
  sig_estr <- if (nrow(sig_e))
    paste(sprintf("%s %s (OR=%.2f, FDR=%.1g)",
                  sig_e$category, sig_e$dir, sig_e$OR, sig_e$FDR),
          collapse = "; ")
  else "no sig etiologies"

  pH2 <- ggplot2::ggplot(eti_lode,
                         ggplot2::aes(x = .data$x, y = .data$y,
                                      stratum  = .data$stratum,
                                      alluvium = .data$alluvium,
                                      fill     = .data$stratum)) +
    ggalluvial::geom_flow(alpha = 0.85, knot.pos = 0.4,
                          color = "white", linewidth = 0.2) +
    ggalluvial::geom_stratum(color = "grey25", linewidth = 0.4,
                             width = 0.32) +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(stratum)),
                       size = 3.4, color = "white", fontface = "bold") +
    ggplot2::scale_fill_manual(values = pal_eti, name = "Indication") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                name = "Within-column fraction",
                                expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::labs(
      title    = sprintf("Fig 5H2: indication occupancy — overall vs enriched (top %d clones)",
                         length(cand_cdr3)),
      subtitle = paste0("Left = all NIU eye T cells. Right = cells in top-",
                        length(cand_cdr3),
                        " clonal candidates (CDR3 ranked by parent motif). Cand vs baseline FDR (BH): ",
                        sig_estr, ".")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   axis.title.x  = ggplot2::element_blank(),
                   axis.text.x   = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  save_pdf_png(pH2, file.path(p$base, "fig5_h2_indication_occupancy"),
               w = 8, h = 7)
  invisible(NULL)
}

# Internal helper for Panel H: takes a frame with `group` + `category` factor
# columns and returns long counts + within-group fractions.
.group_to_frac <- function(df, levels) {
  if (!nrow(df)) {
    return(data.frame(group = character(0), category = factor(character(0),
                                                              levels = levels),
                      n = integer(0), frac = numeric(0)))
  }
  tab <- as.data.frame(table(group = df$group, category = df$category),
                       stringsAsFactors = FALSE)
  names(tab)[3] <- "n"
  tab$category <- factor(tab$category, levels = levels)
  tab <- tab[order(tab$group, tab$category), , drop = FALSE]
  tab$frac <- tab$n / pmax(1, ave(tab$n, tab$group, FUN = sum))
  tab
}

# Panel I (LIANA linkage — which myeloid populations signal into the
# top candidate T-cell clones). Uses every autoimmune-biased
# (NIU-driving) ligand-receptor pair from the cross-compartment LR
# analysis in R/47, regardless of which T-cell cluster LIANA flagged as
# the target. The displayed gene set is then narrowed to LR genes that
# are actually expressed (non-zero variance) across the candidate clones.
viz_fig5_i_clone_group_lr_with_b27 <- function(cfg) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE) ||
      !requireNamespace("viridisLite", quietly = TRUE)) {
    log_message("  viz_fig5_i: ComplexHeatmap/circlize/viridisLite not installed; skipping.")
    return(invisible(NULL))
  }
  # Render two variants: a focused top-25 + a broad top-100. The two N
  # values are configurable via cfg$...$panel_i_top_n_clones (focused)
  # and cfg$...$panel_i_top_n_clones_broad (broad).
  n_focus <- cfg$tcr_advanced$fig5$novel_tcr_discovery$panel_i_top_n_clones %||% 25L
  n_broad <- cfg$tcr_advanced$fig5$novel_tcr_discovery$panel_i_top_n_clones_broad %||% 100L
  variants <- list(
    list(n = n_focus, suffix = sprintf("top%d", n_focus)),
    list(n = n_broad, suffix = sprintf("top%d", n_broad))
  )
  for (v in variants) {
    .viz_fig5_i_one(cfg, v$n, v$suffix)
  }
  invisible(NULL)
}

# Worker for Panel I — renders one ComplexHeatmap variant for a given
# top-N clone count. Used twice from the wrapper (focused + broad views).
.viz_fig5_i_one <- function(cfg, top_n_clones, suffix) {
  p <- .tcra_paths(cfg)
  lia_csv <- file.path(cfg$paths$results_tables, "cross_compartment",
                       "liana_myeloid_to_tcell_combined.csv")
  if (!file.exists(lia_csv)) {
    log_message("  viz_fig5_i: LIANA combined table missing; skipping.")
    return(invisible(NULL))
  }
  lia <- utils::read.csv(lia_csv, stringsAsFactors = FALSE)

  pval_thresh <- cfg$tcr_advanced$fig5$novel_tcr_discovery$panel_i_niu_pvalue_threshold %||% 0.05
  bias_thresh <- cfg$tcr_advanced$fig5$novel_tcr_discovery$panel_i_disease_bias_threshold %||% 0

  # Keep every autoimmune-biased (NIU-driving) myeloid->T-cell LR pair,
  # regardless of which T-cell cluster LIANA assigned as the target. The
  # downstream zero-variance gene drop trims genes with no expression
  # signal across the displayed clones.
  lia_sig <- lia[grepl("^tcell_", lia$target) &
                 !is.na(lia$disease_bias_logfc) &
                 lia$disease_bias_logfc > bias_thresh &
                 !is.na(lia$cellphonedb.pvalue_NIU) &
                 lia$cellphonedb.pvalue_NIU < pval_thresh, , drop = FALSE]
  if (!nrow(lia_sig)) {
    log_message(sprintf("  viz_fig5_i (%s): no LIANA pairs pass NIU p<%g & bias>%g.",
                        suffix, pval_thresh, bias_thresh))
    return(invisible(NULL))
  }
  # Pull every gene from the NIU-biased LR pairs (ligand_complex AND
  # receptor_complex). The intersection with rownames(obj) below naturally
  # retains only T-cell-expressed genes, so dual-expressed ligands (e.g.,
  # HLA-DR) are kept on the T-cell side too.
  rec_genes <- unique(c(
    unlist(strsplit(lia_sig$receptor_complex, "_")),
    unlist(strsplit(lia_sig$ligand_complex,   "_"))
  ))
  rec_genes <- rec_genes[nzchar(rec_genes)]

  # Oversample candidate clones so that after dropping Unknown-disease
  # clones we still hit the requested top_n_clones row count.
  oversample <- cfg$tcr_advanced$fig5$novel_tcr_discovery$panel_i_oversample %||% 5L
  clones_all <- .fig5_top_clones(cfg, top_n_clones * oversample)
  if (is.null(clones_all) || !nrow(clones_all)) {
    log_message("  viz_fig5_i: no clonal candidates; skipping.")
    return(invisible(NULL))
  }

  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  if (!file.exists(tcell_rds)) {
    log_message("  viz_fig5_i: tcell object missing; skipping.")
    return(invisible(NULL))
  }
  # First pass: load metadata only, compute per-clone disease predominance
  # across the oversampled candidate pool. Used to filter Unknown clones
  # BEFORE picking the top-N rows so the displayed count matches request.
  predom_all <- local({
    obj <- readRDS(tcell_rds)
    md  <- obj@meta.data
    trb <- stringr::str_split(as.character(md$CTaa), "_", simplify = TRUE)
    md$CDR3b <- if (ncol(trb) >= 2) trb[, 2] else NA_character_
    rm(obj); invisible(gc())
    md_niu <- md[!is.na(md$CDR3b) & md$CDR3b != "" &
                 md$Phenotype_2 == "NIU" &
                 md$CDR3b %in% clones_all$CDR3b &
                 !is.na(md$Etiology) & md$Etiology != "", , drop = FALSE]
    if (!nrow(md_niu)) return(NULL)
    md_niu |>
      dplyr::count(.data$CDR3b, .data$Etiology, name = "n") |>
      dplyr::group_by(.data$CDR3b) |>
      dplyr::mutate(dom_frac = .data$n / sum(.data$n)) |>
      dplyr::slice_max(.data$n, n = 1L, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::rename(disease_predom  = "Etiology",
                    predom_dom_frac = "dom_frac") |>
      dplyr::select("CDR3b", "disease_predom", "predom_dom_frac")
  })
  clones_all <- dplyr::left_join(clones_all, predom_all, by = "CDR3b")
  before_n <- nrow(clones_all)
  clones_all <- clones_all[!is.na(clones_all$disease_predom), , drop = FALSE]
  log_message(sprintf(
    "  viz_fig5_i (%s): %d / %d candidates retained after Unknown-disease filter.",
    suffix, nrow(clones_all), before_n))
  if (!nrow(clones_all)) {
    log_message("  viz_fig5_i: no clones after Unknown filter; skipping.")
    return(invisible(NULL))
  }
  clones <- utils::head(clones_all, top_n_clones)
  if (nrow(clones) < top_n_clones)
    log_message(sprintf(
      "  viz_fig5_i (%s): only %d of requested %d clones available post-filter; consider raising panel_i_oversample.",
      suffix, nrow(clones), top_n_clones))
  clones$is_public <- clones$n_subjects > 1L

  # Second pass: now load expression only for the final clone set.
  loaded <- local({
    obj <- readRDS(tcell_rds)
    md  <- obj@meta.data
    trb <- stringr::str_split(as.character(md$CTaa), "_", simplify = TRUE)
    md$CDR3b <- if (ncol(trb) >= 2) trb[, 2] else NA_character_
    md$barcode <- colnames(obj)
    in_clone <- !is.na(md$CDR3b) & md$CDR3b %in% clones$CDR3b &
                md$Phenotype_2 == "NIU"
    bcs_have <- md$barcode[in_clone]
    have_genes <- intersect(rec_genes, rownames(obj))
    missing_g  <- setdiff(rec_genes, have_genes)
    if (length(missing_g))
      log_message(sprintf("  viz_fig5_i: %d LR genes missing from T cell object (kept %d).",
                          length(missing_g), length(have_genes)))
    if (!length(have_genes) || !length(bcs_have)) {
      rm(obj); invisible(gc())
      return(NULL)
    }
    em <- Seurat::FetchData(obj[, bcs_have], vars = have_genes,
                            layer = "data")
    em$barcode <- rownames(em)
    md_sub <- md[in_clone, c("barcode", "CDR3b", "Etiology"), drop = FALSE]
    rm(obj); invisible(gc())
    list(em = em, md = md_sub)
  })
  if (is.null(loaded)) {
    log_message("  viz_fig5_i: no expression data; skipping.")
    return(invisible(NULL))
  }
  em <- loaded$em
  md_sub <- loaded$md

  # Long expression frame -> per-clone mean per receptor.
  expr_long <- tidyr::pivot_longer(em, cols = -.data$barcode,
                                    names_to = "receptor",
                                    values_to = "expr")
  expr_long$CDR3b <- md_sub$CDR3b[match(expr_long$barcode, md_sub$barcode)]
  clone_expr <- expr_long |>
    dplyr::group_by(.data$CDR3b, .data$receptor) |>
    dplyr::summarise(mean_expr = mean(.data$expr, na.rm = TRUE),
                     .groups = "drop")
  clone_score <- clone_expr |>
    dplyr::group_by(.data$CDR3b) |>
    dplyr::summarise(liana_engagement = mean(.data$mean_expr, na.rm = TRUE),
                     .groups = "drop")
  clones <- dplyr::left_join(clones, clone_score, by = "CDR3b")

  utils::write.csv(
    clones[, c("CDR3b", "parent_cluster_id", "motif",
               "composite_rank", "clone_rank",
               "liana_engagement", "n_cells", "n_subjects", "is_public",
               "dominant_tcell_cluster", "disease_predom",
               "predom_dom_frac")],
    file.path(p$tables_rep,
              sprintf("fig5_i_clone_liana_engagement_%s.csv", suffix)),
    row.names = FALSE)

  # Build wide matrix and column-z-score so per-receptor patterns are
  # visible regardless of absolute expression magnitude.
  mat_wide <- tidyr::pivot_wider(clone_expr,
                                 names_from = "receptor",
                                 values_from = "mean_expr",
                                 values_fill = 0)
  mat <- as.matrix(mat_wide[, -1, drop = FALSE])
  rownames(mat) <- mat_wide$CDR3b
  missing_rows <- setdiff(clones$CDR3b, rownames(mat))
  if (length(missing_rows)) {
    add <- matrix(0, nrow = length(missing_rows), ncol = ncol(mat),
                  dimnames = list(missing_rows, colnames(mat)))
    mat <- rbind(mat, add)
  }
  mat <- mat[clones$CDR3b, , drop = FALSE]
  # Drop zero-variance genes — these contribute no signal and just dilute
  # the column dendrogram + axis labels.
  col_var <- apply(mat, 2, stats::var)
  keep_cols <- !is.na(col_var) & col_var > 0
  n_dropped <- sum(!keep_cols)
  if (n_dropped) {
    log_message(sprintf("  viz_fig5_i (%s): dropping %d zero-variance genes (kept %d).",
                        suffix, n_dropped, sum(keep_cols)))
    mat <- mat[, keep_cols, drop = FALSE]
  }
  if (!ncol(mat)) {
    log_message("  viz_fig5_i: no genes left after zero-variance filter; skipping.")
    return(invisible(NULL))
  }
  mat_z <- scale(mat)
  mat_z[!is.finite(mat_z)] <- 0
  # Cap z extremes for readable color scaling.
  z_clip <- 3
  mat_z[mat_z >  z_clip] <-  z_clip
  mat_z[mat_z < -z_clip] <- -z_clip

  # ---- Annotation palettes -------------------------------------------
  cat_pal_n <- function(n) {
    base <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
              "#8c564b","#e377c2","#17becf","#bcbd22","#7f7f7f")
    if (n <= length(base)) base[seq_len(n)]
    else grDevices::colorRampPalette(base)(n)
  }
  # Dominant cluster uses RColorBrewer "Paired" so the categorical chips
  # carry the same palette family used elsewhere in the manuscript.
  paired_pal_n <- function(n) {
    raw <- RColorBrewer::brewer.pal(max(3L, min(12L, n)), "Paired")
    if (n <= length(raw)) raw[seq_len(n)]
    else grDevices::colorRampPalette(raw)(n)
  }
  clu_levels <- sort(unique(as.character(clones$dominant_tcell_cluster)))
  clu_pal <- setNames(paired_pal_n(length(clu_levels)), clu_levels)
  eti_levels <- sort(unique(as.character(clones$disease_predom)))
  eti_pal <- setNames(cat_pal_n(length(eti_levels)), eti_levels)
  pub_pal <- c(`Public` = "#222222", `Private` = "grey80")

  # Viridis maps z-scores from -z_clip (dark purple) to +z_clip (yellow).
  expr_col <- circlize::colorRamp2(
    seq(-z_clip, z_clip, length.out = 11),
    viridisLite::viridis(11))
  # Clone size: inferno sequential ramp at 0%, 50%, 100% quantiles.
  size_col <- circlize::colorRamp2(
    stats::quantile(clones$n_cells, c(0, 0.5, 1), na.rm = TRUE),
    viridisLite::inferno(3))
  # Overall rank: low (best) = black, high (worst) = white.
  rank_col <- circlize::colorRamp2(
    range(clones$clone_rank, na.rm = TRUE),
    c("black", "white"))

  # Row labels on the LEFT; row annotation block on the LEFT too so the
  # chips sit next to the CDR3 labels.
  la <- ComplexHeatmap::rowAnnotation(
    `Disease predom.`    = clones$disease_predom,
    `Public clone`       = ifelse(clones$is_public, "Public", "Private"),
    `Dominant cluster`   = clones$dominant_tcell_cluster,
    `Clone size (cells)` = clones$n_cells,
    `Overall rank`       = clones$clone_rank,
    col = list(
      `Disease predom.`    = eti_pal,
      `Public clone`       = pub_pal,
      `Dominant cluster`   = clu_pal,
      `Clone size (cells)` = size_col,
      `Overall rank`       = rank_col),
    annotation_name_gp   = grid::gpar(fontsize = 8, fontface = "bold"),
    annotation_name_side = "top",
    annotation_width = grid::unit(c(0.45, 0.45, 0.45, 0.45, 0.45), "cm"),
    gap = grid::unit(1.2, "mm")
  )

  n_col_split <- cfg$tcr_advanced$fig5$novel_tcr_discovery$panel_i_column_split %||% 4L
  n_col_split <- max(1L, min(as.integer(n_col_split), ncol(mat_z)))

  # Compute the column clustering ourselves so we can map each gene to a
  # module ID, then label each module by gene-level functional families.
  col_dist <- stats::dist(t(mat_z))
  col_hc   <- stats::hclust(col_dist, method = "ward.D2")
  col_groups <- stats::cutree(col_hc, k = n_col_split)  # named: gene -> group_id

  # Classify each gene symbol into a functional family. Uses HGNC-style
  # prefixes plus curated allowlists. Allowlists are checked first so
  # specific functional markers (PDCD1, CD28, GZMB) override generic
  # prefix matches (e.g., a gene starting with "CD" defaulting to "Other").
  classify_gene <- function(g) {
    checkpoint <- c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","BTLA","VSIR",
                    "CD274","PDCD1LG2","CD276","VTCN1","SIRPG","SIRPA")
    costim     <- c("CD28","CD80","CD86","CD40","CD40LG","ICOS","ICOSLG",
                    "CD27","CD70","TNFRSF4","TNFSF4","TNFRSF9","TNFSF9",
                    "TNFRSF14","TNFSF14","CD226","CD244","CD2","CD58",
                    "SLAMF1","SLAMF6","CD96","CRTAM","CD48")
    cytotox    <- c("GZMA","GZMB","GZMH","GZMK","PRF1","GNLY","NKG7",
                    "KLRD1","KLRK1","KLRG1","FASLG","FAS","TNFSF10",
                    "NCR3","CD160","KIR2DL1","KIR2DL3","KIR3DL1")
    tcr_core   <- c("CD3D","CD3E","CD3G","CD247","CD4","CD8A","CD8B","CD7",
                    "CD5","CD6","CD45","PTPRC")
    out <- character(length(g)); out[] <- "Other"
    out[g %in% checkpoint] <- "Checkpoint"
    out[g %in% costim]     <- "Costim"
    out[g %in% cytotox]    <- "Cytotox"
    out[g %in% tcr_core]   <- "TCR/CD-core"
    is_o <- out == "Other"
    out[is_o & grepl("^HLA-D", g)]                <- "HLA-II"
    out[is_o & grepl("^HLA-[ABCEFG]$", g)]        <- "HLA-I"
    is_o <- out == "Other"
    # Cytokine receptors before cytokines so IL2RA / IL10RA / TNFRSF1B
    # land as receptors rather than as the parent cytokine family.
    out[is_o & grepl("^IL[0-9]+R[A-Z]?[0-9]?$|^IFNAR[12]?$|^IFNGR[12]?$|^TGFBR[123]?$|^TNFRSF[0-9]+[A-Z]?$",
                     g)] <- "Cytokine-R"
    is_o <- out == "Other"
    out[is_o & grepl("^IL[0-9]+$|^IFN[ABGL]?[0-9]?$|^TNF$|^TNFSF[0-9]+[A-Z]?$|^TGFB[123]?$",
                     g)] <- "Cytokine"
    is_o <- out == "Other"
    out[is_o & grepl("^CCL[0-9]|^CXCL[0-9]|^CCR[0-9]|^CXCR[0-9]|^XCL|^XCR|^CX3C",
                     g)] <- "Chemokine"
    is_o <- out == "Other"
    out[is_o & grepl("^ITG[AB]|^SELL$|^SELPLG$|^ICAM[0-9]|^VCAM[0-9]|^PECAM[0-9]|^CADM[0-9]|^F11R$|^ALCAM$|^JAM[123]?$",
                     g)] <- "Adhesion"
    is_o <- out == "Other"
    out[is_o & grepl("^S100A[0-9]|^LGALS[0-9]|^LYZ$|^MMP[0-9]|^MIF$|^EBI3$|^SPP1$",
                     g)] <- "Inflammatory"
    out
  }

  # Per module: count gene-level families. Pick the top family; if a
  # secondary family contributes >=40% of the top, show both. Drop
  # "Other" from the dominance vote so it never wins a module name, but
  # if a module is *entirely* Other, fall back to that label.
  module_family <- vapply(seq_len(n_col_split), function(g) {
    genes_in_mod <- names(col_groups)[col_groups == g]
    fams <- classify_gene(genes_in_mod)
    fams_named <- fams[fams != "Other"]
    if (!length(fams_named)) return("Other")
    tbl <- sort(table(fams_named), decreasing = TRUE)
    top1 <- names(tbl)[1]
    if (length(tbl) > 1L && tbl[2] >= tbl[1] * 0.40) {
      sprintf("%s + %s", top1, names(tbl)[2])
    } else top1
  }, character(1))
  log_message(sprintf("  viz_fig5_i (%s): module-family map = %s",
                      suffix,
                      paste(sprintf("%s:%s", LETTERS[seq_len(n_col_split)],
                                    module_family), collapse = "; ")))

  # Order modules by first appearance in dendrogram leaf order (so visual
  # left-to-right matches the dendrogram traversal).
  leaf_genes <- colnames(mat_z)[col_hc$order]
  group_by_leaf <- col_groups[leaf_genes]
  first_pos <- vapply(seq_len(n_col_split),
                      function(g) match(TRUE, group_by_leaf == g),
                      integer(1))
  group_visual_order <- order(first_pos)

  # Build pretty labels in visual order: "FAMILY (Module L)".
  visual_letters <- LETTERS[seq_len(n_col_split)]
  module_label <- setNames(character(n_col_split), seq_len(n_col_split))
  for (i in seq_along(group_visual_order)) {
    g <- group_visual_order[i]
    module_label[as.character(g)] <- sprintf(
      "%s\n(Module %s)", module_family[g], visual_letters[i])
  }
  split_factor <- factor(module_label[as.character(col_groups)],
                         levels = module_label[as.character(group_visual_order)])
  split_labels <- if (n_col_split > 1L) levels(split_factor) else NULL

  hm <- ComplexHeatmap::Heatmap(
    mat_z,
    name  = "Column z-score\n(mean log-norm expr)",
    col   = expr_col,
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows    = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows    = "ward.D2",
    clustering_method_columns = "ward.D2",
    show_row_dend    = TRUE,
    show_column_dend = TRUE,
    row_dend_side    = "right",
    column_dend_side = "top",
    row_dend_width   = grid::unit(1.4, "cm"),
    column_dend_height = grid::unit(1.4, "cm"),
    row_names_side   = "left",
    show_row_names   = TRUE,
    row_names_gp     = grid::gpar(fontsize = if (top_n_clones > 50) 6 else 8),
    column_names_gp  = grid::gpar(fontsize = if (ncol(mat_z) > 150) 4 else
                                            if (ncol(mat_z) > 80)  5 else 7),
    column_names_rot = 45,
    column_names_side = "bottom",
    show_column_names = TRUE,
    left_annotation  = la,
    border = TRUE,
    column_split = if (n_col_split > 1L) split_factor else NULL,
    column_gap   = grid::unit(2.5, "mm"),
    column_title = split_labels,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    row_title = NULL,
    heatmap_legend_param = list(
      legend_height = grid::unit(3, "cm"),
      title_gp = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8))
  )

  w_in <- max(11, ncol(mat_z) * 0.13 + 7)
  h_in <- max(6,  nrow(mat_z) * (if (top_n_clones > 50) 0.18 else 0.32) + 3)

  overall_title <- sprintf(
    "Fig 5I (%s clones): T-cell expression of all NIU-enriched LR genes (any target cluster), column z-scored",
    suffix)

  base_path <- file.path(p$base,
                         sprintf("fig5_i_clone_liana_engagement_%s",
                                 suffix))
  grDevices::pdf(paste0(base_path, ".pdf"), width = w_in, height = h_in)
  ComplexHeatmap::draw(hm, merge_legend = TRUE,
                       column_title = overall_title,
                       column_title_gp = grid::gpar(fontsize = 12,
                                                    fontface = "bold"))
  grDevices::dev.off()
  grDevices::png(paste0(base_path, ".png"), width = w_in, height = h_in,
                 units = "in", res = 300)
  ComplexHeatmap::draw(hm, merge_legend = TRUE,
                       column_title = overall_title,
                       column_title_gp = grid::gpar(fontsize = 12,
                                                    fontface = "bold"))
  grDevices::dev.off()
  log_message("  Saved: ", basename(base_path), ".{pdf,png} (",
              nrow(mat_z), " clones x ", ncol(mat_z), " receptors).")
  invisible(NULL)
}

# Panel D: eye-to-blood alluvial colored by antigen class (extends
# .tcell_alluvial() from R/75 by adding the VDJdb annotation + HLA-B27 flag
# as the antigen-class color rather than etiology).
viz_fig5_d_alluvial_by_antigen_class <- function(cfg) {
  p <- .tcra_paths(cfg)
  ct_path <- "outputs/tables/repertoire/TCR_top_expanded_eye_celltype.csv"
  ann_csv <- file.path(p$tables_rep, "vdjdb_annotations.csv")
  if (!file.exists(ct_path)) {
    log_message("  viz_fig5_d: ", ct_path, " missing; skipping.")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    log_message("  viz_fig5_d: ggalluvial not installed; skipping.")
    return(invisible(NULL))
  }
  ct <- utils::read.csv(ct_path, stringsAsFactors = FALSE)
  # Two-stage filter: Eye rows must be tcell-substate (so the alluvial's
  # left axis is the T cell substate); Blood rows can be any substate label
  # (the table merges across the broad atlas, so blood rows often use
  # full-atlas labels rather than tcell_X). We only need blood for the
  # "shared with blood?" axis, so substate identity there is not critical.
  ct <- ct[ct$tissue %in% c("Eye", "Blood"), , drop = FALSE]
  eye_keep <- ct$tissue == "Eye" & grepl("^tcell_", ct$substate)
  blood_keep <- ct$tissue == "Blood"
  ct <- ct[eye_keep | blood_keep, , drop = FALSE]
  if (!nrow(ct)) return(invisible(NULL))

  # Pull HLA-B27 pathogenic clone IDs + per-clone VDJdb species (max-priority
  # per clone, same order as Panel A).
  top_species <- cfg$tcr_advanced$fig5$panel_a_top_species %||%
                  c("CMV", "EBV", "VZV", "InfluenzaA", "SARS-CoV-2", "HTLV-1")
  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  clone_antigen <- character(0)
  b27_clone_ids <- character(0)
  if (file.exists(tcell_rds)) {
    obj <- readRDS(tcell_rds)
    meta <- obj@meta.data
    meta$barcode <- colnames(obj)
    # B27 pathogenic clones
    flag <- flag_hla_b27_pathogenic_tcr(meta$CTgene, meta$CTaa)
    b27_clone_ids <- unique(meta$CTstrict[flag])
    # VDJdb per-clone
    if (file.exists(ann_csv)) {
      ann <- utils::read.csv(ann_csv, stringsAsFactors = FALSE)
      ann <- ann[ann$antigen_species %in% top_species, , drop = FALSE]
      ann$prio <- match(ann$antigen_species, top_species)
      ann <- ann[order(ann$barcode, ann$prio), ]
      ann <- ann[!duplicated(ann$barcode), ]
      bc_to_clone <- setNames(meta$CTstrict, meta$barcode)
      ann$CTstrict <- bc_to_clone[ann$barcode]
      ann <- ann[!is.na(ann$CTstrict), ]
      ann <- ann[order(ann$CTstrict, ann$prio), ]
      ann <- ann[!duplicated(ann$CTstrict), ]
      clone_antigen <- setNames(ann$antigen_species, ann$CTstrict)
    }
    rm(obj); invisible(gc())
  }

  ct$antigen_class <- "Unannotated"
  has_v <- ct$clone_id %in% names(clone_antigen)
  ct$antigen_class[has_v] <- unname(clone_antigen[ct$clone_id[has_v]])
  ct$antigen_class[ct$clone_id %in% b27_clone_ids] <- "HLA-B27"

  agg <- ct |>
    dplyr::group_by(.data$clone_id, .data$antigen_class, .data$tissue) |>
    dplyr::slice_max(.data$n_cells, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(.data$clone_id, .data$antigen_class, .data$tissue,
                  .data$substate) |>
    tidyr::pivot_wider(names_from = .data$tissue,
                       values_from = .data$substate)
  if (!"Eye" %in% colnames(agg)) return(invisible(NULL))
  # Blood column may be absent entirely when none of the top expanded eye
  # clones have a paired-blood row (small cohort case). Synthesize it as NA
  # so the alluvial still renders with everything in "eye-only".
  if (!"Blood" %in% colnames(agg)) agg$Blood <- NA_character_
  agg$blood_status <- ifelse(is.na(agg$Blood), "eye-only", "shared with blood")
  agg <- agg[!is.na(agg$Eye), , drop = FALSE]
  flow <- agg |>
    dplyr::count(.data$Eye, .data$antigen_class, .data$blood_status,
                 name = "n_clones")
  if (!nrow(flow)) return(invisible(NULL))

  # Panel D palette matches Panel A: same species hex via PATHOGEN_COLORS,
  # same grey88 for Unannotated, plus the HLA-B27 subtype color so the two
  # panels read as one antigen-color system.
  pal <- vapply(top_species, function(s) {
    if (s %in% names(PATHOGEN_COLORS)) PATHOGEN_COLORS[[s]] else
      unname(PATHOGEN_COLORS["Other"])
  }, character(1))
  names(pal) <- top_species
  pal["HLA-B27"]    <- unname(ETIOLOGY_SUBTYPE_COLORS["HLA_B27"])
  pal["Unannotated"]<- "grey88"

  pD <- ggplot2::ggplot(flow,
                        ggplot2::aes(axis1 = .data$Eye,
                                     axis2 = .data$antigen_class,
                                     axis3 = .data$blood_status,
                                     y = .data$n_clones)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data$antigen_class),
                              width = 1/6, alpha = 0.85) +
    ggalluvial::geom_stratum(width = 1/6, alpha = 0.55,
                             fill = "grey95", color = "grey25") +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(stratum)),
                       size = 2.4) +
    ggplot2::scale_x_discrete(limits = c("Eye substate",
                                          "Antigen class",
                                          "Blood detection"),
                              expand = c(.04, .04)) +
    ggplot2::scale_fill_manual(values = pal, na.value = "grey75",
                               name = "Antigen class") +
    ggplot2::labs(
      title    = "Fig 5D: top expanded eye T-cell clones, antigen class -> blood detection",
      subtitle = paste0("Clone-level traces (per CTstrict). Antigen class ",
                        "from VDJdb (>=10 hits/species) + HLA-B27 pathogenic ",
                        "motif (TRAV21 + [YF]S[TS]). 'Unannotated' = neither."),
      y = "Clones") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   axis.title.x  = ggplot2::element_blank())
  save_pdf_png(pD, file.path(p$base,
                             "fig5_d_alluvial_by_antigen_class"),
               w = 11, h = 7.5)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Supplemental: blood<->eye sharing of expanded TCR clones (alluvial).
# Strata are the two tissues (Blood, Eye); each ribbon is one expanded eye TCR
# clone, coloured by whether it is also detected in the paired blood (shared)
# or eye-restricted. Each disease group is normalised to 100% so NIU and Viral
# are comparable on one axis, with a per-subject Wilcoxon test on the shared
# fraction. Reads outputs/tables/repertoire/TCR_top_expanded_eye.csv.
# ---------------------------------------------------------------------------

# Per-subject Phenotype_2 lookup from the sample metadata (one row/subject).
.fig5_supp_subject_pheno <- function(cfg) {
  md_path <- file.path(cfg$paths$inputs %||% "inputs", "data", "metadata.csv")
  if (!file.exists(md_path)) md_path <- "inputs/data/metadata.csv"
  if (!file.exists(md_path)) return(NULL)
  md <- utils::read.csv(md_path, stringsAsFactors = FALSE, check.names = TRUE)
  if (!all(c("Subject", "Phenotype_2") %in% colnames(md))) return(NULL)
  dplyr::distinct(md[, c("Subject", "Phenotype_2")])
}

viz_fig5_s3b_eye_blood_clone_alluvial <- function(cfg) {
  p <- .tcra_paths(cfg)
  exp_csv <- file.path(p$tables_rep, "TCR_top_expanded_eye.csv")
  groups_keep <- cfg$tcr_advanced$fig5_supp$groups_keep %||% c("NIU", "Viral")
  min_cells <- cfg$tcr_advanced$fig5_supp$expanded_clone_min_cells %||% 3L
  if (!file.exists(exp_csv)) {
    log_message("  viz_fig5_s3b: TCR_top_expanded_eye.csv missing; skipping.")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    log_message("  viz_fig5_s3b: ggalluvial not installed; skipping.")
    return(invisible(NULL))
  }
  ex <- utils::read.csv(exp_csv, stringsAsFactors = FALSE)
  ex <- ex[!is.na(ex$n_cells_eye) & ex$n_cells_eye >= min_cells, , drop = FALSE]
  pheno <- .fig5_supp_subject_pheno(cfg)
  if (!is.null(pheno))
    ex <- dplyr::left_join(ex, pheno, by = c("subject" = "Subject"))
  else ex$Phenotype_2 <- NA_character_
  ex <- ex[ex$Phenotype_2 %in% groups_keep, , drop = FALSE]
  if (nrow(ex) == 0) {
    log_message("  viz_fig5_s3b: no expanded clones map to NIU/Viral; skipping.")
    return(invisible(NULL))
  }
  ex$found    <- as.logical(as.character(ex$found_in_blood))
  ex$found[is.na(ex$found)] <- FALSE
  ex$clone_uid <- paste(ex$subject, ex$clone_id, sep = "|")
  ex$sharing   <- ifelse(ex$found, "Shared (eye + blood)", "Eye-restricted")

  # Per-group totals + shared fraction, and a NIU-vs-Viral test on the
  # per-subject shared fraction (same statistic as the S3 boxplot).
  ann <- ex |>
    dplyr::group_by(.data$Phenotype_2) |>
    dplyr::summarise(n = dplyr::n(), n_shared = sum(.data$found),
                     pct = round(100 * mean(.data$found)), .groups = "drop")
  subj <- ex |>
    dplyr::group_by(.data$Phenotype_2, .data$subject) |>
    dplyr::summarise(frac = mean(.data$found), .groups = "drop")
  pval <- tryCatch(stats::wilcox.test(frac ~ Phenotype_2, data = subj)$p.value,
                   error = function(e) NA_real_)
  n_subj <- table(subj$Phenotype_2)
  grp_order <- intersect(groups_keep, ann$Phenotype_2)
  ex$facet  <- factor(ex$Phenotype_2, levels = grp_order)

  # Relative sizing: normalise each disease group so its column sums to 100%,
  # putting NIU and Viral on the same axis (matches the Fig 6G compartment
  # panel style). Each clone's weight = 100 / (clones in its group).
  n_by_grp <- stats::setNames(ann$n, ann$Phenotype_2)
  ex$w <- 100 / n_by_grp[ex$Phenotype_2]

  # Long (lodes) form: one row per (clone, tissue).
  long <- rbind(
    data.frame(clone_uid = ex$clone_uid, facet = ex$facet, sharing = ex$sharing,
               w = ex$w, tissue = "Blood",
               present = ifelse(ex$found, "In blood", "Eye only"),
               stringsAsFactors = FALSE),
    data.frame(clone_uid = ex$clone_uid, facet = ex$facet, sharing = ex$sharing,
               w = ex$w, tissue = "Eye", present = "In eye",
               stringsAsFactors = FALSE))
  long$tissue  <- factor(long$tissue, levels = c("Blood", "Eye"))
  long$present <- factor(long$present, levels = c("In blood", "Eye only", "In eye"))
  long$sharing <- factor(long$sharing,
                         levels = c("Shared (eye + blood)", "Eye-restricted"))

  # Absolute-count labels on each stratum, positioned by the normalised stack.
  # geom_stratum default reverse=TRUE puts the first present-level on TOP, so at
  # Blood "In blood" sits on top (height = pct), "Eye only" below it.
  lab_df <- do.call(rbind, lapply(seq_len(nrow(ann)), function(i) {
    g <- ann$Phenotype_2[i]; pct <- 100 * ann$n_shared[i] / ann$n[i]
    data.frame(
      facet = factor(g, levels = grp_order),
      tissue = factor(c("Blood", "Blood", "Eye"), levels = c("Blood", "Eye")),
      y = c(100 - pct / 2, (100 - pct) / 2, 50),
      label = c(sprintf("In blood\nn=%d", ann$n_shared[i]),
                sprintf("Eye only\nn=%d", ann$n[i] - ann$n_shared[i]),
                sprintf("In eye\nn=%d", ann$n[i])),
      stringsAsFactors = FALSE)
  }))

  # Bold shared-fraction callout, centred above each facet.
  callout <- data.frame(
    facet = factor(ann$Phenotype_2, levels = grp_order), y = 110,
    label = sprintf("%d%% shared (%d / %d)", ann$pct, ann$n_shared, ann$n),
    stringsAsFactors = FALSE)

  # Highlight shared (indigo) against eye-restricted (sage green), echoing the
  # Fig 6G "Mixed" highlight over "Eye only".
  pal <- c(`Shared (eye + blood)` = "#3B4CA0", `Eye-restricted` = "#A6CE9B")
  g <- ggplot2::ggplot(long,
                       ggplot2::aes(x = .data$tissue, stratum = .data$present,
                                    alluvium = .data$clone_uid, y = .data$w,
                                    fill = .data$sharing)) +
    ggalluvial::geom_flow(width = 0.36, alpha = 0.65,
                          color = "white", linewidth = 0.15) +
    ggalluvial::geom_stratum(width = 0.36, fill = "grey96", color = "black",
                             linewidth = 0.5) +
    ggplot2::geom_text(data = lab_df, inherit.aes = FALSE,
                       ggplot2::aes(x = .data$tissue, y = .data$y, label = .data$label),
                       size = 2.9, lineheight = 0.9) +
    ggplot2::geom_text(data = callout, inherit.aes = FALSE, hjust = 0.5,
                       ggplot2::aes(x = 1.5, y = .data$y, label = .data$label),
                       fontface = "bold", size = 4.1) +
    ggplot2::facet_wrap(~ .data$facet) +
    ggplot2::scale_fill_manual(values = pal, name = NULL, na.translate = FALSE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.3, 0.3))) +
    ggplot2::scale_y_continuous(limits = c(0, 116),
                                breaks = c(0, 25, 50, 75, 100),
                                labels = function(b) paste0(b, "%")) +
    ggplot2::labs(
      title = "Blood-to-eye sharing of expanded TCR clones",
      subtitle = sprintf(paste0("Each group normalised to 100%% of its expanded eye clones ",
                                "(>= %d cells). NIU vs Viral shared fraction: Wilcoxon p = %s."),
                         min_cells, formatC(pval, format = "g", digits = 2)),
      x = NULL, y = "Expanded clones (within-group fraction)") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 8),
                   strip.text = ggplot2::element_text(face = "bold", size = 12),
                   strip.background = ggplot2::element_rect(fill = "grey92",
                                                            color = NA),
                   legend.position = "top")
  save_pdf_png(g, file.path(p$base, "fig5_s3b_eye_blood_clone_alluvial"),
               w = 8.5, h = 6)
  log_message(sprintf("  viz_fig5_s3b: wrote blood-eye clone alluvial (Wilcoxon p=%s).",
                      formatC(pval, format = "g", digits = 2)))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Top-level dispatcher
# ---------------------------------------------------------------------------
run_visualizations_tcr_advanced <- function(cfg) {
  if (!isTRUE(cfg$steps$viz_tcr_advanced)) {
    log_message("viz_tcr_advanced disabled. Skipping.")
    return(invisible(TRUE))
  }
  log_message("Starting viz_tcr_advanced...")
  # Legacy panels (kept for back-compat / supplemental).
  try(viz_gliph_tcrdist_correlation(cfg), silent = FALSE)
  try(viz_gliph_network_metrics(cfg),     silent = FALSE)
  try(viz_vdjdb_overlay(cfg),             silent = FALSE)
  try(viz_antigen_substate_composition(cfg), silent = FALSE)
  try(viz_hla_b27_substate_composition(cfg), silent = FALSE)
  try(viz_tcr_genex_signatures(cfg),      silent = FALSE)
  try(viz_meta_motifs(cfg),               silent = FALSE)
  try(viz_gliph_motif_directional(cfg),   silent = FALSE)
  try(viz_clone_group_antigen_signatures(cfg), silent = FALSE)
  # Figure 5 panels (Acts 1-3).
  try(viz_fig5_a_vdjdb_umap_composite(cfg),    silent = FALSE)
  try(viz_fig5_c_hla_b27_per_subject(cfg),     silent = FALSE)
  try(viz_fig5_d_alluvial_by_antigen_class(cfg), silent = FALSE)
  try(viz_fig5_f_tcrdist_niu_vs_vdjdb(cfg),    silent = FALSE)
  try(viz_fig5_g_gliph_logos_with_b27(cfg),    silent = FALSE)
  try(viz_fig5_h_top_niu_motif_overlap(cfg),   silent = FALSE)
  try(viz_fig5_i_clone_group_lr_with_b27(cfg), silent = FALSE)
  try(viz_fig5h_topn_sensitivity(cfg),         silent = FALSE)
  # Supplemental blood<->eye expanded-clone sharing alluvial.
  try(viz_fig5_s3b_eye_blood_clone_alluvial(cfg), silent = FALSE)
  log_message("viz_tcr_advanced complete.")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Figure 5H top-N sensitivity sweep (folded in from scripts/fig5h_topN_sensitivity.R).
# Diagnostic that the chosen top-N cutoff on the novel-TCR candidate ranking
# does not bias the (cluster, etiology) enrichment. Reads the ranked candidates,
# GLIPH clusters, and the T cell compartment object; writes one CSV and two
# trajectory figures under outputs/viz/eye/tcell/tcr_advanced/.
# ---------------------------------------------------------------------------
viz_fig5h_topn_sensitivity <- function(cfg) {
  cand_path <- "outputs/tables/repertoire/novel_tcr_candidates_ranked.csv"
  cls_path  <- "outputs/tables/repertoire/gliph_clusters.csv"
  obj_path  <- "outputs/objects/eye/tcell/IntegratedSeuratObject.rds"
  if (!all(file.exists(c(cand_path, cls_path, obj_path)))) {
    log_message("  fig5h sensitivity: inputs not on disk; skipping.")
    return(invisible(NULL))
  }
  cand <- read.csv(cand_path, stringsAsFactors = FALSE)
  cls  <- read.csv(cls_path,  stringsAsFactors = FALSE)

  exclude_b27 <- isTRUE(cfg$tcr_advanced$fig5$novel_tcr_discovery$exclude_b27_in_viz)
  if (exclude_b27 && any(cand$is_b27_known, na.rm = TRUE))
    cand <- cand[!cand$is_b27_known, , drop = FALSE]
  cand <- cand[order(cand$composite_rank), , drop = FALSE]

  obj <- readRDS(obj_path)
  md  <- obj@meta.data
  trb <- stringr::str_split(as.character(md$CTaa), "_", simplify = TRUE)
  md$CDR3b <- if (ncol(trb) >= 2) trb[, 2] else NA_character_
  md$knn_cluster <- as.character(md$knn.leiden.cluster)
  md <- md[md$Phenotype_2 == "NIU", , drop = FALSE]
  if ("Tissue_1" %in% colnames(md)) md <- md[md$Tissue_1 == "Eye", , drop = FALSE]
  rm(obj); invisible(gc())

  fisher_per_cat <- function(meta_cells, is_cand, category_col) {
    cat_vec <- meta_cells[[category_col]]
    keep <- !is.na(cat_vec) & cat_vec != ""
    cat_v <- cat_vec[keep]; cand_v <- is_cand[keep]
    cats <- sort(unique(cat_v))
    rows <- lapply(cats, function(ct) {
      a <- sum(cand_v &  cat_v == ct); b <- sum(cand_v &  cat_v != ct)
      c <- sum(!cand_v & cat_v == ct); d <- sum(!cand_v & cat_v != ct)
      ft <- tryCatch(stats::fisher.test(matrix(c(a, c, b, d), nrow = 2)),
                     error = function(e) NULL)
      data.frame(category = ct, n_cand = a, n_baseline = c,
                 OR = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
                 p  = if (!is.null(ft)) ft$p.value          else NA_real_,
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    out$FDR <- stats::p.adjust(out$p, method = "BH")
    out
  }

  N_grid <- sort(unique(c(10, 25, 50, 75, 100, 150, 200, 300, 400, 500,
                          nrow(cand))))
  N_grid <- N_grid[N_grid <= nrow(cand)]
  sweep_results <- list()
  for (N in N_grid) {
    top_ids <- cand$cluster_id[seq_len(N)]
    cand_cdr3 <- unique(cls$CDR3b[cls$cluster_id %in% top_ids])
    is_cand <- !is.na(md$CDR3b) & md$CDR3b != "" & md$CDR3b %in% cand_cdr3
    if (!any(is_cand)) next
    fc <- fisher_per_cat(md, is_cand, "knn_cluster"); fc$N <- N; fc$type <- "cluster"
    fe <- fisher_per_cat(md, is_cand, "Etiology");    fe$N <- N; fe$type <- "etiology"
    sweep_results[[as.character(N)]] <- rbind(fc, fe)
  }
  if (length(sweep_results) == 0L) {
    log_message("  fig5h sensitivity: no candidate cells matched; skipping.")
    return(invisible(NULL))
  }
  sweep <- dplyr::bind_rows(sweep_results)

  ensure_dir(file.path(cfg$paths$results_tables, "repertoire"))
  write.csv(sweep, file.path(cfg$paths$results_tables,
            "repertoire", "fig5h_topN_sensitivity.csv"), row.names = FALSE)

  plot_dir <- .tcra_paths(cfg)$base
  ensure_dir(plot_dir)
  pal_cat <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
               "#8c564b","#e377c2","#17becf","#bcbd22","#7f7f7f")
  mk_traj <- function(d, title) {
    d$log2_OR <- log2(pmax(d$OR, 1e-6))
    d$sig <- !is.na(d$FDR) & d$FDR < 0.05
    cats <- sort(unique(d$category))
    pal <- setNames(pal_cat[seq_along(cats)], cats)
    ggplot(d, aes(x = N, y = log2_OR, color = category, group = category)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
      geom_line(linewidth = 0.7, alpha = 0.85) +
      geom_point(aes(shape = sig), size = 2.4, stroke = 0.6) +
      scale_color_manual(values = pal, name = title) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                         labels = c(`TRUE` = "FDR<0.05", `FALSE` = "n.s."), name = NULL) +
      scale_x_log10(breaks = N_grid) +
      labs(title = sprintf("Top-N sensitivity: log2(OR) by %s", tolower(title)),
           subtitle = "Cand vs baseline (NIU eye T). Flat lines = signal stable across cutoff.",
           x = "Top-N candidate motifs (log scale)", y = "log2(odds ratio)") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(size = 9),
            panel.grid.minor = element_blank())
  }
  for (ty in c("cluster", "etiology")) {
    lab <- if (ty == "cluster") "Cluster" else "Indication"
    p <- mk_traj(sweep[sweep$type == ty, ], lab)
    ggsave(file.path(plot_dir, paste0("fig5h_topN_sensitivity_", ty, ".pdf")),
           p, width = 8, height = 5)
    ggsave(file.path(plot_dir, paste0("fig5h_topN_sensitivity_", ty, ".png")),
           p, width = 8, height = 5, dpi = 300)
  }
  log_message("  fig5h sensitivity: wrote CSV + trajectory figures.")
  invisible(sweep)
}
