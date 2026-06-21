# R/83_viz_full.R
# Figure 1 (full atlas) visualizations + F1-only viz_* helpers.
# Shared viz_* functions used by both F1 and F2 live in 82_viz_dispatch.R.
suppressPackageStartupMessages({
  library(Seurat)
  library(scplotter)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(patchwork)
})

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

# ============================================================================
# COMPARTMENT CONTRIBUTION TO THE NIU-vs-VIRAL BULK SIGNAL
# ============================================================================
# For each leading program gene, what fraction of its total counts each lineage
# contributes (abundance-weighted decomposition of the pseudobulk contrast).
# Count-accurate, run in-pipeline off the full integrated object — replaces the
# former scripts/fig2g_*.R standalones. Two curated programs are rendered:
#   antiviral effector (up in Viral) and autoimmune (up in NIU), plus a combined
#   panel. Genes absent from the object are dropped.
viz_program_contribution <- function(obj, cfg, paths = cfg$paths,
                                     cluster_col = "merged.celltype.cluster",
                                     assay = "RNA") {
  log_message("Visualizing: compartment contribution (antiviral / autoimmune programs)")
  if (!cluster_col %in% colnames(obj[[]])) {
    log_message("  ", cluster_col, " absent; skipping program contribution.")
    return(invisible(NULL))
  }
  vdir <- file.path(paths$viz_dir, "program_contribution")
  ensure_dir(vdir)

  submodules_viral <- list(
    ISG               = c("IFI27","IFI6","IFI44L","IFIT2","IFIT3","ISG15","MX1","RSAD2"),
    `CXCR3 chemokine` = c("CXCL9","CXCL10","CXCL11"),
    `Mono chemokine`  = c("CCL2","CCL7","CCL8","CCL18"),
    Complement        = c("C1QC","CFB","SERPING1"),
    Checkpoint        = c("LAG3","IL10","CD274"))
  submodules_auto <- list(
    `DC / cDC2`      = c("FCER1A","CD1C","CD1E","CD1A","CLEC10A"),
    `Classical mono` = c("FCN1","S100A12","FOLR3","CD300H","ADGRE1","CTSG","VCAN","SELL"),
    `Th17 / IL-23`   = c("RORC","IL23A","CCR6","IL22"),
    `T costim`       = c("TNFRSF8","CD40LG"),
    `B lineage`      = c("CD79B","MS4A1","TCL1A"))
  submod_levels <- c(names(submodules_viral), names(submodules_auto))

  lineage_levels <- c("Myeloid","CD8 T","CD4 T","Other T","NK","B/Plasma","Mixed")
  lineage_pal <- c(Myeloid="#E4572E","CD8 T"="#2E86AB","CD4 T"="#5BC0BE",
                   "Other T"="#9BC53D", NK="#7B2CBF","B/Plasma"="#F2A65A",
                   Mixed="#B0B0B0")
  lineage_of <- function(x) dplyr::case_when(
    grepl("Mono|Macro|DC|cDC|pDC|Langerhans|Myeloid", x) ~ "Myeloid",
    grepl("^CD8", x)                                     ~ "CD8 T",
    grepl("^CD4|Treg", x)                                ~ "CD4 T",
    grepl("gdT|MAIT|dnT", x)                             ~ "Other T",
    grepl("NK", x)                                       ~ "NK",
    grepl("B/Plasma|Plasma|^B ", x)                      ~ "B/Plasma",
    TRUE                                                 ~ "Mixed")

  present <- rownames(obj)
  build_long <- function(submods) {
    submods <- lapply(submods, intersect, present)
    submods <- submods[lengths(submods) >= 1]
    if (length(submods) == 0) return(NULL)
    g2m <- utils::stack(submods) |>
      dplyr::transmute(gene = as.character(values), submodule = as.character(ind))
    genes <- unique(g2m$gene)
    ord <- genes[order(match(g2m$submodule[match(genes, g2m$gene)], names(submods)))]
    ag <- tryCatch(
      AggregateExpression(obj, features = genes, group.by = cluster_col,
                          assays = assay, layer = "counts")[[assay]],
      error = function(e)
        AggregateExpression(obj, features = genes, group.by = cluster_col,
                            assays = assay, slot = "counts")[[assay]])
    ag <- as.matrix(ag)
    cb <- sweep(ag, 1, pmax(rowSums(ag), 1), "/")
    tibble::as_tibble(cb, rownames = "gene") |>
      tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "frac") |>
      dplyr::mutate(lineage = factor(lineage_of(cluster), lineage_levels)) |>
      dplyr::left_join(g2m, by = "gene") |>
      dplyr::group_by(gene, submodule, lineage) |>
      dplyr::summarise(frac = sum(frac), .groups = "drop") |>
      dplyr::mutate(submodule = factor(submodule, submod_levels),
                    gene = factor(gene, rev(ord)))
  }

  panel <- function(d, ttl) {
    ggplot(d, aes(frac, gene, fill = lineage)) +
      geom_col(width = 0.8, color = "white", linewidth = 0.2) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
      scale_fill_manual(values = lineage_pal, name = "Lineage", drop = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
      facet_grid(submodule ~ ., scales = "free_y", space = "free_y", switch = "y") +
      theme_minimal(base_size = 10) +
      theme(strip.text.y.left = element_text(angle = 0, face = "bold"),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank()) +
      labs(title = ttl, x = "Fraction of total expression", y = NULL)
  }

  viral <- build_long(submodules_viral)
  auto  <- build_long(submodules_auto)
  if (is.null(viral) && is.null(auto)) {
    log_message("  No program genes present in object; skipping.")
    return(invisible(NULL))
  }

  viral_ttl <- "Antiviral effector program: compartment contribution (up in Viral)"
  auto_ttl  <- "Autoimmune program: compartment contribution (up in NIU)"

  if (!is.null(viral)) {
    ggsave(file.path(vdir, "antiviral_program_compartment_contribution.pdf"),
           panel(viral, viral_ttl), device = grDevices::cairo_pdf,
           width = 8, height = 0.3 * nlevels(droplevels(viral$gene)) + 2)
    utils::write.csv(viral, file.path(paths$results_tables,
                     "antiviral_program_contribution.csv"), row.names = FALSE)
  }
  if (!is.null(auto)) {
    ggsave(file.path(vdir, "autoimmune_program_compartment_contribution.pdf"),
           panel(auto, auto_ttl), device = grDevices::cairo_pdf,
           width = 8, height = 0.3 * nlevels(droplevels(auto$gene)) + 2)
    utils::write.csv(auto, file.path(paths$results_tables,
                     "autoimmune_program_contribution.csv"), row.names = FALSE)
  }
  if (!is.null(viral) && !is.null(auto)) {
    combined <- (panel(viral, viral_ttl) / panel(auto, auto_ttl)) +
      plot_layout(guides = "collect",
                  heights = c(nlevels(droplevels(viral$gene)),
                              nlevels(droplevels(auto$gene)))) +
      plot_annotation(
        title = "Compartment contribution to the NIU vs Viral bulk signal",
        subtitle = "Fraction of each gene's total counts contributed by each lineage")
    ggsave(file.path(vdir, "niu_vs_viral_program_contribution_combined.pdf"),
           combined, device = grDevices::cairo_pdf, width = 8.5,
           height = 0.26 * (nlevels(droplevels(viral$gene)) +
                            nlevels(droplevels(auto$gene))) + 2.5)
  }
  log_message("  Wrote program-contribution panels to ", vdir)
  invisible(TRUE)
}

# ============================================================================
# F1 ENTRY POINT
# ============================================================================
# run_visualizations_full coordinates Figure 1 panels: shared viz_* helpers
# (defined in 82_viz_dispatch.R) plus the F1-only repertoire / BCR-lineage
# blocks above. Reads the full atlas IntegratedSeuratObject.rds.
run_visualizations_full <- function(cfg) {
  paths <- get_target_paths(cfg, "all")
  log_message("=== Figure 1 (full atlas) visualizations ===")
  ensure_dir(paths$viz_dir)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Full atlas IntegratedSeuratObject.rds not found. Skipping F1 viz.")
    return(invisible(TRUE))
  }
  obj <- readRDS(obj_path)

  # Shared blocks
  viz_qc_summary(obj, cfg, paths)
  viz_integration(obj, cfg, paths)
  viz_celltypes(obj, cfg, paths)
  viz_markers(obj, cfg, paths)
  viz_dge(obj, cfg, paths, target = "all")
  viz_milo(obj, cfg, paths, target = "all")
  viz_composition(cfg, paths, target = "all")
  viz_program_contribution(obj, cfg, paths)

  # F1-only blocks
  if (isTRUE(cfg$steps$bcr_lineage)) viz_bcr_lineage(cfg, paths)

  log_message("=== Figure 1 visualizations complete ===")
  invisible(TRUE)
}
