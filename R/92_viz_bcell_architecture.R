# R/92_viz_bcell_architecture.R
# Plots for Figure 6 panels F-i, F-ii, F-iii, G-i, G-ii.
# Reads CSVs written by R/58_bcell_lineage_architecture.R.
# One ggplot per panel, one PDF per panel.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(ggalluvial)
})

.disease_palette <- c(Viral = "#3B6FB6", Autoimmune = "#C0362C")

.disease_factor <- function(x) {
  factor(dplyr::recode(x, NIU = "Autoimmune"),
         levels = c("Viral", "Autoimmune"))
}

.theme_fig6 <- function() {
  ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

.write_pdf <- function(p, path, w = 3.0, h = 2.8) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(path, p, width = w, height = h, units = "in", device = "pdf")
}

run_viz_bcell_architecture <- function(cfg) {
  paths   <- get_target_paths(cfg, "bcell")
  set.seed(cfg$seed %||% 1L)
  # paths$results_tables and paths$viz_dir already end in "bcell"; mirror the
  # R/58 convention and write under <results_tables>/architecture and
  # <viz_dir>/10_lineage_arch/architecture.
  in_dir  <- file.path(paths$results_tables, "architecture")
  out_dir <- file.path(viz_subdir(paths, "lineage_arch"), "architecture")

  pc_path <- file.path(in_dir, "clone_architecture_metrics.csv")
  if (!file.exists(pc_path)) {
    log_message("[fig6 viz F/G] missing inputs at ", pc_path, "; skipping.")
    return(invisible(FALSE))
  }
  log_message("[fig6 viz F/G] reading ", pc_path)
  per_clone <- readr::read_csv(pc_path, show_col_types = FALSE) %>%
    dplyr::filter(phenotype %in% c("Viral", "NIU")) %>%
    dplyr::mutate(phenotype = .disease_factor(phenotype))

  # ---- Panel F-i: clusters spanned per clone -----------------------------
  p_fi <- ggplot2::ggplot(
    per_clone, ggplot2::aes(x = phenotype, y = n_clusters,
                            fill = phenotype, colour = phenotype)) +
    ggplot2::geom_violin(alpha = 0.35, scale = "width",
                         draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_jitter(width = 0.12, size = 0.4, alpha = 0.4) +
    ggplot2::scale_fill_manual(values = .disease_palette) +
    ggplot2::scale_colour_manual(values = .disease_palette) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    ggplot2::labs(x = NULL, y = "B-cell clusters per clone") +
    .theme_fig6()
  .write_pdf(p_fi, file.path(out_dir, "fig6F_i_clusters_spanned.pdf"))

  # ---- Panel F-ii: SHM by shared/non-shared, faceted by phenotype --------
  # Mirrors the S6c layout: 2 boxes per facet (Non-shared vs Shared) on
  # per-clone shm_mean, faceted by Viral / Autoimmune. Per-disease
  # Wilcoxon contrasts live in stats_summary.csv.
  fii_dat <- per_clone %>%
    dplyr::filter(!is.na(shm_mean)) %>%
    dplyr::mutate(shared_label = factor(
      ifelse(is_shared, "Shared", "Non-shared"),
      levels = c("Non-shared", "Shared")))
  p_fii <- ggplot2::ggplot(
    fii_dat,
    ggplot2::aes(x = shared_label, y = shm_mean,
                 fill = shared_label, colour = shared_label)) +
    ggplot2::geom_boxplot(width = 0.5, alpha = 0.3,
                          outlier.shape = NA, linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.15, size = 0.4, alpha = 0.45) +
    ggplot2::facet_wrap(~ phenotype, nrow = 1) +
    ggplot2::scale_fill_manual(values = c(`Non-shared` = "#888888",
                                          Shared = "#3B6FB6")) +
    ggplot2::scale_colour_manual(values = c(`Non-shared` = "#555555",
                                            Shared = "#264E84")) +
    ggplot2::labs(x = NULL, y = "Mean SHM frequency per clone") +
    .theme_fig6() +
    ggplot2::theme(legend.position = "none")
  .write_pdf(p_fii, file.path(out_dir, "fig6F_ii_shm_load.pdf"),
             w = 4.0, h = 2.8)

  # ---- Panel F-iii: CDR-R selection ratio (log10), faceted by n_clusters --
  fiii_dat <- per_clone %>%
    dplyr::filter(!is.na(sel_cdr_r_ratio), sel_cdr_r_ratio > 0)
  # Keep facets with >=5 clones in both diseases; drop sparse ones (<5).
  facet_ok <- fiii_dat %>%
    dplyr::count(phenotype, n_clusters) %>%
    tidyr::pivot_wider(names_from = phenotype, values_from = n,
                       values_fill = 0) %>%
    dplyr::filter(.data[["Viral"]] >= 5, .data[["Autoimmune"]] >= 5) %>%
    dplyr::pull(n_clusters)
  log_message("[fig6 viz F-iii] kept n_clusters facets (>=5 in both): ",
              paste(facet_ok, collapse = ", "))
  fiii_dat <- dplyr::filter(fiii_dat, n_clusters %in% facet_ok)

  p_fiii <- ggplot2::ggplot(
    fiii_dat,
    ggplot2::aes(x = phenotype, y = sel_cdr_r_ratio,
                 fill = phenotype, colour = phenotype)) +
    ggplot2::geom_violin(alpha = 0.35, scale = "width",
                         draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_jitter(width = 0.12, size = 0.4, alpha = 0.35) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
    ggplot2::facet_wrap(~ n_clusters, nrow = 1,
                        labeller = ggplot2::labeller(
                          n_clusters = function(x) paste0(x, " clusters"))) +
    ggplot2::scale_fill_manual(values = .disease_palette) +
    ggplot2::scale_colour_manual(values = .disease_palette) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL,
                  y = "CDR-R selection ratio (log10)") +
    .theme_fig6()
  .write_pdf(p_fiii, file.path(out_dir, "fig6F_iii_selection_ratio.pdf"),
             w = 4.2, h = 2.8)

  # ---- Panel G-i: shared-compartment alluvial flow (per-phenotype frac) ---
  # Two-column stacked proportion bar: x = phenotype (Viral, Autoimmune),
  # within-column fill = comp_class (eye_only / blood_only / mixed),
  # connected by translucent ggalluvial flow ribbons so same-stratum
  # proportions read continuously across the two diseases. Operates on
  # the *shared* clone subset to mirror the design's "shared compartment"
  # framing.
  shared_gi <- per_clone %>% dplyr::filter(is_shared)
  if (nrow(shared_gi) == 0L) {
    log_message("[fig6 viz G-i] no shared clones; skipping panel G-i.")
  } else {
    comp_levels <- c("eye_only", "blood_only", "mixed")
    # tidyr::complete pads (phenotype, comp_class) so missing strata get
    # zero-height ribbons that taper visibly to zero on the empty side.
    # Without this, geom_alluvium has no lode to connect and the ribbon
    # is invisible (e.g. all Autoimmune shared clones are eye_only).
    g_i_df <- shared_gi %>%
      dplyr::filter(!is.na(comp_class)) %>%
      dplyr::count(phenotype, comp_class) %>%
      tidyr::complete(phenotype, comp_class = comp_levels,
                      fill = list(n = 0L)) %>%
      dplyr::group_by(phenotype) %>%
      dplyr::mutate(frac = n / sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(comp_class = factor(comp_class, levels = comp_levels))

    p_gi <- ggplot2::ggplot(
      g_i_df,
      ggplot2::aes(x = phenotype, y = frac,
                   stratum = comp_class, alluvium = comp_class,
                   fill = comp_class)) +
      ggalluvial::geom_alluvium(alpha = 0.5, knot.pos = 0.25) +
      ggalluvial::geom_stratum(width = 0.4, colour = "black",
                               linewidth = 0.4) +
      ggplot2::scale_fill_manual(
        values = c(eye_only = "#7AA86F",
                   blood_only = "#C97B3F",
                   mixed = "#4E5BAA"),
        labels = c(eye_only = "Eye only",
                   blood_only = "Blood only",
                   mixed = "Mixed (eye + blood)")) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = c(0, 0)) +
      ggplot2::labs(x = NULL,
                    y = "Shared clones (within-column fraction)",
                    fill = NULL) +
      .theme_fig6()
    .write_pdf(p_gi, file.path(out_dir, "fig6G_i_shared_compartment.pdf"))
  }

  # ---- Panel G-ii: eye-blood Shannon entropy of shared clones -------------
  if (nrow(shared_gi) > 0) {
    gii_dat <- shared_gi %>% dplyr::filter(!is.na(tissue_entropy))
    if (nrow(gii_dat) > 0) {
      p_gii <- ggplot2::ggplot(
        gii_dat,
        ggplot2::aes(x = phenotype, y = tissue_entropy,
                     fill = phenotype, colour = phenotype)) +
        ggplot2::geom_violin(alpha = 0.35, scale = "width",
                             draw_quantiles = c(0.25, 0.5, 0.75)) +
        ggplot2::geom_jitter(width = 0.12, size = 0.5, alpha = 0.45) +
        ggplot2::scale_fill_manual(values = .disease_palette) +
        ggplot2::scale_colour_manual(values = .disease_palette) +
        ggplot2::labs(x = NULL,
                      y = "Eye-blood Shannon entropy") +
        .theme_fig6()
      .write_pdf(p_gii, file.path(out_dir, "fig6G_ii_shared_tissue_entropy.pdf"))
    } else {
      log_message("[fig6 viz G-ii] no non-NA tissue_entropy rows; skipping panel G-ii.")
    }
  } else {
    log_message("[fig6 viz G-ii] no shared clones; skipping panel G-ii.")
  }

  log_message("[fig6 viz F/G] done.")
  invisible(TRUE)
}
