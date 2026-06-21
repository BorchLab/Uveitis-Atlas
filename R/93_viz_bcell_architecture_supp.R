# R/93_viz_bcell_architecture_supp.R
# Fig S6 supplement panels (BASELINe, public clones, isotype dynamics,
# tree topology). One ggplot per panel, one PDF per panel.
# Reads CSVs written by R/58_bcell_lineage_architecture.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(ggalluvial)
})

# Helpers duplicated from R/76 so this script loads in isolation
# (mirrors the R/74b precedent — duplication over coupling for ~5 small constants).
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

run_viz_bcell_architecture_supp <- function(cfg) {
  paths   <- get_target_paths(cfg, "bcell")
  set.seed(cfg$seed %||% 1L)
  # Mirror R/76's deviation: paths$results_tables and paths$viz_dir already
  # end in "bcell"; write architecture-supplement outputs under
  # <results_tables>/architecture/ and <viz_dir>/10_lineage_arch/architecture_supp/.
  in_dir  <- file.path(paths$results_tables, "architecture")
  out_dir <- file.path(viz_subdir(paths, "lineage_arch"), "architecture_supp")

  pc_path <- file.path(in_dir, "clone_architecture_metrics.csv")
  if (!file.exists(pc_path)) {
    log_message("[figS6] missing per-clone metrics at ", pc_path,
                "; skipping all supplement panels.")
    return(invisible(FALSE))
  }
  log_message("[figS6] reading ", pc_path)
  per_clone <- readr::read_csv(pc_path, show_col_types = FALSE) %>%
    dplyr::filter(phenotype %in% c("Viral", "NIU"))

  # ---- Panel S6a: BASELINe selection sigma per subject x substate x region
  # R/59 writes baseline_selection.csv with columns:
  #   subject, phenotype, substate, region (cdr/fwr), baseline_sigma,
  #   baseline_ci_lower, baseline_ci_upper, baseline_ci_pvalue.
  # Plot baseline_sigma as violin + jitter per disease, faceted by region
  # (rows) and substate (cols), with a dashed zero line for neutral
  # selection. Skips gracefully if R/59 hasn't been run.
  bl_path <- file.path(in_dir, "baseline_selection.csv")
  if (file.exists(bl_path)) {
    log_message("[figS6a] reading BASELINe summary from ", bl_path)
    bl <- readr::read_csv(bl_path, show_col_types = FALSE) %>%
      dplyr::filter(phenotype %in% c("Viral", "NIU"),
                    !is.na(baseline_sigma), is.finite(baseline_sigma)) %>%
      dplyr::mutate(
        phenotype = .disease_factor(phenotype),
        region    = factor(toupper(region), levels = c("CDR", "FWR")))

    if (nrow(bl) > 0L) {
      p_s6a <- ggplot2::ggplot(
        bl, ggplot2::aes(x = phenotype, y = baseline_sigma,
                         fill = phenotype, colour = phenotype)) +
        ggplot2::geom_violin(alpha = 0.35, scale = "width",
                             draw_quantiles = 0.5) +
        ggplot2::geom_jitter(width = 0.12, size = 0.8, alpha = 0.85) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2,
                            colour = "grey60") +
        ggplot2::facet_grid(region ~ substate,
                            labeller = ggplot2::labeller(
                              substate = ggplot2::label_wrap_gen(width = 14))) +
        ggplot2::scale_fill_manual(values = .disease_palette) +
        ggplot2::scale_colour_manual(values = .disease_palette) +
        ggplot2::labs(x = NULL,
                      y = expression(BASELINe~sigma~"(focused)")) +
        .theme_fig6()
      .write_pdf(p_s6a, file.path(out_dir, "figS6_a_baseline_selection.pdf"),
                 w = 5)
    } else {
      log_message("[figS6a] no Viral/NIU sigma rows; skipping.")
    }
  } else {
    log_message("[figS6a] baseline_selection.csv not at ", bl_path,
                "; skipping (run R/59).")
  }

  # ---- Panel S6b: public clones per subject -------------------------------
  # R/53 writes outputs/tables/repertoire/BCR_public_clones.csv. Columns:
  #   cell_id_unique, sequence_id, subject_id, sample_id, v_gene, j_gene,
  #   cdr3_len, junction_aa, c_call, cluster_id, stratum.
  # We collapse to unique (subject_id, cluster_id) pairs to count public
  # clones per subject (cluster_id is the shared public-clone identifier).
  paths_rep <- tryCatch(get_target_paths(cfg, "repertoire"),
                        error = function(e) NULL)
  pub_path <- if (!is.null(paths_rep) &&
                  !is.null(paths_rep$results_tables)) {
    file.path(paths_rep$results_tables, "BCR_public_clones.csv")
  } else {
    file.path("outputs", "tables", "repertoire", "BCR_public_clones.csv")
  }
  if (file.exists(pub_path)) {
    log_message("[figS6b] reading public clones from ", pub_path)
    pub <- readr::read_csv(pub_path, show_col_types = FALSE)
    if (all(c("subject_id", "cluster_id", "stratum") %in% colnames(pub))) {
      subj_pub <- pub %>%
        dplyr::filter(stratum %in% c("Viral", "NIU")) %>%
        dplyr::distinct(subject_id, cluster_id, stratum) %>%
        dplyr::count(subject_id, stratum, name = "n_public") %>%
        dplyr::rename(phenotype = stratum) %>%
        dplyr::mutate(phenotype = .disease_factor(phenotype))
      if (nrow(subj_pub) > 0) {
        # Wilcoxon + Cliff's delta annotation on per-subject counts.
        subtitle <- NULL
        if (dplyr::n_distinct(subj_pub$phenotype) >= 2L &&
            nrow(subj_pub) >= 4L) {
          w <- tryCatch(
            stats::wilcox.test(n_public ~ phenotype, data = subj_pub,
                               exact = FALSE),
            error = function(e) NULL)
          cd <- NA_real_
          if (requireNamespace("effsize", quietly = TRUE)) {
            cd <- tryCatch(
              effsize::cliff.delta(n_public ~ phenotype,
                                   data = subj_pub)$estimate,
              error = function(e) NA_real_)
          }
          if (!is.null(w)) {
            subtitle <- sprintf(
              "W = %.0f, p = %.2g, Cliff's delta = %.2f",
              unname(w$statistic), w$p.value, cd)
          }
        }
        p_s6b <- ggplot2::ggplot(
          subj_pub, ggplot2::aes(x = phenotype, y = n_public,
                                 fill = phenotype, colour = phenotype)) +
          ggplot2::geom_boxplot(width = 0.4, alpha = 0.25,
                                outlier.shape = NA) +
          ggplot2::geom_jitter(width = 0.12, size = 1.2, alpha = 0.85) +
          ggplot2::scale_fill_manual(values = .disease_palette) +
          ggplot2::scale_colour_manual(values = .disease_palette) +
          ggplot2::labs(x = NULL, y = "Public clones per subject",
                        subtitle = subtitle) +
          .theme_fig6()
        .write_pdf(p_s6b, file.path(out_dir, "figS6_b_public_clones.pdf"))
      } else {
        log_message("[figS6b] no Viral/NIU subjects in public table; skipping.")
      }
    } else {
      log_message("[figS6b] BCR_public_clones.csv missing expected ",
                  "columns (subject_id, cluster_id, stratum); skipping.")
    }
  } else {
    log_message("[figS6b] public clone CSV not at ", pub_path, "; skipping.")
  }

  # ---- Panel S6c: shared vs non-shared x isotype alluvial flow ------------
  # Within-column stacked proportion bars (Shared / Non-shared) of clone
  # dominant isotype, connected by translucent ggalluvial flow ribbons,
  # faceted by disease. Operates on per_clone (one row per clone) so the
  # "Shared" label matches the lineage architecture definition (n_clusters
  # >= min_clusters_for_shared).
  have_gg <- requireNamespace("ggalluvial", quietly = TRUE)
  if (have_gg) {
    suppressPackageStartupMessages(
      tryCatch({
        if (!"package:ggalluvial" %in% search())
          attachNamespace("ggalluvial")
        TRUE
      }, error = function(e) FALSE))
  }
  if (!have_gg) {
    log_message("[figS6c] ggalluvial not installed; skipping panel S6c.")
  } else if (!"dominant_isotype" %in% colnames(per_clone)) {
    log_message("[figS6c] dominant_isotype column absent; skipping panel S6c.")
  } else {
    iso_levels <- c("IgM", "IgD", "IgG", "IgA", "IgE", "Other")
    # tidyr::complete pads (phenotype, shared_label, iso) so missing isotypes
    # in either Shared or Non-shared column get zero-height ribbons that
    # taper visibly. Without this, geom_alluvium has no matching lode and
    # the ribbon is invisible.
    iso_df <- per_clone %>%
      dplyr::filter(!is.na(dominant_isotype)) %>%
      dplyr::mutate(
        iso = dplyr::case_when(
          grepl("^IGHM", dominant_isotype) ~ "IgM",
          grepl("^IGHD", dominant_isotype) ~ "IgD",
          grepl("^IGHG", dominant_isotype) ~ "IgG",
          grepl("^IGHA", dominant_isotype) ~ "IgA",
          grepl("^IGHE", dominant_isotype) ~ "IgE",
          TRUE                              ~ "Other"),
        shared_label = factor(ifelse(is_shared, "Shared", "Non-shared"),
                              levels = c("Non-shared", "Shared")),
        phenotype = .disease_factor(phenotype)) %>%
      dplyr::count(phenotype, shared_label, iso, name = "n") %>%
      tidyr::complete(phenotype, shared_label, iso = iso_levels,
                      fill = list(n = 0L)) %>%
      dplyr::group_by(phenotype, shared_label) %>%
      dplyr::mutate(frac = n / sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(iso = factor(iso, levels = iso_levels))

    if (nrow(iso_df) == 0L) {
      log_message("[figS6c] no clones with non-NA dominant_isotype; skipping.")
    } else {
      p_s6c <- ggplot2::ggplot(
        iso_df,
        ggplot2::aes(x = shared_label, y = frac,
                     stratum = iso, alluvium = iso, fill = iso)) +
        ggalluvial::geom_alluvium(alpha = 0.5, knot.pos = 0.25) +
        ggalluvial::geom_stratum(width = 0.4, colour = "black",
                                 linewidth = 0.4) +
        ggplot2::facet_wrap(~ phenotype, nrow = 1) +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                    expand = c(0, 0)) +
        ggplot2::labs(x = NULL, y = "Clones (within-column fraction)",
                      fill = "Isotype") +
        .theme_fig6() +
        ggplot2::theme(legend.title = ggplot2::element_text())
      .write_pdf(p_s6c,
                 file.path(out_dir, "figS6_c_isotype_dynamics.pdf"),
                 w = 4.5)
    }
  }

  # ---- Panel S6d: tree topology (Sackin + max depth) ----------------------
  # Cached tree object is a NAMED LIST of dowser tibbles, keyed
  # "<subject>|<clone_id>", each row holding one ape::phylo in $trees[[1]].
  # No phenotype column on the tibble; we join phenotype from the
  # per-clone metrics (which carries (subject, clone_id, phenotype)).
  tree_path <- file.path("outputs", "objects", "ibex",
                         "bcell_lineage_top_trees.rds")
  if (file.exists(tree_path) &&
      requireNamespace("ape", quietly = TRUE) &&
      requireNamespace("purrr", quietly = TRUE)) {
    log_message("[figS6d] reading tree cache ", tree_path)
    trees <- readRDS(tree_path)
    if (is.list(trees) && length(trees) > 0L) {
      pheno_map <- per_clone %>%
        dplyr::mutate(clone_id = as.character(clone_id)) %>%
        dplyr::select(subject, clone_id, phenotype) %>%
        dplyr::distinct()

      topo <- purrr::map_dfr(names(trees), function(key) {
        tib <- trees[[key]]
        if (!is.data.frame(tib) || nrow(tib) == 0L) return(NULL)
        phy <- tib$trees[[1]]
        if (is.null(phy) || !inherits(phy, "phylo")) return(NULL)
        parts <- strsplit(key, "|", fixed = TRUE)[[1]]
        if (length(parts) != 2L) return(NULL)
        depths <- ape::node.depth.edgelength(phy)
        n_tip  <- ape::Ntip(phy)
        # Sackin imbalance proxy: sum of root-to-leaf depths over all tips.
        # On uniform branch lengths this is the classic Sackin index; with
        # real branch lengths it is the cumulative tip depth — still a valid
        # between-cohort comparison.
        tibble::tibble(
          subject   = parts[1],
          clone_id  = parts[2],
          sackin    = sum(depths[seq_len(n_tip)]),
          max_depth = max(depths)
        )
      })

      if (!is.null(topo) && nrow(topo) > 0L) {
        topo <- topo %>%
          dplyr::left_join(pheno_map, by = c("subject", "clone_id")) %>%
          dplyr::filter(phenotype %in% c("Viral", "NIU")) %>%
          dplyr::mutate(phenotype = .disease_factor(phenotype))

        if (nrow(topo) > 0L) {
          topo_long <- topo %>%
            tidyr::pivot_longer(c(sackin, max_depth),
                                names_to = "metric", values_to = "value") %>%
            dplyr::filter(is.finite(value))
          p_s6d <- ggplot2::ggplot(
            topo_long,
            ggplot2::aes(x = phenotype, y = value,
                         fill = phenotype, colour = phenotype)) +
            ggplot2::geom_violin(alpha = 0.35, scale = "width",
                                 draw_quantiles = c(0.25, 0.5, 0.75)) +
            ggplot2::geom_jitter(width = 0.12, size = 0.8, alpha = 0.8) +
            ggplot2::facet_wrap(
              ~ metric, scales = "free_y",
              labeller = ggplot2::labeller(
                metric = c(sackin    = "Sackin imbalance",
                           max_depth = "Max root-to-tip depth"))) +
            ggplot2::scale_fill_manual(values = .disease_palette) +
            ggplot2::scale_colour_manual(values = .disease_palette) +
            ggplot2::labs(x = NULL, y = NULL) +
            .theme_fig6()
          .write_pdf(p_s6d,
                     file.path(out_dir, "figS6_d_tree_topology.pdf"),
                     w = 4.5)
        } else {
          log_message("[figS6d] no Viral/NIU trees after phenotype join; skipping.")
        }
      } else {
        log_message("[figS6d] tree cache produced 0 topology rows; skipping.")
      }
    } else {
      log_message("[figS6d] tree cache empty; skipping.")
    }
  } else {
    log_message("[figS6d] tree cache missing or ape unavailable; skipping.")
  }

  log_message("[figS6] done.")
  invisible(TRUE)
}
