# R/91_viz_cross_compartment.R
# Cross-compartment viz for Figure 4 panels F-H.
#
# Panel F  — Subject-level myeloid PC1 vs T cell PC1 scatter with within-
#            etiology regression, bootstrap CI on the legend (weighted
#            primary, unweighted as supplement).
# Panel G  — LIANA chord diagram (NIU + Viral side-by-side) for top-N LR
#            pairs by aggregate_rank in each condition.
# Panel G2 — Companion dot plot: y = ligand-receptor, x = source x target,
#            dot size = -log10(consensus_rank), dot color = disease_bias.
# Panel H  — NicheNet heatmap of top ligands x myeloid substates per T cell
#            substate x pole (only rendered when nichenet step ran).
#
# Driver run_visualizations_cross_compartment(cfg) is dispatched by
# R/82_viz_dispatch.R when target = "cross_compartment".
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

.cc_paths <- function(cfg) {
  cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/",
    viz    = "outputs/viz/cross_compartment/")
}

# Panel F driver. Reads pc1_bridge_scores.csv + pc1_bridge_correlation.csv
# and writes three PDFs per weighting plus a per-substate heatmap:
#   myeloid_tcell_pc1_within_etiology_<w>     primary panel - faceted within-etiology scatter
#                                with subject labels and per-arm Pearson r.
#                                This is the honest within-patient claim.
#   myeloid_tcell_pc1_partial_correlation_<w>             residual-residual scatter after regressing
#                                each compartment's PC1 on Phenotype_2.
#                                Partial-r controlling for etiology.
#   myeloid_tcell_pc1_pooled_correlation_<w>   pooled overplot kept as supplement so the
#                                between-group separation is visible but not
#                                presented as the headline.
#   myeloid_tcell_pc1_per_substate_coupling_heatmap  per-(myeloid x tcell) substate
#                                partial-r heatmap (reads pc1_bridge_per_substate.csv).
viz_cross_compartment_bridge <- function(cfg) {
  cc <- .cc_paths(cfg)
  ensure_dir(cc$viz)
  scores_csv <- file.path(cc$tables, "pc1_bridge_scores.csv")
  corr_csv   <- file.path(cc$tables, "pc1_bridge_correlation.csv")
  if (!file.exists(scores_csv) || !file.exists(corr_csv)) {
    log_message("viz_cross_compartment_bridge: scores or correlation CSV ",
                "missing; run cross_compartment_bridge first.")
    return(invisible(FALSE))
  }
  scores <- utils::read.csv(scores_csv, stringsAsFactors = FALSE)
  corr   <- utils::read.csv(corr_csv,   stringsAsFactors = FALSE)
  if (nrow(scores) == 0L) {
    log_message("viz_cross_compartment_bridge: scores empty; nothing to plot.")
    return(invisible(FALSE))
  }

  pal <- if (exists("ETIOLOGY_GROUP_COLORS", inherits = TRUE))
            ETIOLOGY_GROUP_COLORS else c(NIU = "#E21F26", Viral = "#397FB9")

  for (w in unique(scores$weighting)) {
    df <- dplyr::filter(scores, .data$weighting == w)
    if (nrow(df) < 3L) next

    # ---- Pull per-stratum stats from the correlation table -----------------
    pooled <- corr |>
      dplyr::filter(.data$weighting == w, .data$stratum == "pooled") |>
      dplyr::slice_head(n = 1)
    partial <- corr |>
      dplyr::filter(.data$weighting == w,
                    .data$stratum == "partial_controlling_Phenotype_2") |>
      dplyr::slice_head(n = 1)
    by_eti <- corr |>
      dplyr::filter(.data$weighting == w,
                    grepl("^etiology:", .data$stratum)) |>
      dplyr::mutate(etiology = sub("^etiology:", "", .data$stratum))

    # ============================================================
    # Primary panel: within-etiology faceted scatter
    # ============================================================
    eti_subtitle_lines <- vapply(seq_len(nrow(by_eti)), function(i) {
      sprintf("%s  n=%d  r=%.2f  95%% CI [%.2f, %.2f]  perm p=%.3g",
              by_eti$etiology[i], by_eti$n[i], by_eti$pearson_r[i],
              by_eti$pearson_ci_lo[i], by_eti$pearson_ci_hi[i],
              by_eti$permutation_p[i])
    }, character(1))
    partial_line <- if (nrow(partial) > 0L &&
                        !is.na(partial$pearson_r))
      sprintf("Partial r (controlling Phenotype_2)  n=%d  r=%.2f  95%% CI [%.2f, %.2f]  perm p=%.3g",
              partial$n, partial$pearson_r, partial$pearson_ci_lo,
              partial$pearson_ci_hi, partial$permutation_p)
    else ""
    primary_subtitle <- paste(c(eti_subtitle_lines, partial_line),
                              collapse = "\n")

    p_within <- ggplot(df, aes(.data$myeloid_pc1, .data$tcell_pc1,
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
      labs(title = "Within-etiology coupling of myeloid + T cell PC1",
           subtitle = primary_subtitle,
           x = "Myeloid PC1 (NIU <-> Viral)",
           y = "T cell PC1 (NIU <-> Viral)") +
      theme_bw(base_size = 11) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(size = 9, lineheight = 1.2),
            strip.text    = element_text(face = "bold"),
            aspect.ratio  = 1)

    save_pdf_png(p_within,
                 file.path(cc$viz,
                           paste0("myeloid_tcell_pc1_within_etiology_", w)),
                 w = 11, h = 6.5)

    # ============================================================
    # Partial residual-residual scatter (controls for Phenotype_2)
    # ============================================================
    rx <- tryCatch(stats::resid(stats::lm(df$myeloid_pc1 ~ df$Phenotype_2)),
                   error = function(e) NULL)
    ry <- tryCatch(stats::resid(stats::lm(df$tcell_pc1   ~ df$Phenotype_2)),
                   error = function(e) NULL)
    if (!is.null(rx) && !is.null(ry)) {
      df_p <- df
      df_p$resid_myeloid <- rx
      df_p$resid_tcell   <- ry
      partial_subtitle <- if (nrow(partial) > 0L)
        sprintf("Residuals after lm(PC1 ~ Phenotype_2).  Partial r=%.2f  95%% CI [%.2f, %.2f]  perm p=%.3g",
                partial$pearson_r, partial$pearson_ci_lo,
                partial$pearson_ci_hi, partial$permutation_p)
      else "Partial correlation (stats unavailable)"
      p_partial <- ggplot(df_p, aes(.data$resid_myeloid, .data$resid_tcell,
                                    color = .data$Phenotype_2)) +
        geom_hline(yintercept = 0, linetype = "dashed",
                   linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.3, color = "grey70") +
        geom_smooth(method = "lm", se = TRUE, linewidth = 0.6,
                    formula = y ~ x,
                    aes(group = 1), color = "black", alpha = 0.15) +
        geom_point(size = 5.2, alpha = 0.9) +
        ggrepel::geom_text_repel(aes(label = .data$subject), size = 5.2,
                                  box.padding = 0.4, max.overlaps = Inf,
                                  show.legend = FALSE) +
        scale_color_manual(values = pal, name = NULL) +
        labs(title = "Partial coupling: residual myeloid vs T cell PC1",
             subtitle = partial_subtitle,
             x = "Myeloid PC1 residual (after Phenotype_2)",
             y = "T cell PC1 residual (after Phenotype_2)") +
        theme_bw(base_size = 11) +
        theme(plot.title    = element_text(face = "bold"),
              plot.subtitle = element_text(size = 9, lineheight = 1.2),
              aspect.ratio  = 1)
      save_pdf_png(p_partial,
                   file.path(cc$viz,
                             paste0("myeloid_tcell_pc1_partial_correlation_", w)),
                   w = 7.5, h = 7.5)
    }

    # ============================================================
    # Pooled scatter (kept as supplement - between-group separation
    # visible but flagged as confounded)
    # ============================================================
    pooled_subtitle <- if (nrow(pooled) > 0L)
      sprintf("Pooled (CONFOUNDED by NIU/Viral PC1 orientation) n=%d  r=%.2f  95%% CI [%.2f, %.2f]  perm p=%.3g",
              pooled$n, pooled$pearson_r, pooled$pearson_ci_lo,
              pooled$pearson_ci_hi, pooled$permutation_p)
    else "Pooled scatter (CONFOUNDED)"
    p_pooled <- ggplot(df, aes(.data$myeloid_pc1, .data$tcell_pc1,
                               color = .data$Phenotype_2)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 linewidth = 0.3, color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed",
                 linewidth = 0.3, color = "grey70") +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.6,
                  formula = y ~ x, aes(group = 1),
                  color = "grey40", linetype = "dashed") +
      geom_point(size = 2.6, alpha = 0.9) +
      scale_color_manual(values = pal, name = NULL) +
      labs(title = "Pooled myeloid vs T cell PC1 (supplement only)",
           subtitle = pooled_subtitle,
           x = "Myeloid PC1 (NIU <-> Viral)",
           y = "T cell PC1 (NIU <-> Viral)") +
      theme_bw(base_size = 11) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(size = 9, lineheight = 1.2),
            aspect.ratio  = 1)
    save_pdf_png(p_pooled,
                 file.path(cc$viz,
                           paste0("myeloid_tcell_pc1_pooled_correlation_", w)),
                 w = 7.5, h = 7.5)
  }

  # ============================================================
  # Per-(myeloid substate x tcell substate) heatmap. Uses partial r so
  # within-patient coupling is the visible signal. Cells with permutation
  # p < 0.05 are marked with an asterisk and a thick black border so the
  # significant pairs pop visually.
  # ============================================================
  per_sub_csv <- file.path(cc$tables, "pc1_bridge_per_substate.csv")
  if (file.exists(per_sub_csv)) {
    ps <- utils::read.csv(per_sub_csv, stringsAsFactors = FALSE)
    if (nrow(ps) > 0L) {
      ps$myeloid_display <- vapply(as.character(ps$myeloid_substate),
                                   function(id) get_substate_display(cfg, "myeloid", id),
                                   character(1))
      ps$tcell_display   <- vapply(as.character(ps$tcell_substate),
                                   function(id) get_substate_display(cfg, "tcell", id),
                                   character(1))
      ps$myeloid_display <- factor(ps$myeloid_display,
                                   levels = sort(unique(ps$myeloid_display)))
      ps$tcell_display   <- factor(ps$tcell_display,
                                   levels = sort(unique(ps$tcell_display)))
      ps$sig <- !is.na(ps$partial_permutation_p) & ps$partial_permutation_p < 0.05
      ps$label <- ifelse(
        is.na(ps$partial_pearson_r), "",
        sprintf("%s%.2f\n(n=%d)",
                ifelse(ps$sig, "* ", ""),
                ps$partial_pearson_r, ps$n))
      lim <- max(abs(ps$partial_pearson_r), na.rm = TRUE)
      lim <- if (is.finite(lim)) lim else 1
      ps_sig <- dplyr::filter(ps, .data$sig)
      p_heat <- ggplot(ps, aes(.data$tcell_display, .data$myeloid_display,
                               fill = .data$partial_pearson_r)) +
        geom_tile(color = "grey95", linewidth = 0.4) +
        # Overlay a thick black border on significant cells. Use geom_tile
        # again restricted to sig rows, with fill = NA so the original color
        # shows through.
        geom_tile(data = ps_sig, fill = NA, color = "black",
                  linewidth = 1.0) +
        geom_text(aes(label = .data$label,
                      fontface = ifelse(.data$sig, "bold", "plain")),
                  size = 2.8) +
        viridis::scale_fill_viridis(option = "viridis",
                                    limits = c(-lim, lim),
                                    name = "Partial Pearson r\n(within etiology)") +
        coord_fixed(ratio = 1) +
        labs(title = "Per-substate coupling (partial r controlling Phenotype_2)",
             subtitle = paste0("Subject pseudobulk PC1, myeloid sender (y) vs ",
                               "T cell receiver (x). * marks permutation p < 0.05. ",
                               sum(ps$sig, na.rm = TRUE), " of ",
                               nrow(ps), " cells significant."),
             x = "T cell substate", y = "Myeloid substate") +
        theme_minimal(base_size = 11) +
        theme(plot.title    = element_text(face = "bold"),
              plot.subtitle = element_text(size = 9),
              panel.grid    = element_blank(),
              axis.text.x   = element_text(angle = 30, hjust = 1))
      # Square cells via coord_fixed; canvas dimension follows the grid
      # aspect so the plot is visually balanced.
      n_x <- length(levels(ps$tcell_display))
      n_y <- length(levels(ps$myeloid_display))
      cell_in <- 0.85
      save_pdf_png(p_heat, file.path(cc$viz, "myeloid_tcell_pc1_per_substate_coupling_heatmap"),
                   w = n_x * cell_in + 4.5,
                   h = n_y * cell_in + 3.5)

    }
  }
  invisible(TRUE)
}

# Driver for panel G3 (ligand-family bar). Reads the per-condition
# family-summary tables. Writes:
#   liana_myeloid_tcell_ligand_family_by_arm        bar plot of LR-family signal per etiology,
#                                using config myeloid_programs categories
viz_cross_compartment_liana <- function(cfg) {
  cc <- .cc_paths(cfg)
  ensure_dir(cc$viz)
  lcfg <- cfg$liana %||% list()
  conditions <- as.character(lcfg$conditions %||% c("NIU", "Viral"))

  pal <- if (exists("ETIOLOGY_GROUP_COLORS", inherits = TRUE))
            ETIOLOGY_GROUP_COLORS else c(NIU = "#E21F26", Viral = "#397FB9")

  # --- G3: LR-family aggregated bar per arm ---
  fam_rows <- list()
  for (cond in conditions) {
    fam_csv <- file.path(cc$tables,
                         paste0("liana_family_summary_", cond, ".csv"))
    if (!file.exists(fam_csv)) next
    fam_rows[[cond]] <- utils::read.csv(fam_csv, stringsAsFactors = FALSE)
  }
  if (length(fam_rows) > 0L) {
    fams <- dplyr::bind_rows(fam_rows)
    if (nrow(fams) > 0L) {
      p_fam <- ggplot(fams,
                      aes(.data$ligand_family, .data$n_top5pct,
                          fill = .data$condition)) +
        geom_col(position = position_dodge(width = 0.75),
                 width = 0.7) +
        scale_fill_manual(values = pal, name = NULL) +
        labs(title = "Strong LR pairs per ligand family, by etiology",
             subtitle = paste0("LR pair counts with aggregate_rank <= 0.05 ",
                               "(per LIANA RRA). Families from cfg$myeloid_programs."),
             x = "Ligand family", y = "Number of LR pairs (rank <= 0.05)") +
        theme_bw(base_size = 11) +
        theme(plot.title    = element_text(face = "bold"),
              plot.subtitle = element_text(size = 9),
              axis.text.x   = element_text(angle = 30, hjust = 1))
      save_pdf_png(p_fam, file.path(cc$viz, "liana_myeloid_tcell_ligand_family_by_arm"),
                   w = 9, h = 6)
    }
  }
  invisible(TRUE)
}

# Classify an LR pair into a functional category. Operates on the parsed
# subunits (LIANA complexes use "_" / "&" / "+" as separators) so e.g.
# CD86_CD80 -> CTLA4 is recognized as Costim/Checkpoint by either subunit.
# Rule order is intentional — more specific receptor pairings are matched
# before broader ligand-family rules, so HLA-DR -> LAG3 sorts under
# Costim/Checkpoint rather than MHC-II antigen presentation.
.liana_lr_functional_group <- function(ligand, receptor) {
  stopifnot(length(ligand) == length(receptor))
  split_sub <- function(s) {
    if (is.na(s) || nchar(s) == 0L) return(character(0))
    unlist(strsplit(as.character(s), "[_&+]"))
  }
  classify_one <- function(L_parts, R_parts) {
    has_L <- function(re) any(grepl(re, L_parts))
    has_R <- function(re) any(grepl(re, R_parts))
    in_L  <- function(s) any(L_parts %in% s)
    in_R  <- function(s) any(R_parts %in% s)
    # Checkpoint / costimulation (specific receptor pairings first).
    if (has_L("^HLA-D[RPQ]") && in_R("LAG3"))            return("Costim/Checkpoint")
    if (in_L(c("CD80","CD86"))   && in_R(c("CTLA4","PDCD1")))
                                                          return("Costim/Checkpoint")
    if (in_L(c("CD274","PDCD1LG2","LGALS9","HHLA2","VTCN1")) &&
        in_R(c("PDCD1","HAVCR2","TIGIT","BTLA")))         return("Costim/Checkpoint")
    if (in_L(c("CD80","CD86","CD40","ICOSLG","TNFSF4","TNFSF9","TNFSF18")) &&
        in_R(c("CD28","ICOS","CD40LG","TNFRSF4","TNFRSF9","TNFRSF18")))
                                                          return("Costim/Checkpoint")
    # MHC class I antigen presentation.
    if (in_L(c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","B2M")))
                                                          return("MHC-I")
    # MHC class II antigen presentation.
    if (has_L("^HLA-D[RPQ]") || in_L("CD74"))             return("MHC-II")
    # Chemokine signaling.
    if (has_L("^(CCL|CXCL|XCL|CX3CL)") ||
        has_R("^(CCR|CXCR|XCR|CX3CR)"))                   return("Chemokine")
    # Alarmins / DAMPs.
    if (has_L("^S100A"))                                  return("Alarmin/DAMP")
    # Cytokine signaling.
    if (has_L("^IL[0-9]") || has_L("^IFN[A-Z]?[GBL]?[0-9]?$") ||
        has_L("^TGFB[0-9]?$") || has_L("^TNF$|^TNFSF[0-9]"))
                                                          return("Cytokine")
    if (has_R("^IL[0-9].*R[A-Z]?$") || has_R("^IFN.*R[12]?$") ||
        has_R("^TNFR"))                                   return("Cytokine")
    # Complement.
    if (has_L("^C1Q") || in_L(c("C3","C5")))              return("Complement")
    # Cytotoxic effectors.
    if (has_L("^GZM") || in_L(c("PRF1","GNLY")))          return("Cytotoxic effector")
    # Adhesion / integrins / ECM.
    if (has_R("^ITG[AB]") ||
        in_R(c("CD44","SELL","SELE","SELP","ICAM1","ICAM2","ICAM3","VCAM1","CD2")))
                                                          return("Adhesion/ECM")
    if (in_L(c("ICAM1","ICAM2","VCAM1","VCAN","FN1","TNC")) ||
        has_L("^COL[0-9]") || has_L("^F13A1$") ||
        has_L("^FBLN") || has_L("^LAM[ABC]") || has_L("^THBS"))
                                                          return("Adhesion/ECM")
    "Other"
  }
  fams <- vapply(seq_along(ligand), function(i) {
    classify_one(split_sub(ligand[i]), split_sub(receptor[i]))
  }, character(1))
  factor(fams,
         levels = c("MHC-I", "MHC-II", "Costim/Checkpoint",
                    "Chemokine", "Cytokine", "Adhesion/ECM",
                    "Alarmin/DAMP", "Complement",
                    "Cytotoxic effector", "Other"))
}

# ---------------------------------------------------------------------------
# Panel G synthesis — bridge-conditioned L-R interaction summary.
#
# The other G panels show the full LIANA landscape; Panel G in the figure is
# the take-home. This reads the per-(myeloid x tcell) substate bridge table,
# keeps only cluster pairs significantly coupled in Panel F
# (partial_permutation_p < sig_p), then for each one renders the top-N LR
# pairs (by best aggregate_rank across NIU/Viral) as a single dotplot:
#   x = significant cluster pair (with partial r annotated)
#   y = LR pair, grouped by functional family
#   dot size  = -log10(best aggregate_rank across conditions)
#   dot color = disease_bias_logfc (NIU red, shared grey, Viral blue —
#               matches ETIOLOGY_GROUP_COLORS used in Panel C / Panel F)
#
# Output: liana_myeloid_tcell_summary_lr.{pdf,png}
# Audit:  liana_myeloid_tcell_summary_lr_pairs.csv
.liana_viz_panel_g_summary <- function(cc, cfg,
                                       sig_p = 0.05,
                                       candidate_pool_per_pair = 15L,
                                       n_lr_pairs = 25L,
                                       split_shared_private = FALSE,
                                       facet_by = c("tcell_target",
                                                    "wiring_group"),
                                       tcell_facet_order = c("4","3","1","5"),
                                       w = 5, h = 1.25,
                                       lr_pair_whitelist = NULL,
                                       collapse_singletons = TRUE,
                                       output_suffix = "",
                                       hide_y_axis = FALSE) {
  facet_by <- match.arg(facet_by)
  compact  <- isTRUE(h <= 2 || w <= 6)
  bridge_csv   <- file.path(cc$tables, "pc1_bridge_per_substate.csv")
  combined_csv <- file.path(cc$tables, "liana_myeloid_to_tcell_combined.csv")
  if (!file.exists(bridge_csv) || !file.exists(combined_csv)) {
    log_message("  Panel G summary: missing bridge or LIANA combined CSV.")
    return(invisible(FALSE))
  }
  bridge   <- utils::read.csv(bridge_csv,   stringsAsFactors = FALSE)
  combined <- utils::read.csv(combined_csv, stringsAsFactors = FALSE)
  if (!nrow(bridge) || !nrow(combined)) {
    log_message("  Panel G summary: empty input table; skipping.")
    return(invisible(FALSE))
  }

  sig <- bridge[!is.na(bridge$partial_permutation_p) &
                bridge$partial_permutation_p < sig_p, , drop = FALSE]
  if (!nrow(sig)) {
    log_message("  Panel G summary: no cluster pairs at p < ", sig_p,
                "; nothing to plot.")
    return(invisible(FALSE))
  }

  sig$source  <- paste0("myeloid_", sig$myeloid_substate)
  sig$target  <- paste0("tcell_",   sig$tcell_substate)
  sig$src_lbl <- vapply(as.character(sig$myeloid_substate),
                        function(id) get_substate_display(cfg, "myeloid", id),
                        character(1))
  sig$tgt_lbl <- vapply(as.character(sig$tcell_substate),
                        function(id) get_substate_display(cfg, "tcell", id),
                        character(1))
  # Per-pair stats (r, n) live in Panel F; the x-axis label here is just the
  # sender -> receiver identity so the columns read cleanly. When faceting by
  # T cell target we also expose the sender alone for the column tick label.
  sig$pair_lbl   <- sprintf("%s ->\n%s", sig$src_lbl, sig$tgt_lbl)
  sig$myel_lbl   <- sig$src_lbl
  sig$tcell_lbl  <- sig$tgt_lbl

  combined$lr_pair <- paste(combined$ligand_complex,
                            combined$receptor_complex, sep = " -> ")
  rk_a <- combined$aggregate_rank_NIU
  rk_b <- combined$aggregate_rank_Viral
  rk_a[is.na(rk_a)] <- Inf
  rk_b[is.na(rk_b)] <- Inf
  combined$best_rank <- pmin(rk_a, rk_b)

  per_pair <- merge(sig[, c("source", "target", "pair_lbl",
                            "myel_lbl", "tcell_lbl", "myeloid_substate",
                            "tcell_substate", "partial_pearson_r")],
                    combined, by = c("source", "target"))
  if (!nrow(per_pair)) {
    log_message("  Panel G summary: no LIANA rows match significant cluster ",
                "pairs; check that LIANA was run on the same substates.")
    return(invisible(FALSE))
  }

  # Candidate pool: top-K LR per cluster pair. Used both to compute the
  # sharing score and as the universe from which the final LR set is drawn.
  per_pair <- per_pair[order(per_pair$pair_lbl, per_pair$best_rank),
                       , drop = FALSE]
  rank_within_pair <- stats::ave(seq_len(nrow(per_pair)), per_pair$pair_lbl,
                                 FUN = function(i) seq_along(i))
  pool <- per_pair[rank_within_pair <= candidate_pool_per_pair &
                   is.finite(per_pair$best_rank), , drop = FALSE]
  if (!nrow(pool)) {
    log_message("  Panel G summary: candidate pool empty.")
    return(invisible(FALSE))
  }

  # Sharing score = number of distinct cluster pairs where this LR is in the
  # candidate pool. Tie-break by best (smallest) aggregate_rank seen for the LR.
  share_tbl <- stats::aggregate(
    list(n_pairs = pool$pair_lbl, best_rank = pool$best_rank),
    by = list(lr_pair = pool$lr_pair),
    FUN = function(x) x)
  # aggregate() with FUN=function(x) x returns a list-column. Compute n / min
  # explicitly to be safe across R versions.
  uniq_lr <- unique(pool$lr_pair)
  share_df <- data.frame(
    lr_pair  = uniq_lr,
    n_pairs  = vapply(uniq_lr, function(g)
                      length(unique(pool$pair_lbl[pool$lr_pair == g])),
                      integer(1)),
    min_rank = vapply(uniq_lr, function(g)
                      min(pool$best_rank[pool$lr_pair == g], na.rm = TRUE),
                      numeric(1)),
    stringsAsFactors = FALSE)
  share_df <- share_df[order(-share_df$n_pairs, share_df$min_rank), ]
  # Curated whitelist overrides the sharing-based selection. Useful for the
  # main-figure panel where pair choice is driven by which axes interact
  # with the manuscript text. Order in the whitelist is preserved as the
  # y-axis order (after row-grouping by family below).
  if (!is.null(lr_pair_whitelist) && length(lr_pair_whitelist)) {
    keep_lr <- as.character(lr_pair_whitelist)
    missing_lr <- setdiff(keep_lr, pool$lr_pair)
    if (length(missing_lr)) {
      log_message("  Panel G summary: ", length(missing_lr),
                  " whitelist LR pair(s) absent from candidate pool: ",
                  paste(missing_lr, collapse = ", "))
      keep_lr <- intersect(keep_lr, pool$lr_pair)
    }
  } else {
    keep_lr <- utils::head(share_df$lr_pair, n_lr_pairs)
  }
  panel_df <- pool[pool$lr_pair %in% keep_lr, , drop = FALSE]
  if (!nrow(panel_df)) {
    log_message("  Panel G summary: no LR rows survived sharing filter.")
    return(invisible(FALSE))
  }

  panel_df$family <- as.character(
    .liana_lr_functional_group(panel_df$ligand_complex,
                               panel_df$receptor_complex))
  panel_df$n_pairs <- share_df$n_pairs[match(panel_df$lr_pair,
                                              share_df$lr_pair)]

  # Collapse singleton families (only one displayed LR pair) into "Other" so
  # the row strips don't waste vertical space on a single-row category.
  # Skipped when a curated whitelist is in use — there the family signal
  # itself is the message and singletons (e.g. one Costim row, one MHC-II
  # row) should keep their own strip.
  if (isTRUE(collapse_singletons) && is.null(lr_pair_whitelist)) {
    fam_per_lr <- unique(panel_df[, c("lr_pair", "family")])
    fam_n      <- table(fam_per_lr$family)
    singletons <- setdiff(names(fam_n)[fam_n <= 1L], "Other")
    if (length(singletons))
      panel_df$family[panel_df$family %in% singletons] <- "Other"
  }

  fam_levels_full <- c("MHC-I", "MHC-II", "Costim/Checkpoint",
                       "Chemokine", "Cytokine", "Adhesion/ECM",
                       "Alarmin/DAMP", "Complement",
                       "Cytotoxic effector", "Other")
  surviving <- intersect(fam_levels_full, unique(panel_df$family))
  panel_df$family <- factor(panel_df$family, levels = surviving)

  # X-axis facet variable. Default: T cell target so each facet collects all
  # myeloid senders feeding into one T cell cluster (the biological grouping
  # the user landed on after the wiring-group prototype). Wiring-group (cut
  # of an hclust on LR-presence) is kept as an alternative behind the
  # `facet_by` argument for the supplement.
  if (facet_by == "tcell_target") {
    panel_df$tcell_facet <- sprintf("%s: %s",
                                    panel_df$tcell_substate,
                                    panel_df$tcell_lbl)
    tcell_keys <- unique(data.frame(
      id    = as.character(panel_df$tcell_substate),
      label = panel_df$tcell_facet,
      stringsAsFactors = FALSE))
    ord_ids <- c(intersect(tcell_facet_order, tcell_keys$id),
                 setdiff(tcell_keys$id, tcell_facet_order))
    panel_df$tcell_facet <- factor(panel_df$tcell_facet,
                                    levels = tcell_keys$label[match(ord_ids,
                                                                    tcell_keys$id)])
    # X tick label = myeloid sender only (T cell identity is in the strip).
    myel_keys <- unique(data.frame(
      id    = as.integer(panel_df$myeloid_substate),
      label = as.character(panel_df$myel_lbl),
      stringsAsFactors = FALSE))
    myel_keys <- myel_keys[order(myel_keys$id), ]
    panel_df$x_label <- factor(as.character(panel_df$myel_lbl),
                                levels = myel_keys$label)
  } else {
    # Wiring-group fallback: hclust on binary LR-presence, cut at k=3.
    pair_levels <- unique(panel_df$pair_lbl)
    pair_group  <- NULL
    if (length(pair_levels) >= 3L) {
      presence <- table(factor(panel_df$lr_pair, levels = keep_lr),
                         factor(panel_df$pair_lbl, levels = pair_levels))
      presence[presence > 0] <- 1L
      presence_mat <- t(as.matrix(presence))
      hc <- tryCatch(stats::hclust(stats::dist(presence_mat, method = "binary"),
                                    method = "average"),
                     error = function(e) NULL)
      if (!is.null(hc)) {
        pair_order <- rownames(presence_mat)[hc$order]
        k <- min(3L, length(pair_levels) - 1L); k <- max(2L, k)
        clu <- stats::cutree(hc, k = k)
        group_letters <- character(0)
        for (g in clu[pair_order]) {
          gk <- as.character(g)
          if (!gk %in% names(group_letters))
            group_letters[gk] <- LETTERS[length(group_letters) + 1L]
        }
        pair_group <- setNames(
          paste0("Wiring group ", unname(group_letters[as.character(clu)])),
          names(clu))
      } else {
        pair_order <- pair_levels
      }
    } else {
      pair_order <- pair_levels
    }
    panel_df$x_label <- factor(panel_df$pair_lbl, levels = pair_order)
    if (!is.null(pair_group)) {
      panel_df$tcell_facet <- factor(
        pair_group[as.character(panel_df$pair_lbl)],
        levels = unique(pair_group[pair_order]))
    }
  }

  # Row ordering: functional family (existing facet), then either the
  # whitelist input order (curated panel) or sharing count DESC + alphabetic.
  lr_meta <- unique(panel_df[, c("lr_pair", "family", "n_pairs")])
  if (!is.null(lr_pair_whitelist) && length(lr_pair_whitelist)) {
    lr_meta$wl_idx <- match(lr_meta$lr_pair, as.character(lr_pair_whitelist))
    lr_meta <- lr_meta[order(lr_meta$family, lr_meta$wl_idx), ]
  } else {
    lr_meta <- lr_meta[order(lr_meta$family,
                              -lr_meta$n_pairs,
                              lr_meta$lr_pair), ]
  }
  panel_df$lr_pair <- factor(panel_df$lr_pair, levels = rev(lr_meta$lr_pair))

  # Optional two-tier split: shared (n_pairs >= 2) vs private (n_pairs == 1).
  if (isTRUE(split_shared_private)) {
    panel_df$tier <- factor(
      ifelse(panel_df$n_pairs >= 2L,
             "Shared across pairs (n>=2)",
             "Private to one pair"),
      levels = c("Shared across pairs (n>=2)", "Private to one pair"))
  }

  bias_lim <- max(abs(panel_df$disease_bias_logfc), na.rm = TRUE)
  if (!is.finite(bias_lim) || bias_lim == 0) bias_lim <- 1

  pal <- if (exists("ETIOLOGY_GROUP_COLORS", inherits = TRUE))
            ETIOLOGY_GROUP_COLORS else c(NIU = "#E21F26", Viral = "#397FB9")

  n_cols <- length(unique(panel_df$pair_lbl))
  n_rows <- length(unique(panel_df$lr_pair))
  n_shared <- sum(lr_meta$n_pairs >= 2L)

  has_facet <- "tcell_facet" %in% colnames(panel_df)
  facet_spec <- if (isTRUE(hide_y_axis)) {
    # Drop row faceting; y-axis ordering still respects family but no row
    # strips render. y-axis tick labels are preserved so readers can still
    # identify each LR pair.
    if (has_facet)
      facet_grid(cols = vars(.data$tcell_facet),
                 scales = "free_x", space = "free_x")
    else
      NULL
  } else if (isTRUE(split_shared_private)) {
    if (has_facet)
      facet_grid(rows = vars(.data$tier, .data$family),
                 cols = vars(.data$tcell_facet),
                 scales = "free", space = "free", switch = "y")
    else
      facet_grid(rows = vars(.data$tier, .data$family),
                 scales = "free_y", space = "free_y", switch = "y")
  } else if (has_facet) {
    facet_grid(rows = vars(.data$family),
               cols = vars(.data$tcell_facet),
               scales = "free", space = "free", switch = "y")
  } else {
    facet_grid(rows = vars(.data$family),
               scales = "free_y", space = "free_y", switch = "y")
  }

  p <- ggplot(panel_df,
              aes(.data$x_label, .data$lr_pair,
                  size  = -log10(.data$best_rank + 1e-6),
                  color = .data$disease_bias_logfc)) +
    geom_point()
  if (!is.null(facet_spec)) p <- p + facet_spec
  p <- p +
    scale_color_gradient2(low = unname(pal[["Viral"]]),
                          mid = "grey92",
                          high = unname(pal[["NIU"]]),
                          midpoint = 0,
                          limits = c(-bias_lim, bias_lim),
                          oob = scales::squish,
                          name = "NIU - Viral\nlogFC bias") +
    scale_size_continuous(name = "-log10(best rank)",
                          range = c(0.6, if (compact) 3.5 else 6)) +
    labs(x = NULL, y = NULL)

  if (!compact) {
    p <- p +
      labs(title = "Ligand-receptor wiring of significantly coupled myeloid -> T cell pairs",
           subtitle = paste0("Columns: T cell receivers (significant pairs from Panel F). ",
                             "Within each facet: myeloid senders. ",
                             "Rows: top ", n_lr_pairs,
                             " LR interactions prioritized by sharing across ",
                             "cluster pairs (", n_shared, "/", n_rows,
                             " in 2+ pairs). Color: shared (grey) vs ",
                             "disease-private (NIU red / Viral blue)."))
  }

  base_sz <- if (compact) 6 else 10
  p <- p +
    theme_bw(base_size = base_sz) +
    theme(plot.title         = element_text(face = "bold"),
          plot.subtitle      = element_text(size = 9, lineheight = 1.2),
          strip.text.y.left  = element_text(angle = 0, face = "bold",
                                            size = if (compact) 5 else 9,
                                            hjust = 1),
          strip.text.x       = element_text(face = "bold",
                                            size = if (compact) 5 else 9),
          strip.placement    = "outside",
          strip.background.y = element_rect(fill = "grey92", color = NA),
          strip.background.x = element_rect(fill = "grey85", color = NA),
          axis.text.x        = element_blank(),
          axis.ticks.x       = element_blank(),
          axis.text.y        = element_text(size = if (compact) 4.5 else 8),
          panel.grid.minor   = element_blank(),
          panel.spacing.x    = grid::unit(if (compact) 0.25 else 0.8, "lines"),
          legend.position    = if (compact) "none" else "right",
          plot.margin        = grid::unit(rep(if (compact) 1 else 5, 4), "pt"))

  base_name <- paste0("liana_myeloid_tcell_summary_lr", output_suffix)
  save_pdf_png(p, file.path(cc$viz, base_name), w = w, h = h)

  audit <- panel_df[, c("pair_lbl", "partial_pearson_r", "lr_pair",
                        "family", "n_pairs",
                        "ligand_complex", "receptor_complex",
                        "aggregate_rank_NIU", "aggregate_rank_Viral",
                        "best_rank",
                        "logfc.logfc_comb_NIU", "logfc.logfc_comb_Viral",
                        "disease_bias_logfc", "unique_to")]
  utils::write.csv(audit,
                   file.path(cc$tables,
                             paste0(base_name, "_pairs.csv")),
                   row.names = FALSE)
  log_message("  Panel G summary: wrote ", nrow(audit),
              " LR rows across ", n_cols, " cluster pairs, ",
              n_shared, "/", n_rows, " LR pairs shared (n>=2).")
  invisible(TRUE)
}

# Panel H driver. Heatmap of NicheNet ligand activity x myeloid sender
# substate, faceted by T cell substate x pole. No-op when the NicheNet CSVs
# aren't on disk.
viz_cross_compartment_nichenet <- function(cfg) {
  cc <- .cc_paths(cfg)
  ensure_dir(cc$viz)
  niu_csv   <- file.path(cc$tables, "nichenet_NIU_pole_ligands.csv")
  viral_csv <- file.path(cc$tables, "nichenet_Viral_pole_ligands.csv")
  if (!file.exists(niu_csv) || !file.exists(viral_csv)) {
    log_message("viz_cross_compartment_nichenet: NicheNet CSVs missing; ",
                "skipping panel H.")
    return(invisible(FALSE))
  }
  niu   <- utils::read.csv(niu_csv,   stringsAsFactors = FALSE)
  viral <- utils::read.csv(viral_csv, stringsAsFactors = FALSE)
  df <- dplyr::bind_rows(niu, viral)
  if (nrow(df) == 0L) {
    log_message("viz_cross_compartment_nichenet: nothing to plot.")
    return(invisible(FALSE))
  }

  p <- ggplot(df, aes(.data$sending_myeloid_substate_top, .data$ligand,
                      fill = .data$pearson)) +
    geom_tile(color = "grey90", linewidth = 0.2) +
    viridis::scale_fill_viridis(option = "viridis",
                                name = "Pearson activity") +
    facet_grid(.data$tcell_substate ~ .data$pole, scales = "free_y",
               space = "free_y") +
    labs(title = "NicheNet ligand activity per T cell substate x disease pole",
         x = "Top myeloid sender substate",
         y = "Predicted ligand") +
    theme_bw(base_size = 10) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          strip.text.y = element_text(angle = 0))
  save_pdf_png(p, file.path(cc$viz, "nichenet_myeloid_tcell_ligand_activity"),
               w = max(7, length(unique(df$sending_myeloid_substate_top)) * 0.6 + 4),
               h = max(6, length(unique(df$ligand)) * 0.18 + 3))
  invisible(TRUE)
}

# Orchestrator dispatched by R/82_viz_dispatch.R when target = "cross_compartment".
run_visualizations_cross_compartment <- function(cfg) {
  log_message("=== Cross-compartment (Figure 4 F-H) visualizations ===")
  viz_cross_compartment_bridge(cfg)
  viz_cross_compartment_liana(cfg)
  # Panel G take-home (main figure): curated 10-pair condensed dotplot, sized
  # for the figure slot (5 x 1.25 in). Pairs are hand-picked to interact with
  # the manuscript text — HLA-II/costim (CD86->CTLA4, HLA-DQA1->CD4), the
  # CXCR3-axis chemokine recruitment (CXCL10/11, CCL7), Treg activation
  # (ICAM1->IL2RA), the alarmin and adhesion hubs (S100A8->CD69, VCAN->CD44,
  # CD14->ITGA4), and the DC priming axis MRC1->PTPRC. Whitelist also disables
  # singleton-family collapse so each axis keeps its own row strip.
  panel_g_main_whitelist <- c(
    "CD86 -> CTLA4",
    "HLA-DQA1 -> CD4",
    "CXCL10 -> CXCR3",
    "CXCL11 -> CXCR3",
    "CCL7 -> CXCR3",
    "ICAM1 -> IL2RA",
    "S100A8 -> CD69",
    "VCAN -> CD44",
    "CD14 -> ITGA4",
    "MRC1 -> PTPRC")
  .liana_viz_panel_g_summary(.cc_paths(cfg), cfg,
                             lr_pair_whitelist = panel_g_main_whitelist,
                             w = 5, h = 1.25,
                             hide_y_axis = TRUE)
  # Expanded supplemental: top-40 LR pairs by sharing across cluster pairs,
  # singleton families collapsed for legibility, sized to render at full size.
  .liana_viz_panel_g_summary(.cc_paths(cfg), cfg,
                             n_lr_pairs = 40L,
                             candidate_pool_per_pair = 20L,
                             w = 12, h = 12,
                             output_suffix = "_expanded")
  viz_cross_compartment_nichenet(cfg)
  log_message("=== Cross-compartment visualizations complete ===")
  invisible(TRUE)
}
