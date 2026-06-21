# R/70_tcr_genex_signatures.R
# Link three orthogonal TCR-derived clone groups — expanded, public, and
# GLIPH-enriched — to T cell gene-expression programs, and then bridge the
# enriched T cell substates to their candidate myeloid partners using the
# pc1_bridge_per_substate.csv coupling table written by R/46.
#
# Clone groups (boolean per cell, computed on the eye T cell object):
#   expanded        : cloneSize %in% c("Large", "Hyperexpanded")
#   public          : CTstrict is observed in >= 2 distinct Subjects
#   gliph_enriched  : CTaa-derived CDR3b is a member of a GLIPH cluster with
#                     gliph_enrichment_viral_vs_niu$FDR < cfg$tcr_advanced$gliph_enrichment_fdr
#
# Outputs:
#   tcr_signature_scores_by_clone_group.csv  — UCell scores per (Subject, substate, group, panel)
#   tcr_pseudobulk_DE_{expanded,public,gliph}.csv
#   tcr_myeloid_coupling.csv

`%||%` <- function(x, y) if (is.null(x)) y else x

.tgs_label_clones <- function(obj, gliph_enriched_cdr3s) {
  meta <- obj@meta.data
  # expanded: scRepertoire stores cloneSize as a factor whose levels carry
  # the proportion range suffix (e.g. "Large (0.01 < X <= 0.1)"), so an
  # exact %in% c("Large","Hyperexpanded") never matches. Use prefix match.
  cs <- as.character(meta$cloneSize)
  expanded <- grepl("^Large|^Hyperexpanded", cs)
  expanded[is.na(expanded)] <- FALSE

  # public: count distinct Subjects per CTstrict
  ct <- as.character(meta$CTstrict)
  sb <- as.character(meta$Subject)
  ct_subj <- tapply(sb, ct, function(x) length(unique(x[!is.na(x)])))
  public  <- !is.na(ct) & ct_subj[ct] >= 2L
  public[is.na(public)] <- FALSE

  # gliph_enriched: CTaa -> TRB CDR3 -> in any significant GLIPH cluster
  trb <- stringr::str_split(as.character(meta$CTaa), "_", simplify = TRUE)
  trb_aa <- if (ncol(trb) >= 2) trb[, 2] else rep(NA_character_, nrow(meta))
  gliph_e <- !is.na(trb_aa) & trb_aa %in% gliph_enriched_cdr3s
  gliph_e[is.na(gliph_e)] <- FALSE

  obj$grp_expanded <- expanded
  obj$grp_public   <- public
  obj$grp_gliph    <- gliph_e
  obj
}

.tgs_pseudobulk_DE <- function(obj, group_col, cfg, target = "tcell") {
  # group_col is a boolean column already stamped on obj. Run within-substate
  # contrasts: TRUE vs FALSE per substate_key, returning the same long-form
  # data.frame run_pseudobulk_deseq2 produces.
  obj@meta.data[[group_col]] <- ifelse(obj@meta.data[[group_col]],
                                       "in_group", "out_group")
  cluster_col <- if ("substate_key" %in% colnames(obj@meta.data))
                   "substate_key" else "knn.leiden.cluster"
  tryCatch(
    run_pseudobulk_deseq2(obj,
                          group_col   = group_col,
                          group1      = "in_group",
                          group2      = "out_group",
                          cluster_col = cluster_col,
                          cfg         = cfg,
                          target      = target,
                          force_simple = TRUE),
    error = function(e) {
      log_message("  pseudobulk DE failed (", group_col, "): ",
                  conditionMessage(e))
      data.frame()
    })
}

.tgs_ucell <- function(obj, panels) {
  if (!requireNamespace("UCell", quietly = TRUE)) {
    log_message("  UCell not installed; skipping signature scoring.")
    return(NULL)
  }
  obj <- UCell::AddModuleScore_UCell(obj, features = panels, name = "_UCell")
  obj
}

.tgs_aggregate_ucell <- function(obj, panels) {
  if (is.null(obj)) return(data.frame())
  meta <- obj@meta.data
  panel_cols <- paste0(names(panels), "_UCell")
  panel_cols <- intersect(panel_cols, colnames(meta))
  if (length(panel_cols) == 0) return(data.frame())

  # Drop any pre-existing barcode column so rownames_to_column doesn't
  # collide. Upstream airr ingest stashes a `barcode` field on meta.data
  # that's redundant with the rownames (which are the canonical Seurat
  # cell key) — recreate it from rownames for consistency.
  if ("barcode" %in% colnames(meta)) meta$barcode <- NULL
  long <- meta |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(barcode, Subject, substate_key,
                  grp_expanded, grp_public, grp_gliph,
                  dplyr::all_of(panel_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(panel_cols),
                        names_to = "panel", values_to = "ucell")
  long$panel <- sub("_UCell$", "", long$panel)

  group_rows <- list()
  for (gcol in c("grp_expanded", "grp_public", "grp_gliph")) {
    agg <- long |>
      dplyr::group_by(Subject, substate_key, panel,
                      group = .data[[gcol]]) |>
      dplyr::summarise(mean_ucell = mean(ucell, na.rm = TRUE),
                       n = dplyr::n(),
                       .groups = "drop") |>
      dplyr::mutate(clone_group = sub("^grp_", "", gcol))
    group_rows[[gcol]] <- agg
  }
  dplyr::bind_rows(group_rows)
}

.tgs_wilcoxon_signatures <- function(agg) {
  if (nrow(agg) == 0) return(data.frame())
  # Drop the per-group cell count `n` before pivoting; with names_from=group
  # and values_from=mean_ucell, pivot_wider would otherwise treat `n` as an
  # id column. Since `n` differs between group=TRUE and group=FALSE rows for
  # the same (Subject, substate, panel, clone_group), the pivot would put
  # g_TRUE and g_FALSE on separate rows and the subsequent paired filter
  # would drop everything.
  rows <- agg |>
    dplyr::select(-dplyr::any_of("n")) |>
    tidyr::pivot_wider(id_cols = c(Subject, substate_key, panel, clone_group),
                       names_from = group, values_from = mean_ucell,
                       names_prefix = "g_") |>
    dplyr::filter(!is.na(g_TRUE), !is.na(g_FALSE)) |>
    dplyr::group_by(clone_group, substate_key, panel) |>
    dplyr::summarise(
      n_subjects = dplyr::n(),
      median_in  = stats::median(g_TRUE,  na.rm = TRUE),
      median_out = stats::median(g_FALSE, na.rm = TRUE),
      delta      = median_in - median_out,
      p          = tryCatch(stats::wilcox.test(g_TRUE, g_FALSE,
                                               paired = TRUE)$p.value,
                            error = function(e) NA_real_),
      .groups    = "drop") |>
    dplyr::group_by(clone_group) |>
    dplyr::mutate(fdr = stats::p.adjust(p, method = "BH")) |>
    dplyr::ungroup()
  rows
}

.tgs_myeloid_coupling <- function(per_substate_csv, enriched_tcell_substates,
                                  p_cut = 0.05) {
  if (!file.exists(per_substate_csv) || length(enriched_tcell_substates) == 0)
    return(data.frame())
  d <- utils::read.csv(per_substate_csv, stringsAsFactors = FALSE)
  # Expected columns from R/46_cross_compartment_bridge.R::.per_substate_bridge:
  #   myeloid_substate, tcell_substate, n, partial_r, partial_permutation_p, ...
  if (!all(c("myeloid_substate", "tcell_substate",
             "partial_r", "partial_permutation_p") %in% colnames(d)))
    return(data.frame())
  d <- d[d$tcell_substate %in% enriched_tcell_substates, , drop = FALSE]
  d <- d[!is.na(d$partial_permutation_p) &
         d$partial_permutation_p < p_cut &
         d$partial_r > 0, , drop = FALSE]
  if (nrow(d) == 0) return(d)
  d$coupling_strength <- -log10(pmax(d$partial_permutation_p, 1e-300)) *
                         sign(d$partial_r)
  d[order(-d$coupling_strength), ]
}

run_tcr_genex_signatures <- function(cfg) {
  if (!isTRUE(cfg$steps$tcr_genex_signatures)) {
    log_message("TCR-genex signatures disabled. Skipping.")
    return(invisible(TRUE))
  }
  log_message("Starting TCR-genex signature analysis...")

  paths_tcell <- get_target_paths(cfg, "tcell")
  out_tables  <- paths_tcell$results_tables
  ensure_dir(out_tables)

  tcell_rds <- file.path(paths_tcell$results_objects,
                         "IntegratedSeuratObject.rds")
  if (!file.exists(tcell_rds)) {
    log_message("  T cell compartment object not found at ", tcell_rds)
    return(invisible(FALSE))
  }
  obj <- readRDS(tcell_rds)
  # Older compartment objects pre-date the apply_substate_keys() saver, so
  # stamp substate_key from knn.leiden.cluster (the canonical compartment
  # cluster column) when absent. get_substate_key_vector lives in
  # R/23_substate_labels.R which run_pipeline.R sources before us.
  if (!"substate_key" %in% colnames(obj@meta.data)) {
    obj$substate_key <- get_substate_key_vector(obj, "tcell")
  }

  gliph_rds <- file.path(cfg$paths$results_objects, "ImmGLIPHResults.rds")
  gliph_enriched_cdr3s <- character(0)
  fdr_cut <- cfg$tcr_advanced$gliph_enrichment_fdr %||% 0.1
  if (file.exists(gliph_rds)) {
    g <- readRDS(gliph_rds)
    if (!is.null(g$enrich) && !is.null(g$clusters)) {
      sig_ids <- g$enrich$cluster_id[!is.na(g$enrich$FDR) &
                                     g$enrich$FDR < fdr_cut]
      gliph_enriched_cdr3s <- unique(
        g$clusters$CDR3b[g$clusters$cluster_id %in% sig_ids])
    }
  }
  log_message("  GLIPH-enriched CDR3b count (FDR<", fdr_cut, "): ",
              length(gliph_enriched_cdr3s))

  obj <- .tgs_label_clones(obj, gliph_enriched_cdr3s)
  log_message("  Cells: expanded=", sum(obj$grp_expanded),
              ", public=",          sum(obj$grp_public),
              ", gliph=",            sum(obj$grp_gliph))

  panels <- cfg$tcr_advanced$signature_panels
  if (is.null(panels) || length(panels) == 0) {
    log_message("  cfg$tcr_advanced$signature_panels missing; skipping UCell.")
  } else {
    obj <- .tgs_ucell(obj, panels)
    agg <- .tgs_aggregate_ucell(obj, panels)
    sig_tests <- .tgs_wilcoxon_signatures(agg)
    if (nrow(sig_tests) > 0) {
      utils::write.csv(sig_tests,
                       file.path(out_tables,
                                 "tcr_signature_scores_by_clone_group.csv"),
                       row.names = FALSE)
      log_message("  Saved: tcr_signature_scores_by_clone_group.csv (",
                  nrow(sig_tests), " rows)")
    }
  }

  # --- Pseudobulk DE per clone group, per substate ---------------------------
  for (gcol in c("grp_expanded", "grp_public", "grp_gliph")) {
    if (sum(obj@meta.data[[gcol]], na.rm = TRUE) < cfg$tcr_advanced$pseudobulk_min_cells %||% 20) {
      log_message("  Skipping DE for ", gcol, " (n cells too small).")
      next
    }
    log_message("  Pseudobulk DE: ", gcol, " (in_group vs out_group)...")
    de <- .tgs_pseudobulk_DE(obj, group_col = gcol, cfg = cfg, target = "tcell")
    if (nrow(de) > 0) {
      group_short <- sub("^grp_", "", gcol)
      out_path <- file.path(out_tables,
                            paste0("tcr_pseudobulk_DE_", group_short, ".csv"))
      utils::write.csv(de, out_path, row.names = FALSE)
      log_message("  Saved: ", basename(out_path), " (", nrow(de), " gene rows)")
    }
  }

  # --- Myeloid coupling --------------------------------------------------------
  cc_tables <- cfg$paths_cross_compartment$tables %||%
               "outputs/tables/cross_compartment/"
  per_sub_csv <- file.path(cc_tables, "pc1_bridge_per_substate.csv")
  enriched_tcell_substates <- character(0)
  meta <- obj@meta.data
  for (gcol in c("grp_expanded", "grp_public", "grp_gliph")) {
    if (!gcol %in% colnames(meta)) next
    by_sub <- tapply(meta[[gcol]], meta$substate_key, mean, na.rm = TRUE)
    # Enriched substates: any substate where the within-group fraction is
    # >= 2x the global fraction. A simple, interpretable rule of thumb that
    # works for all three clone group definitions.
    global <- mean(meta[[gcol]], na.rm = TRUE)
    if (is.finite(global) && global > 0) {
      enriched_tcell_substates <- union(
        enriched_tcell_substates,
        names(by_sub)[by_sub >= 2 * global])
    }
  }
  enriched_tcell_substates <- unique(enriched_tcell_substates)
  log_message("  Enriched T cell substates: ",
              paste(enriched_tcell_substates, collapse = ", "))

  coupling <- .tgs_myeloid_coupling(per_sub_csv, enriched_tcell_substates,
                                    p_cut = 0.05)
  if (nrow(coupling) > 0) {
    utils::write.csv(coupling,
                     file.path(out_tables, "tcr_myeloid_coupling.csv"),
                     row.names = FALSE)
    log_message("  Saved: tcr_myeloid_coupling.csv (", nrow(coupling), " edges)")
  } else {
    log_message("  No myeloid coupling edges met p<0.05; skipping output.")
  }

  log_message("TCR-genex signature analysis complete.")
  invisible(TRUE)
}
