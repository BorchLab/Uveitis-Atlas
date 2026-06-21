# R/01_setup_utils.R
suppressPackageStartupMessages({
  library(Matrix)
  library(glue)
})

log_message <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(glue("[{ts}] {paste(..., collapse=' ')}"))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}



# --- Dual-mode plot saving: labeled + stripped versions ----------------------
# Saves a ggplot as both a labeled version (standard) and a stripped version
# (black background, no axes/labels/legend). Maintains 1:1 aspect ratio.
dual_save_plot <- function(expr, base_path, width = 8, height = 8,
                           stripped_bg = "black", device = "pdf",
                           envir = parent.frame()) {
  tryCatch({
    p <- eval(expr, envir = envir)
    if (is.null(p)) return(invisible(NULL))

    ensure_dir(dirname(base_path))
    labeled_path  <- paste0(base_path, "_labeled.",  device)
    stripped_path <- paste0(base_path, "_stripped.", device)

    # Labeled version: standard with coord_fixed
    p_labeled <- p + ggplot2::coord_fixed(ratio = 1)
    ggsave(labeled_path, p_labeled, width = width, height = height)
    log_message("  Saved: ", basename(labeled_path))

    # Stripped version: black background, no axes/labels/legend
    p_stripped <- p +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.background  = ggplot2::element_rect(fill = stripped_bg, color = NA),
        panel.background = ggplot2::element_rect(fill = stripped_bg, color = NA)
      ) +
      ggplot2::guides(color = "none", fill = "none") +
      ggplot2::coord_fixed(ratio = 1)
    ggsave(stripped_path, p_stripped, width = width, height = height)
    log_message("  Saved: ", basename(stripped_path))

  }, error = function(e) {
    log_message("  WARN dual_save_plot failed (", basename(base_path), "): ",
                conditionMessage(e))
  })
}

set_global_seed <- function(seed) { set.seed(seed) }

# Save a ggplot in both PDF and PNG. Used by the F3-F5 compartment figure
# files where each panel lands as one file. Different from dual_save_plot
# (which writes labeled + stripped versions of one format) -- save_pdf_png
# writes one labeled version in two formats. PNG is rendered at 300 dpi for
# manuscript use.
save_pdf_png <- function(p, base_path, w = 8, h = 7) {
  if (is.null(p)) return(invisible(NULL))
  ensure_dir(dirname(base_path))
  ggplot2::ggsave(paste0(base_path, ".pdf"), p, width = w, height = h)
  ggplot2::ggsave(paste0(base_path, ".png"), p, width = w, height = h, dpi = 300)
  log_message("  Saved: ", basename(base_path), ".{pdf,png}")
}

# Canonical viz output subfolder scheme. One numbered bucket per analysis
# stage, applied uniformly to full / eye / compartment targets. Numbers are
# contiguous 01-13 so folders sort in figure order. Single source of truth:
# every viz write resolves its subfolder via viz_subdir() rather than a string
# literal, so the scheme can't drift across scripts again. See
# docs/plans/2026-06-16-viz-numbering-cleanup-design.md.
VIZ_BUCKETS <- c(
  qc           = "01_qc",
  integration  = "02_integration",
  celltypes    = "03_celltypes",
  markers      = "04_markers",
  dge          = "05_dge",
  escape       = "06_escape",
  milo         = "07_milo",
  repertoire   = "08_repertoire",
  tcr_motif    = "09_tcr_motif",
  lineage_arch = "10_lineage_arch",
  pca_coupling = "11_pca_coupling",
  composition  = "12_composition",
  bcr_lineage  = "13_bcr_lineage"
)

# Resolve a viz output subfolder under a target's viz_dir by bucket key.
# Does NOT create the directory; the save helpers (save_pdf_png /
# dual_save_plot) create the parent on write, so empty buckets never appear.
viz_subdir <- function(paths, key) {
  if (length(key) != 1L || !key %in% names(VIZ_BUCKETS))
    stop("Unknown viz bucket key: ", key,
         ". Valid: ", paste(names(VIZ_BUCKETS), collapse = ", "))
  file.path(paths$viz_dir, VIZ_BUCKETS[[key]])
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Flag T cells whose TCR carries the published HLA-B27-associated public
# signature: TRAV21 paired with a CDR3-beta matching .{5,9}[YF]S[TS].{4,8}.
# Both criteria must hold. Inputs follow the scRepertoire convention:
#   CTgene = "TRAV.TRAJ.TRAC_TRBV.TRBD.TRBJ.TRBC"
#   CTaa   = "TRAaminoacid_TRBaminoacid"
# Returns a logical vector the same length as the inputs; NAs are FALSE.
flag_hla_b27_pathogenic_tcr <- function(CTgene, CTaa) {
  n <- length(CTgene)
  if (n == 0L) return(logical(0))
  stopifnot(length(CTaa) == n)

  # TRA portion of CTgene -> first dot-separated token -> contains "TRAV21"
  tra_chain  <- stringr::str_split(CTgene, "_",   simplify = TRUE)[, 1]
  trav_token <- stringr::str_split(tra_chain, "\\.", simplify = TRUE)[, 1]
  has_trav21 <- grepl("TRAV21", trav_token, fixed = FALSE)

  # TRB portion of CTaa -> matches the [YF]S[TS] flank-constrained motif
  trb_aa     <- stringr::str_split(CTaa, "_", simplify = TRUE)[, 2]
  has_motif  <- grepl(".{5,9}[YF]S[TS].{4,8}", trb_aa)

  out <- has_trav21 & has_motif
  out[is.na(CTgene) | is.na(CTaa) | trb_aa == "" | trav_token == ""] <- FALSE
  out
}

# Resolve the celltype_broad column name in a Seurat-metadata data.frame.
# Eye-object carry-over renames `celltype_broad` -> `celltype_broad_full`
# whenever `cfg$eye_focus$rename_inherited_labels` is TRUE and the eye
# annotation step is disabled. Returns the column name to use, or NULL
# when neither variant is present.
resolve_celltype_broad <- function(meta) {
  if ("celltype_broad"      %in% colnames(meta)) return("celltype_broad")
  if ("celltype_broad_full" %in% colnames(meta)) return("celltype_broad_full")
  NULL
}

# Same convention for the fine `celltype` column.
resolve_celltype <- function(meta) {
  if ("celltype"      %in% colnames(meta)) return("celltype")
  if ("celltype_full" %in% colnames(meta)) return("celltype_full")
  NULL
}

# Resolve paths for the full ("all"), eye ("eye"), or compartment branches
# ("myeloid", "bcell", "tcell"). Returns a list with the same keys as cfg$paths
# so call sites can use it uniformly. Compartment branches inherit shared dirs
# (parent_dir, processed_dir, qc_dir, imgt_dir) from cfg$paths_eye and then
# from cfg$paths.
get_target_paths <- function(cfg, target = c("all", "eye", "myeloid", "bcell", "tcell")) {
  target <- match.arg(target)
  if (target == "all") return(cfg$paths)

  if (target == "eye") {
    p <- cfg$paths_eye
    if (is.null(p)) {
      stop("get_target_paths: target='eye' but cfg$paths_eye is not defined")
    }
  } else {
    key <- paste0("paths_", target)
    p <- cfg[[key]]
    if (is.null(p)) {
      stop("get_target_paths: target='", target, "' but cfg$", key, " is not defined")
    }
    # Compartment paths inherit from cfg$paths_eye for fields not set locally.
    if (!is.null(cfg$paths_eye)) {
      for (field in names(cfg$paths_eye)) {
        p[[field]] <- p[[field]] %||% cfg$paths_eye[[field]]
      }
    }
  }
  # Inherit project-wide dirs from the main paths block. qc_dir and
  # processed_dir are written once during ingest and shared by all branches.
  p$imgt_dir      <- p$imgt_dir      %||% cfg$paths$imgt_dir
  p$parent_dir    <- p$parent_dir    %||% cfg$paths$parent_dir
  p$qc_dir        <- p$qc_dir        %||% cfg$paths$qc_dir
  p$processed_dir <- p$processed_dir %||% cfg$paths$processed_dir
  p
}

# Make sure you have dplyr installed and loaded
# install.packages("dplyr")
library(dplyr)

#' Perform Enrichment Analysis using Fisher's Exact Test
#'
#' This function calculates the enrichment of a specific feature within clusters
#' using Fisher's Exact Test.
#'
#' @param data A data frame containing cluster assignments and feature annotations.
#' @param cluster_col A string with the name of the column containing cluster labels.
#' @param feature_col A string with the name of the column containing the binary feature
#'   (e.g., epitope specificity, cell type).
#' @param positive_value The value in `feature_col` that is considered the "positive"
#'   case for enrichment testing (e.g., "Yes", "Positive", "Specific").
#' @param output_file An optional string specifying the path to save the results CSV.
#'
#' @return A data frame containing the enrichment statistics (p-value, odds ratio, FDR)
#'   for each cluster, sorted by significance.
#'
cluster_enrichment <- function(data, cluster_col, feature_col, positive_value, output_file = NULL) {
  
  # --- Input Validation ---
  stopifnot(is.data.frame(data))
  if (!cluster_col %in% names(data)) {
    stop(paste("Cluster column '", cluster_col, "' not found in the data frame.", sep = ""))
  }
  if (!feature_col %in% names(data)) {
    stop(paste("Feature column '", feature_col, "' not found in the data frame.", sep = ""))
  }
  
  feature_values <- unique(data[[feature_col]])
  if (!positive_value %in% feature_values) {
    stop(paste("Positive value '", positive_value, "' not found in the feature column '", feature_col, "'.", sep = ""))
  }
  
  # Automatically determine the "negative" value
  negative_value <- setdiff(feature_values, positive_value)[1]
  if (is.na(negative_value)) {
    stop(paste("Could not determine a negative value in '", feature_col, "'. It should contain at least one value other than '", positive_value, "'.", sep=""))
  }
  
  message(paste0("Testing enrichment for '", positive_value, "' against '", negative_value, "'."))
  
  unique_clusters <- unique(data[[cluster_col]])
  
  # --- Loop over clusters to perform the test ---
  enrichment_results <- lapply(unique_clusters, function(cluster) {
    
    # Create the 2x2 contingency table
    in_cluster_yes <- sum(data[[cluster_col]] == cluster & data[[feature_col]] == positive_value)
    in_cluster_no  <- sum(data[[cluster_col]] == cluster & data[[feature_col]] == negative_value)
    
    not_in_cluster_yes <- sum(data[[cluster_col]] != cluster & data[[feature_col]] == positive_value)
    not_in_cluster_no  <- sum(data[[cluster_col]] != cluster & data[[feature_col]] == negative_value)
    
    contingency_matrix <- matrix(c(in_cluster_yes, not_in_cluster_yes,
                                   in_cluster_no, not_in_cluster_no),
                                 nrow = 2)
    
    # Perform the test to see if the odds ratio is greater than 1
    fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
    
    # Store the results for this cluster
    data.frame(
      cluster = cluster,
      p.value = fisher_result$p.value,
      odds.ratio = fisher_result$estimate
    )
  })
  
  # --- Combine and process results ---
  enrichment_df <- bind_rows(enrichment_results)
  
  # Rename the cluster column to match the input
  names(enrichment_df)[names(enrichment_df) == "cluster"] <- cluster_col
  
  # Calculate FDR and add a significance flag
  enrichment_df <- enrichment_df %>%
    mutate(FDR = p.adjust(p.value, method = "BH"),
           is_significant = FDR < 0.05) %>%
    arrange(FDR) # Sort results by FDR
  
  # --- Optionally save to file ---
  if (!is.null(output_file)) {
    write.csv(enrichment_df, output_file, row.names = FALSE)
    message(paste("Enrichment results saved to:", output_file))
  }
  
  return(enrichment_df)
}

# -----------------------------------------------------------------------------
# Per-cluster differential composition: Fisher's exact + binomial GLMM
# -----------------------------------------------------------------------------
# Two complementary tests, run for every cluster:
#   * Fisher's exact on the 2x2 (in-cluster vs not) x (groups[1] vs groups[2])
#     contingency table — treats cells as independent, returns OR + 95% CI.
#   * Binomial GLMM: cbind(n_in, n_out) ~ <contrast_col> + (1 | <sample_col>),
#     fit with lme4::glmer, LRT against the (1 | sample_col) null. Patient-aware
#     test that respects within-donor correlation.
#
# Args:
#   meta         : cell-level data.frame (one row per cell)
#   contrast_col : name of the column with the two-level grouping
#   groups       : c(reference, alternative) — order sets the OR/beta direction
#   cluster_col  : name of the column with cluster labels
#   sample_col   : donor/subject identifier (default "Subject")
#
# Returns: data.frame, one row per cluster, with columns:
#   cluster, group_ref, group_alt,
#   n_g1_in/out, n_g2_in/out,
#   OR_fisher, OR_lo, OR_hi, p_fisher, q_fisher,
#   beta_glmm, se_glmm, p_glmm, q_glmm, glmm_status
#
# Notes:
#   * BH adjustment runs across clusters within a single call.
#   * GLMM may return p_glmm = NA when fits do not converge, or
#     glmm_status = "singular" when the random-effect variance collapses
#     (still usable, but flag for review).
#   * With small donor counts the GLMM is conservative; pair it with Fisher
#     for sanity-checking direction of effect.
run_fisher_glmm_per_cluster <- function(meta,
                                        contrast_col,
                                        groups,
                                        cluster_col,
                                        sample_col = "Subject") {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("lme4 is required for run_fisher_glmm_per_cluster.")
  }
  stopifnot(length(groups) == 2,
            contrast_col %in% colnames(meta),
            cluster_col %in% colnames(meta),
            sample_col  %in% colnames(meta))

  meta <- meta[meta[[contrast_col]] %in% groups &
                 !is.na(meta[[cluster_col]]) &
                 !is.na(meta[[sample_col]]), , drop = FALSE]
  meta[[contrast_col]] <- factor(as.character(meta[[contrast_col]]),
                                 levels = groups)
  meta[[sample_col]]   <- as.character(meta[[sample_col]])
  meta[[cluster_col]]  <- as.character(meta[[cluster_col]])

  totals <- aggregate(rep(1, nrow(meta)),
                      by = list(meta[[sample_col]],
                                meta[[contrast_col]]),
                      FUN = sum)
  names(totals) <- c(sample_col, contrast_col, "n_total")
  totals[[contrast_col]] <- factor(as.character(totals[[contrast_col]]),
                                   levels = groups)

  ctrl <- lme4::glmerControl(optimizer = "bobyqa",
                             optCtrl = list(maxfun = 2e5))
  fit_safe <- function(formula, data) {
    tryCatch(suppressWarnings(suppressMessages(
      lme4::glmer(formula, data = data,
                  family = stats::binomial(), control = ctrl))),
      error = function(e) NULL)
  }

  per_cluster <- function(cl) {
    in_clust <- meta[[cluster_col]] == cl
    grp <- meta[[contrast_col]]
    a <- sum(in_clust  & grp == groups[1])
    b <- sum(!in_clust & grp == groups[1])
    c <- sum(in_clust  & grp == groups[2])
    d <- sum(!in_clust & grp == groups[2])

    fish <- tryCatch(
      fisher.test(matrix(c(a, c, b, d), nrow = 2)),
      error = function(e) NULL)

    in_per_sub <- if (any(in_clust)) {
      tmp <- aggregate(rep(1, sum(in_clust)),
                       by = list(meta[[sample_col]][in_clust],
                                 meta[[contrast_col]][in_clust]),
                       FUN = sum)
      names(tmp) <- c(sample_col, contrast_col, "n_in")
      tmp
    } else {
      data.frame(setNames(list(character(0), factor(character(0),
                                                    levels = groups),
                                integer(0)),
                          c(sample_col, contrast_col, "n_in")))
    }
    sub_df <- merge(totals, in_per_sub, all.x = TRUE,
                    by = c(sample_col, contrast_col))
    sub_df$n_in[is.na(sub_df$n_in)] <- 0
    sub_df$n_out <- sub_df$n_total - sub_df$n_in
    sub_df[[contrast_col]] <- factor(as.character(sub_df[[contrast_col]]),
                                     levels = groups)

    glmm_status <- "ok"
    p_glmm <- NA_real_; beta_glmm <- NA_real_; se_glmm <- NA_real_

    enough <- length(unique(sub_df[[sample_col]])) >= 2 &&
              length(unique(sub_df[[contrast_col]])) == 2 &&
              sum(sub_df$n_in) > 0
    if (!enough) {
      glmm_status <- "insufficient_data"
    } else {
      f1 <- stats::as.formula(paste0("cbind(n_in, n_out) ~ ", contrast_col,
                                     " + (1 | ", sample_col, ")"))
      f0 <- stats::as.formula(paste0("cbind(n_in, n_out) ~ 1 + (1 | ",
                                     sample_col, ")"))
      m1 <- fit_safe(f1, sub_df)
      m0 <- fit_safe(f0, sub_df)
      if (is.null(m1) || is.null(m0)) {
        glmm_status <- "fit_failed"
      } else {
        lrt <- tryCatch(stats::anova(m0, m1, test = "Chisq"),
                        error = function(e) NULL)
        if (!is.null(lrt)) p_glmm <- lrt[["Pr(>Chisq)"]][2]
        cs <- tryCatch(summary(m1)$coefficients, error = function(e) NULL)
        if (!is.null(cs) && nrow(cs) >= 2) {
          beta_glmm <- cs[2, "Estimate"]
          se_glmm   <- cs[2, "Std. Error"]
        }
        if (isTRUE(lme4::isSingular(m1))) glmm_status <- "singular"
      }
    }

    data.frame(
      cluster     = cl,
      group_ref   = groups[1], group_alt = groups[2],
      n_g1_in     = a, n_g1_out = b,
      n_g2_in     = c, n_g2_out = d,
      OR_fisher   = if (!is.null(fish)) unname(fish$estimate) else NA_real_,
      OR_lo       = if (!is.null(fish)) fish$conf.int[1]      else NA_real_,
      OR_hi       = if (!is.null(fish)) fish$conf.int[2]      else NA_real_,
      p_fisher    = if (!is.null(fish)) fish$p.value           else NA_real_,
      beta_glmm   = beta_glmm,
      se_glmm     = se_glmm,
      p_glmm      = p_glmm,
      glmm_status = glmm_status,
      stringsAsFactors = FALSE
    )
  }

  clusters <- sort(unique(meta[[cluster_col]]))
  out <- do.call(rbind, lapply(clusters, per_cluster))
  out$q_fisher <- p.adjust(out$p_fisher, method = "BH")
  out$q_glmm   <- p.adjust(out$p_glmm,   method = "BH")
  out
}

# substate_labels formats a vector of cluster IDs as "<id>: <label>", reading
# the human-readable label from cfg$compartment_substate_labels[[cmp]][[id]].
# When a label is missing, falls back to "<id>: cluster_<id>". Used by the
# per-figure visualization files to keep cluster numbers visible alongside
# hand-curated names.
substate_labels <- function(cfg, cmp, ids) {
  cfg_labels <- cfg$compartment_substate_labels[[cmp]]
  vapply(as.character(ids), function(id) {
    if (!is.null(cfg_labels) && !is.null(cfg_labels[[id]])) {
      paste0(id, ": ", cfg_labels[[id]])
    } else {
      paste0(id, ": cluster_", id)
    }
  }, character(1), USE.NAMES = FALSE)
}
