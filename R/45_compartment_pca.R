# R/45_compartment_pca.R
# Per-compartment pseudobulk PCA. Refactored from the inline implementation
# in R/85_viz_myeloid.R:.myeloid_panelD_pca_facets so the same logic powers
# both F3 panel F (myeloid) and F4 panels D/E (T cell), and so subject-level
# PC1 scores live as a canonical on-disk CSV that the cross-compartment
# bridge (R/46) and LIANA gene-set construction (R/48) can read directly.
#
# Public entry points:
#   run_compartment_pca(cfg, target)         orchestrates one compartment
#   compute_per_substate_pca(obj, ...)       does the actual aggregation +
#                                             vst + prcomp + sign-orient,
#                                             returns scores / loadings
#                                             tibbles ready to write.
#   .pc1_loadings_by_program(target, programs, cfg)
#                                             pulls per-program PC1 loadings
#                                             restricted to disease-separating
#                                             substates and writes
#                                             pc1_loadings_by_program.csv.
#
# Outputs under outputs/tables/eye/<target>/:
#   pca_subject_scores.csv         one row per (substate, pseudobulk sample)
#   pca_gene_loadings.csv          one row per (substate, gene)
#   pca_variance_explained.csv     one row per (substate, PC)
#   pca_pc1_significance.csv       one row per substate (Welch t + BH q)
#   pc1_loadings_by_program.csv    one row per (substate, program, gene)
#
# Compatibility note: matches the F3 inline implementation's conventions
# exactly: DESeq2::vst(blind = FALSE, design = ~ group), HVG top-2000 by
# rowVars, prcomp(scale. = TRUE), PC1 sign-flipped so the Viral centroid is
# at positive PC1. This means the new module's PC1 scores reproduce the
# (formerly inline) F3 panel F values up to numerical noise.
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Asymmetric pseudobulk floor. Myeloid is sparser per (subject, substate) so
# it gets a lower floor than T / B cell. Sourced from cfg$compartment_pca to
# stay overrideable from config without code edits.
.pca_min_cells <- function(cfg, target) {
  cpcfg <- cfg$compartment_pca %||% list()
  key <- paste0("min_cells_per_subject_substate_", target)
  as.integer(cpcfg[[key]] %||% 20L)
}

# Heart of the module. Given a Seurat compartment object, run per-substate
# pseudobulk PCA following the F3 convention. Returns a list with elements
#   scores      tibble (substate, sample, n_cells, group, PC1..PC5, PC1_oriented)
#   loadings    tibble (substate, gene, PC1..PC5, PC1_oriented, loading_rank_within_substate)
#   variance    tibble (substate, PC, var_explained, cum_var)
#   significance  tibble (substate, t_statistic, df, p_value, q_value, n_NIU, n_Viral, separating)
# Callers (run_compartment_pca, .myeloid_panelD_pca_facets) handle the writing
# and the plotting.
compute_per_substate_pca <- function(obj,
                                     cluster_col       = "knn.leiden.cluster",
                                     group_col         = "Phenotype_2",
                                     subject_col       = "Subject",
                                     sample_col        = NULL,
                                     groups            = c("NIU", "Viral"),
                                     min_cells_per_pb  = 30L,
                                     min_gene_count    = 10L,
                                     hvg_n             = 2000L,
                                     n_pcs             = 5L,
                                     vst_blind         = FALSE,
                                     pc1_split_fdr     = 0.05) {
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("compute_per_substate_pca: DESeq2 required.")
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("compute_per_substate_pca: matrixStats required.")

  if (is.null(sample_col)) {
    sample_col <- if ("Subject_Timepoint" %in% colnames(obj[[]]))
                    "Subject_Timepoint" else "orig.ident"
  }

  pbs <- build_per_substate_pseudobulks(
    obj,
    cluster_col      = cluster_col,
    group_col        = group_col,
    groups           = groups,
    min_cells_per_pb = min_cells_per_pb
  )
  if (length(pbs) == 0) {
    log_message("  compute_per_substate_pca: no pseudobulks after floor=",
                min_cells_per_pb, " filter; aborting.")
    return(NULL)
  }

  # Subject lookup keyed by the sample column used in the pseudobulks. When
  # sample_col is Subject_Timepoint, this collapses "<subject>_<visit>" to
  # "<subject>". Downstream bridge (R/46) further averages across timepoints.
  meta <- obj[[]]
  sample_to_subject <- tapply(as.character(meta[[subject_col]]),
                              as.character(meta[[sample_col]]),
                              function(x) x[1])
  sample_to_pheno   <- tapply(as.character(meta[[group_col]]),
                              as.character(meta[[sample_col]]),
                              function(x) x[1])
  sample_to_eti     <- if ("Etiology" %in% colnames(meta))
                         tapply(as.character(meta$Etiology),
                                as.character(meta[[sample_col]]),
                                function(x) x[1]) else NULL
  sample_to_gran    <- if ("Phenotype" %in% colnames(meta))
                         tapply(as.character(meta$Phenotype),
                                as.character(meta[[sample_col]]),
                                function(x) x[1]) else NULL
  sample_to_cohort  <- if ("Cohort" %in% colnames(meta))
                         tapply(as.character(meta$Cohort),
                                as.character(meta[[sample_col]]),
                                function(x) x[1]) else NULL

  scores_rows   <- list()
  loadings_rows <- list()
  variance_rows <- list()
  sig_rows      <- list()

  for (ck in names(pbs)) {
    pb <- pbs[[ck]]
    if (inherits(pb, "SummarizedExperiment")) {
      cd <- SummarizedExperiment::colData(pb)
      m  <- SummarizedExperiment::assay(pb, "counts")
    } else {
      cd <- pb$coldata; m <- pb$counts
    }
    grp <- factor(cd$group, levels = groups)
    if (ncol(m) < 4 || length(unique(grp)) < 2) {
      log_message("  PCA substate ", ck, ": <4 columns or single group; skipping.")
      next
    }

    pca_res <- tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(m),
        colData   = data.frame(group = grp),
        design    = ~ group)
      dds <- dds[rowSums(DESeq2::counts(dds)) > min_gene_count, ]
      vsd <- DESeq2::vst(dds, blind = vst_blind)
      mat <- SummarizedExperiment::assay(vsd)
      vars <- matrixStats::rowVars(mat)
      keep <- order(-vars)[seq_len(min(hvg_n, length(vars)))]
      mat  <- mat[keep, ]
      stats::prcomp(t(mat), scale. = TRUE)
    }, error = function(e) {
      log_message("    PCA failed for substate ", ck, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(pca_res)) next

    # Sign orient so positive PC1 is the Viral centroid (matches F3 convention).
    grp_chr   <- as.character(grp)
    viral_mean <- mean(pca_res$x[grp_chr == "Viral", 1], na.rm = TRUE)
    niu_mean   <- mean(pca_res$x[grp_chr == "NIU",   1], na.rm = TRUE)
    flip <- if (is.finite(viral_mean) && is.finite(niu_mean) &&
                viral_mean < niu_mean) -1 else 1
    pc1_oriented <- pca_res$x[, 1] * flip
    rotation_pc1 <- pca_res$rotation[, 1] * flip

    n_keep_pc <- min(n_pcs, ncol(pca_res$x))
    scores_mat <- pca_res$x[, seq_len(n_keep_pc), drop = FALSE]
    colnames(scores_mat) <- paste0("PC", seq_len(n_keep_pc))
    rot_mat <- pca_res$rotation[, seq_len(n_keep_pc), drop = FALSE]
    colnames(rot_mat) <- paste0("PC", seq_len(n_keep_pc))

    sample_id <- as.character(cd$sample)
    n_cells   <- as.integer(cd$n_cells)
    scores_df <- tibble::tibble(
      substate     = as.character(ck),
      sample       = sample_id,
      subject      = unname(sample_to_subject[sample_id]),
      Phenotype_2  = grp_chr,
      Etiology     = if (!is.null(sample_to_eti))    unname(sample_to_eti[sample_id])    else NA_character_,
      Phenotype    = if (!is.null(sample_to_gran))   unname(sample_to_gran[sample_id])   else NA_character_,
      Cohort       = if (!is.null(sample_to_cohort)) unname(sample_to_cohort[sample_id]) else NA_character_,
      n_cells      = n_cells
    )
    scores_df <- dplyr::bind_cols(scores_df,
                                  tibble::as_tibble(scores_mat),
                                  PC1_oriented = pc1_oriented)

    rot_df <- tibble::tibble(
      substate = as.character(ck),
      gene     = rownames(rot_mat)
    )
    rot_df <- dplyr::bind_cols(rot_df,
                               tibble::as_tibble(rot_mat),
                               PC1_oriented = rotation_pc1)
    rot_df$loading_rank_within_substate <- rank(-abs(rot_df$PC1_oriented),
                                                ties.method = "first")

    var_pct <- (pca_res$sdev[seq_len(n_keep_pc)]^2 / sum(pca_res$sdev^2)) * 100
    var_df <- tibble::tibble(
      substate      = as.character(ck),
      PC            = paste0("PC", seq_len(n_keep_pc)),
      var_explained = round(var_pct, 3),
      cum_var       = round(cumsum(var_pct), 3)
    )

    # Welch t on PC1_oriented across the two groups. Underpowered substates
    # (n < 3 per group) get NA p-values rather than t-test errors.
    pc1_niu   <- pc1_oriented[grp_chr == "NIU"]
    pc1_viral <- pc1_oriented[grp_chr == "Viral"]
    sig_row <- tibble::tibble(
      substate     = as.character(ck),
      n_NIU        = length(pc1_niu),
      n_Viral      = length(pc1_viral),
      mean_NIU     = if (length(pc1_niu)   > 0) mean(pc1_niu,   na.rm = TRUE) else NA_real_,
      mean_Viral   = if (length(pc1_viral) > 0) mean(pc1_viral, na.rm = TRUE) else NA_real_,
      t_statistic  = NA_real_,
      df           = NA_real_,
      p_value      = NA_real_
    )
    if (length(pc1_niu) >= 3 && length(pc1_viral) >= 3) {
      tt <- tryCatch(stats::t.test(pc1_viral, pc1_niu, var.equal = FALSE),
                     error = function(e) NULL)
      if (!is.null(tt)) {
        sig_row$t_statistic <- unname(tt$statistic)
        sig_row$df          <- unname(tt$parameter)
        sig_row$p_value     <- tt$p.value
      }
    }

    scores_rows[[ck]]   <- scores_df
    loadings_rows[[ck]] <- rot_df
    variance_rows[[ck]] <- var_df
    sig_rows[[ck]]      <- sig_row
  }

  if (length(scores_rows) == 0) return(NULL)

  scores   <- dplyr::bind_rows(scores_rows)
  loadings <- dplyr::bind_rows(loadings_rows)
  variance <- dplyr::bind_rows(variance_rows)
  sig      <- dplyr::bind_rows(sig_rows)
  sig$q_value    <- stats::p.adjust(sig$p_value, method = "BH")
  sig$separating <- !is.na(sig$q_value) & sig$q_value < pc1_split_fdr

  list(scores = scores, loadings = loadings,
       variance = variance, significance = sig)
}

# Entry point called from run_pipeline.R Phase 1d. Reads the compartment
# Seurat object from disk, runs compute_per_substate_pca with the
# compartment-specific floor, writes four CSVs, then optionally calls
# .pc1_loadings_by_program when cfg$<target>_programs is defined.
run_compartment_pca <- function(cfg,
                                target = c("myeloid", "tcell", "bcell")) {
  target <- match.arg(target)
  paths <- get_target_paths(cfg, target)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("compartment_pca[", target, "]: object missing at ",
                obj_path, "; skipping.")
    return(invisible(FALSE))
  }
  log_message("=== compartment_pca [", target, "] ===")
  obj <- readRDS(obj_path)

  min_cells <- .pca_min_cells(cfg, target)
  cpcfg     <- cfg$compartment_pca %||% list()
  res <- compute_per_substate_pca(
    obj,
    min_cells_per_pb = min_cells,
    min_gene_count   = as.integer(cpcfg$min_gene_count %||% 10L),
    hvg_n            = as.integer(cpcfg$hvg_n %||% 2000L),
    n_pcs            = as.integer(cpcfg$n_pcs %||% 5L),
    vst_blind        = isTRUE(cpcfg$vst_blind),
    pc1_split_fdr    = as.numeric(cpcfg$pc1_split_fdr %||% 0.05)
  )
  if (is.null(res)) {
    log_message("compartment_pca[", target, "]: PCA returned nothing.")
    return(invisible(FALSE))
  }

  ensure_dir(paths$results_tables)
  utils::write.csv(res$scores,
                   file.path(paths$results_tables, "pca_subject_scores.csv"),
                   row.names = FALSE)
  utils::write.csv(res$loadings,
                   file.path(paths$results_tables, "pca_gene_loadings.csv"),
                   row.names = FALSE)
  utils::write.csv(res$variance,
                   file.path(paths$results_tables, "pca_variance_explained.csv"),
                   row.names = FALSE)
  utils::write.csv(res$significance,
                   file.path(paths$results_tables, "pca_pc1_significance.csv"),
                   row.names = FALSE)
  log_message("compartment_pca[", target, "]: wrote ",
              nrow(res$scores), " sample rows, ",
              nrow(res$loadings), " gene rows, ",
              sum(res$significance$separating, na.rm = TRUE),
              " separating substate(s).")

  # Per-program loadings (Figure 4 panel E). Skipped when the program list
  # is not defined in config — the bridge / LIANA / NicheNet steps don't
  # depend on it.
  programs <- cfg[[paste0(target, "_programs")]]
  if (!is.null(programs) && length(programs) > 0) {
    .pc1_loadings_by_program(target, programs, cfg,
                             scores = res$scores,
                             loadings = res$loadings,
                             significance = res$significance)
  }

  invisible(TRUE)
}

# Restrict PC1 loadings to a curated set of program gene panels. Used for
# F3 panel E (myeloid programs: HLA_I, HLA_II, CD1, ...) and F4 panel E
# (tcell programs: TCR_signal, Checkpoint, ...). One row per (substate,
# program, gene) keeping only separating substates and only genes present
# in the substate's loadings table.
.pc1_loadings_by_program <- function(target, programs, cfg,
                                     scores, loadings, significance) {
  paths <- get_target_paths(cfg, target)
  sep_subs <- significance$substate[significance$separating]
  if (length(sep_subs) == 0L) {
    log_message("  pc1_loadings_by_program[", target,
                "]: no separating substates; writing empty stub.")
    out <- tibble::tibble(substate = character(), program = character(),
                          gene = character(), PC1_oriented = double(),
                          loading_rank_within_program = integer())
    utils::write.csv(out,
                     file.path(paths$results_tables,
                               "pc1_loadings_by_program.csv"),
                     row.names = FALSE)
    return(invisible(out))
  }
  rows <- list()
  for (prog in names(programs)) {
    genes <- as.character(programs[[prog]])
    for (ck in sep_subs) {
      df <- dplyr::filter(loadings,
                          .data$substate == ck,
                          .data$gene %in% genes)
      if (nrow(df) == 0L) next
      df <- df |>
        dplyr::mutate(program = prog,
                      loading_rank_within_program =
                        rank(-abs(.data$PC1_oriented),
                             ties.method = "first")) |>
        dplyr::select(substate, program, gene, PC1_oriented,
                      loading_rank_within_program)
      rows[[paste(ck, prog, sep = "::")]] <- df
    }
  }
  out <- if (length(rows) == 0L) {
    tibble::tibble(substate = character(), program = character(),
                   gene = character(), PC1_oriented = double(),
                   loading_rank_within_program = integer())
  } else dplyr::bind_rows(rows)

  utils::write.csv(out,
                   file.path(paths$results_tables,
                             "pc1_loadings_by_program.csv"),
                   row.names = FALSE)
  log_message("  pc1_loadings_by_program[", target, "]: wrote ",
              nrow(out), " (substate, program, gene) rows across ",
              length(sep_subs), " separating substate(s).")
  invisible(out)
}
