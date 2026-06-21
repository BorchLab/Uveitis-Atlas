# R/48_nichenet_myeloid_tcell.R
# NicheNet ligand-activity prediction (Figure 4 panel H, supplementary by default).
#
# LIANA tests ligand-receptor coexpression; NicheNet tests whether downstream
# transcriptional changes in receivers are explained by predicted ligand
# activity through prior signaling networks. The two answer different questions;
# the methods-paper-strength claim is when both converge on the same ligands.
#
# This step is gated off by default (`cfg$steps$nichenet_myeloid_tcell: false`)
# — flip on after LIANA results are reviewed so the panel earns its place.
#
# Inputs:
#   outputs/tables/eye/tcell/pca_subject_scores.csv
#   outputs/tables/eye/tcell/pca_gene_loadings.csv
#   outputs/tables/eye/tcell/pca_pc1_significance.csv
#   compartment Seurat objects (for expression filtering)
# Outputs (under outputs/tables/cross_compartment/):
#   nichenet_NIU_pole_ligands.csv          top ligands explaining NIU-pole T cell genes
#   nichenet_Viral_pole_ligands.csv        top ligands explaining Viral-pole T cell genes
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Per-substate fraction of cells expressing each gene. Returns a named numeric
# vector with row gene names. Filters to genes that are present on the assay.
.fraction_expressed_per_substate <- function(obj, substate_col = "substate_key",
                                             threshold = 0.10) {
  DefaultAssay(obj) <- "RNA"
  e <- GetAssayData(obj, assay = "RNA", layer = "data")
  substates <- as.character(obj[[substate_col, drop = TRUE]])
  if (length(unique(substates)) == 0L) return(list(by_substate = NULL,
                                                   any_substate = character(0)))
  per_sub <- lapply(unique(substates), function(s) {
    cells <- which(substates == s)
    Matrix::rowMeans(e[, cells, drop = FALSE] > 0)
  })
  names(per_sub) <- unique(substates)
  any_sub <- Reduce("pmax", per_sub)
  list(by_substate = per_sub,
       any_substate = names(any_sub)[any_sub >= threshold])
}

# Per-substate-and-pole geneset construction. Returns a tibble (substate,
# pole, gene, PC1_oriented).
.nichenet_geneset <- function(loadings, separating_substates, top_n = 50L) {
  rows <- list()
  for (ck in separating_substates) {
    sub <- dplyr::filter(loadings, .data$substate == ck)
    if (nrow(sub) == 0L) next
    pos <- sub |> dplyr::slice_max(.data$PC1_oriented, n = top_n,
                                   with_ties = FALSE)
    neg <- sub |> dplyr::slice_min(.data$PC1_oriented, n = top_n,
                                   with_ties = FALSE)
    rows[[paste0(ck, "_Viral")]] <- pos |>
      dplyr::mutate(substate = ck, pole = "Viral") |>
      dplyr::select("substate", "pole", "gene", "PC1_oriented")
    rows[[paste0(ck, "_NIU")]] <- neg |>
      dplyr::mutate(substate = ck, pole = "NIU") |>
      dplyr::select("substate", "pole", "gene", "PC1_oriented")
  }
  if (length(rows) == 0L) return(tibble::tibble())
  dplyr::bind_rows(rows)
}

# Main entry called from run_pipeline.R Phase 1d. Skips silently when
# nichenetr is missing or the input CSVs aren't on disk.
run_nichenet_myeloid_to_tcell <- function(cfg) {
  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    log_message("nichenet: package not installed; skipping.")
    return(invisible(FALSE))
  }
  ncfg <- cfg$nichenet %||% list()
  top_n_ligands  <- as.integer(ncfg$top_n_ligands %||% 20L)
  geneset_size   <- as.integer(ncfg$geneset_size %||% 50L)
  min_expr_frac  <- as.numeric(ncfg$min_expression_fraction %||% 0.10)

  paths_myel <- get_target_paths(cfg, "myeloid")
  paths_tcel <- get_target_paths(cfg, "tcell")
  cc_paths   <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/",
    viz    = "outputs/viz/cross_compartment/")
  ensure_dir(cc_paths$tables)

  load_csv <- file.path(paths_tcel$results_tables, "pca_gene_loadings.csv")
  sig_csv  <- file.path(paths_tcel$results_tables, "pca_pc1_significance.csv")
  if (!file.exists(load_csv) || !file.exists(sig_csv)) {
    log_message("nichenet: PCA loadings or significance CSV missing on T cell ",
                "(", load_csv, "); run compartment_pca first.")
    return(invisible(FALSE))
  }
  loadings <- utils::read.csv(load_csv, stringsAsFactors = FALSE)
  sig      <- utils::read.csv(sig_csv,  stringsAsFactors = FALSE)
  sep_flag <- as.logical(sig$separating)
  separating <- as.character(sig$substate[!is.na(sep_flag) & sep_flag])
  if (length(separating) == 0L) {
    log_message("nichenet: no T cell substates flagged as separating; ",
                "skipping. Inspect pca_pc1_significance.csv.")
    return(invisible(FALSE))
  }
  geneset <- .nichenet_geneset(loadings, separating, top_n = geneset_size)
  if (nrow(geneset) == 0L) {
    log_message("nichenet: no genes resolved from PCA loadings; skipping.")
    return(invisible(FALSE))
  }

  myel_path <- file.path(paths_myel$results_objects, "IntegratedSeuratObject.rds")
  tcel_path <- file.path(paths_tcel$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(myel_path) || !file.exists(tcel_path)) {
    log_message("nichenet: compartment objects missing; skipping.")
    return(invisible(FALSE))
  }
  myel <- readRDS(myel_path)
  tcel <- readRDS(tcel_path)
  myel$substate_key <- get_substate_key_vector(myel, "myeloid")
  tcel$substate_key <- get_substate_key_vector(tcel, "tcell")

  log_message("nichenet: pulling prior matrix + LR network")
  ltm <- tryCatch(nichenetr::ligand_target_matrix,
                  error = function(e) NULL)
  lrn <- tryCatch(nichenetr::lr_network,
                  error = function(e) NULL)
  if (is.null(ltm) || is.null(lrn)) {
    log_message("  nichenet: prior matrices not available from nichenetr; ",
                "ensure the package data is downloaded (data(...) call) ",
                "or set cfg$nichenet$prior_path to local copies.")
    return(invisible(FALSE))
  }

  myel_expr <- .fraction_expressed_per_substate(myel, "substate_key",
                                                threshold = min_expr_frac)
  tcel_expr <- .fraction_expressed_per_substate(tcel, "substate_key",
                                                threshold = min_expr_frac)
  candidate_ligands <- intersect(myel_expr$any_substate, unique(lrn$from))
  if (length(candidate_ligands) == 0L) {
    log_message("nichenet: no candidate ligands cleared expression filter.")
    return(invisible(FALSE))
  }

  # Helper to compute the dominant myeloid sender per ligand (highest
  # per-substate fraction) so the output table flags the most likely source.
  myel_per_sub <- myel_expr$by_substate
  dominant_sender <- function(g) {
    if (is.null(myel_per_sub) || length(myel_per_sub) == 0L) return(c(NA, NA))
    fr <- vapply(myel_per_sub, function(v) if (g %in% names(v)) v[[g]] else 0,
                 numeric(1))
    if (max(fr, na.rm = TRUE) <= 0) return(c(NA, NA))
    nm <- names(fr)[which.max(fr)]
    c(nm, fr[[nm]])
  }

  out_rows <- list()
  for (ck in separating) {
    bg_genes <- tcel_expr$by_substate[[ck]]
    if (is.null(bg_genes)) next
    background <- intersect(names(bg_genes)[bg_genes >= min_expr_frac],
                            rownames(ltm))
    for (pole in c("NIU", "Viral")) {
      gs <- geneset |>
        dplyr::filter(.data$substate == ck, .data$pole == pole) |>
        dplyr::pull("gene")
      gs <- intersect(gs, rownames(ltm))
      if (length(gs) < 5L || length(background) < 50L) {
        log_message("  nichenet[", ck, "/", pole, "]: geneset=", length(gs),
                    " background=", length(background),
                    " — too small to score.")
        next
      }
      activities <- tryCatch(
        nichenetr::predict_ligand_activities(
          geneset            = gs,
          background_expressed_genes = background,
          ligand_target_matrix       = ltm,
          potential_ligands          = intersect(candidate_ligands,
                                                  colnames(ltm))),
        error = function(e) {
          log_message("    nichenet[", ck, "/", pole, "] failed: ",
                      conditionMessage(e)); NULL })
      if (is.null(activities) || nrow(activities) == 0L) next
      activities <- activities |>
        tibble::as_tibble() |>
        dplyr::arrange(dplyr::desc(.data$pearson)) |>
        dplyr::slice_head(n = top_n_ligands)
      sender <- t(vapply(activities$test_ligand, dominant_sender, character(2)))
      activities$sending_myeloid_substate_top      <- sender[, 1]
      activities$sending_myeloid_substate_fraction <- as.numeric(sender[, 2])
      activities$tcell_substate <- ck
      activities$pole           <- pole
      out_rows[[paste(ck, pole, sep = "::")]] <- activities |>
        dplyr::rename(ligand = "test_ligand") |>
        dplyr::select("tcell_substate", "pole", "ligand", "pearson",
                      "auroc", "aupr",
                      "sending_myeloid_substate_top",
                      "sending_myeloid_substate_fraction")
    }
  }

  if (length(out_rows) == 0L) {
    log_message("nichenet: no successful scoring; skipping write.")
    return(invisible(FALSE))
  }
  combined <- dplyr::bind_rows(out_rows)
  niu_path   <- file.path(cc_paths$tables, "nichenet_NIU_pole_ligands.csv")
  viral_path <- file.path(cc_paths$tables, "nichenet_Viral_pole_ligands.csv")
  utils::write.csv(dplyr::filter(combined, .data$pole == "NIU"),
                   niu_path, row.names = FALSE)
  utils::write.csv(dplyr::filter(combined, .data$pole == "Viral"),
                   viral_path, row.names = FALSE)
  log_message("nichenet: wrote ",
              nrow(combined), " ligand rows across ",
              length(separating), " separating T cell substate(s).")
  invisible(TRUE)
}
