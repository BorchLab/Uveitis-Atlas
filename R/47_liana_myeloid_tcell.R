# R/47_liana_myeloid_tcell.R
# Myeloid -> T cell cell-cell communication via LIANA consensus rank.
# Replaces the prior CellChat module (R/64_cellchat.R, retired 2026-05-19).
# Why LIANA: rank-aggregating across NATMI, Connectome, logFC, CellPhoneDB,
# and SCA reduces single-method false positives. CellChat alone overcalls.
#
# Outputs under outputs/tables/cross_compartment/:
#   liana_myeloid_to_tcell_NIU.csv         per-condition LR rank table
#   liana_myeloid_to_tcell_Viral.csv       per-condition LR rank table
#   liana_myeloid_to_tcell_combined.csv    joined NIU + Viral with disease_bias
#   liana_cluster_n_audit.csv              skipped clusters / cell counts
#
# Naming:
#   ccc_label is built from `substate_key` (myeloid_<id> / tcell_<id>) when
#   present. If a compartment object pre-dates the R/22 stamp pass, the
#   module falls back to "<target>_<cluster_id>" reconstructed from
#   knn.leiden.cluster.
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Substate-key resolver lives in R/23_substate_labels.R::get_substate_key_vector.
# Same helper is used by R/48_nichenet_myeloid_tcell.

# Classify an LR pair by which `myeloid_programs` family the ligand belongs to.
# Returns one of the family names from cfg$myeloid_programs (HLA_I, HLA_II,
# CD1, Costim_ligand, Checkpoint_ligand, Cytokine, Chemokine) or "Other" when
# the ligand isn't in any curated panel. ligand_complex may be a single gene
# (CD80) or a complex ("CD86_CD80" / "CD80&CD86" depending on LIANA version) —
# split on common separators and match any subunit.
.liana_classify_ligand_family <- function(ligands, cfg) {
  programs <- cfg$myeloid_programs %||% list()
  if (length(programs) == 0L) return(rep(NA_character_, length(ligands)))
  # Reverse map: gene -> family. If a gene is in multiple families, first
  # family wins (config order: HLA_I, HLA_II, CD1, Costim, Checkpoint,
  # Cytokine, Chemokine — biologically distinct so collisions are rare).
  fam_lookup <- character(0)
  for (fam in names(programs)) {
    g <- as.character(programs[[fam]])
    new <- setNames(rep(fam, length(g)), g)
    fam_lookup <- c(fam_lookup, new[!names(new) %in% names(fam_lookup)])
  }
  split_subunits <- function(s) {
    if (is.na(s) || nchar(s) == 0L) return(character(0))
    unlist(strsplit(s, "[_&+]"))
  }
  vapply(as.character(ligands), function(s) {
    parts <- split_subunits(s)
    hit <- fam_lookup[parts]
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0L) "Other" else unname(hit[1])
  }, character(1), USE.NAMES = FALSE)
}

# Audit per-(condition, ccc_label) cell counts and return the labels that
# clear the min-cells threshold in BOTH conditions (so LIANA has them on
# both sides of the bias comparison).
.liana_cell_count_audit <- function(meta, condition_col, conditions,
                                    min_cells, audit_path) {
  audit <- meta |>
    dplyr::filter(.data[[condition_col]] %in% conditions) |>
    dplyr::count(.data[[condition_col]], .data$ccc_label,
                 name = "n_cells") |>
    dplyr::mutate(passes = .data$n_cells >= min_cells)
  utils::write.csv(audit, audit_path, row.names = FALSE)
  log_message("  liana cluster audit: ", audit_path)
  wide <- audit |>
    tidyr::pivot_wider(id_cols = "ccc_label",
                       names_from = !!rlang::sym(condition_col),
                       values_from = "passes",
                       values_fill = FALSE)
  keep <- wide$ccc_label[apply(as.matrix(wide[, -1]), 1, all)]
  as.character(keep)
}

# Per-condition LIANA wrapper. Pulls LIANA's modular interface and aggregates
# to consensus rank. Returns the long-format tibble keyed by (source, target,
# ligand_complex, receptor_complex) with `aggregate_rank` and `mean_rank`.
#
# Note (2026-05-19): LIANA's Seurat code path calls
# `Seurat::GetAssayData(..., slot = ...)`, which is defunct in
# SeuratObject 5.0+. We convert to SCE before wrap so LIANA's SCE branch
# is used; that branch reads from `assay(sce, "logcounts")` and never
# touches the deprecated slot argument. `colLabels(sce)` carries the
# grouping LIANA uses as idents.
.liana_run_one <- function(sub, methods, resource) {
  if (!requireNamespace("liana", quietly = TRUE))
    stop("liana package required; install via remotes::install_github('saezlab/liana').")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("SingleCellExperiment required for the LIANA SCE path.")

  # Make sure the active assay's data layer is populated. LIANA expects
  # log-normalized counts; if a freshly-merged Seurat object still has only
  # raw counts in its data slot, normalize first. Also JoinLayers as a
  # safety net in case the parent merge wasn't consolidated.
  Seurat::DefaultAssay(sub) <- "RNA"
  sub <- tryCatch(SeuratObject::JoinLayers(sub, assay = "RNA"),
                  error = function(e) sub)
  sub <- tryCatch(Seurat::NormalizeData(sub, verbose = FALSE),
                  error = function(e) sub)
  sce <- Seurat::as.SingleCellExperiment(sub, assay = "RNA")
  SingleCellExperiment::colLabels(sce) <-
    as.factor(SummarizedExperiment::colData(sce)$ccc_label)

  res <- liana::liana_wrap(sce,
                           method     = methods,
                           resource   = resource,
                           idents_col = "label",
                           verbose    = FALSE)
  agg <- liana::liana_aggregate(res, verbose = FALSE)
  tibble::as_tibble(agg)
}

# Main entry called from run_pipeline.R Phase 1d.
run_liana_myeloid_tcell <- function(cfg) {
  lcfg <- cfg$liana %||% list()
  methods <- as.character(lcfg$methods %||%
                          c("natmi", "connectome", "logfc",
                            "sca", "cellphonedb"))
  resource <- as.character(lcfg$resource %||% "Consensus")
  conditions <- as.character(lcfg$conditions %||% c("NIU", "Viral"))
  min_cells <- as.integer(lcfg$min_cells_per_cluster %||% 10L)
  source_regex <- as.character(lcfg$source_regex %||% "^myeloid_")
  target_regex <- as.character(lcfg$target_regex %||% "^tcell_")

  paths_myel <- get_target_paths(cfg, "myeloid")
  paths_tcel <- get_target_paths(cfg, "tcell")
  cc_paths   <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/",
    viz    = "outputs/viz/cross_compartment/")
  ensure_dir(cc_paths$tables)

  myel_path <- file.path(paths_myel$results_objects, "IntegratedSeuratObject.rds")
  tcel_path <- file.path(paths_tcel$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(myel_path) || !file.exists(tcel_path)) {
    log_message("liana: compartment objects missing; skipping.")
    return(invisible(FALSE))
  }
  log_message("=== liana_myeloid_tcell ===")
  myel <- readRDS(myel_path)
  tcel <- readRDS(tcel_path)
  myel$ccc_label <- get_substate_key_vector(myel, "myeloid")
  tcel$ccc_label <- get_substate_key_vector(tcel, "tcell")

  merged <- tryCatch(merge(myel, y = tcel,
                           add.cell.ids = c("myeloid", "tcell"),
                           merge.data = TRUE),
                     error = function(e) {
                       log_message("  liana: merge() failed: ", conditionMessage(e))
                       NULL })
  if (is.null(merged)) return(invisible(FALSE))
  rm(myel, tcel); gc()

  # Seurat v5 keeps one layer per source object after merge() — LIANA's
  # `GetAssayData` call chokes on multi-layer v5 assays. Consolidate to a
  # single counts/data layer before anything else touches the RNA assay.
  # NOTE: JoinLayers is exported by SeuratObject (not Seurat) — use
  # SeuratObject::JoinLayers explicitly so the namespace lookup succeeds
  # even when Seurat doesn't re-export it.
  Seurat::DefaultAssay(merged) <- "RNA"
  merged <- tryCatch(SeuratObject::JoinLayers(merged, assay = "RNA"),
                     error = function(e) {
                       log_message("  liana: JoinLayers failed: ",
                                   conditionMessage(e)); merged })

  if (!"Phenotype_2" %in% colnames(merged[[]])) {
    log_message("  liana: Phenotype_2 missing on merged object; aborting.")
    return(invisible(FALSE))
  }
  Seurat::Idents(merged) <- "ccc_label"

  audit_path <- file.path(cc_paths$tables, "liana_cluster_n_audit.csv")
  keep_labels <- .liana_cell_count_audit(merged[[]], "Phenotype_2",
                                         conditions, min_cells, audit_path)
  if (length(keep_labels) < 2L) {
    log_message("  liana: <2 labels pass min_cells=", min_cells,
                "; nothing to test.")
    return(invisible(FALSE))
  }

  per_condition <- list()
  for (cond in conditions) {
    cells_keep <- colnames(merged)[merged$Phenotype_2 == cond &
                                   merged$ccc_label %in% keep_labels]
    if (length(cells_keep) < 200L) {
      log_message("  liana: <200 cells for ", cond, " condition; skipping arm.")
      next
    }
    sub <- subset(merged, cells = cells_keep)
    Seurat::Idents(sub) <- "ccc_label"
    log_message("  liana: running ", cond, " on ",
                length(unique(sub$ccc_label)), " labels x ",
                ncol(sub), " cells")
    df <- tryCatch(.liana_run_one(sub, methods, resource),
                   error = function(e) {
                     log_message("  liana (", cond, ") failed: ",
                                 conditionMessage(e)); NULL })
    if (is.null(df) || nrow(df) == 0L) next

    # LIANA column-name harmonization. Newer versions use ligand.complex /
    # receptor.complex (with periods); some older versions use underscores or
    # plain ligand / receptor. Normalize to ligand_complex / receptor_complex
    # so downstream join + viz read a single canonical name.
    df <- df |>
      dplyr::rename_with(~ "ligand_complex",
                         dplyr::any_of(c("ligand.complex", "ligand_complex"))) |>
      dplyr::rename_with(~ "receptor_complex",
                         dplyr::any_of(c("receptor.complex", "receptor_complex")))
    if (!"ligand_complex" %in% colnames(df) && "ligand" %in% colnames(df)) {
      df$ligand_complex <- df$ligand
    }
    if (!"receptor_complex" %in% colnames(df) && "receptor" %in% colnames(df)) {
      df$receptor_complex <- df$receptor
    }

    df <- df |>
      dplyr::filter(grepl(source_regex, .data$source),
                    grepl(target_regex, .data$target))
    if (!"aggregate_rank" %in% colnames(df)) {
      log_message("  liana: aggregate_rank missing from ", cond,
                  " output; skipping write.")
      next
    }
    if (!all(c("ligand_complex", "receptor_complex") %in% colnames(df))) {
      log_message("  liana: ligand/receptor columns missing after ",
                  "normalization (got: ",
                  paste(colnames(df), collapse = ", "),
                  "); skipping write.")
      next
    }
    # Preserve per-method scores so downstream viz can use a sample-size-
    # invariant logFC-based bias instead of the rank-based one (which is
    # confounded by per-arm cell counts). Key score columns:
    #   logfc.logfc_comb         L+R combined logFC vs other source-target
    #                            pairs in the same condition
    #   natmi.edge_specificity   per-edge specificity (NATMI)
    #   sca.LRscore              specificity weighting (SingleCellSignalR)
    #   cellphonedb.pvalue       permutation p-value
    #   connectome.weight_sc     scaled product (Connectome)
    keep_cols <- c("source", "target", "ligand_complex", "receptor_complex",
                   "aggregate_rank", "mean_rank",
                   "logfc.logfc_comb", "natmi.edge_specificity",
                   "sca.LRscore", "cellphonedb.pvalue",
                   "connectome.weight_sc")
    df <- df |> dplyr::select(dplyr::any_of(keep_cols))
    # Add a ligand-family tag using the existing myeloid_programs config so
    # downstream viz can aggregate / facet by family.
    df$ligand_family <- .liana_classify_ligand_family(df$ligand_complex, cfg)
    per_path <- file.path(cc_paths$tables,
                          paste0("liana_myeloid_to_tcell_",
                                 cond, ".csv"))
    utils::write.csv(df, per_path, row.names = FALSE)
    log_message("  liana[", cond, "]: wrote ", nrow(df), " LR rows to ",
                per_path)
    per_condition[[cond]] <- df
  }
  if (length(per_condition) < 2L) {
    log_message("  liana: only ", length(per_condition),
                " condition succeeded; combined table not built.")
    return(invisible(TRUE))
  }

  # Per-arm top-N CSVs (most-significant pairs in each etiology by aggregate
  # rank). These are the cleanest tables for the manuscript supplement: no
  # rank-difference comparison, no sample-size confound.
  for (cond in names(per_condition)) {
    df <- per_condition[[cond]] |>
      dplyr::arrange(.data$aggregate_rank) |>
      dplyr::slice_head(n = 100L)
    out <- file.path(cc_paths$tables,
                     paste0("liana_top100_", cond, ".csv"))
    utils::write.csv(df, out, row.names = FALSE)
    log_message("  liana top100[", cond, "]: wrote ", out)
  }

  # Combine with outer join. Two disease-bias metrics computed:
  #   disease_bias_rank   = rank_Viral_imp - rank_NIU_imp (positive = NIU-biased,
  #                         negative = Viral-biased). Cohort-size confounded.
  #   disease_bias_logfc  = logfc_NIU - logfc_Viral (positive = NIU-driving in
  #                         absolute expression terms, negative = Viral-driving).
  #                         Sample-size invariant when logfc.logfc_comb is the
  #                         per-cell-type mean log-fold-change against the
  #                         compartment background — the absolute biology number.
  shared_keys <- c("source", "target", "ligand_complex",
                   "receptor_complex", "ligand_family")
  score_cols <- intersect(
    c("aggregate_rank", "mean_rank", "logfc.logfc_comb",
      "natmi.edge_specificity", "sca.LRscore", "cellphonedb.pvalue",
      "connectome.weight_sc"),
    Reduce(intersect, lapply(per_condition, colnames)))

  rename_with_suffix <- function(df, suffix) {
    nms <- colnames(df)
    new_nms <- ifelse(nms %in% shared_keys, nms, paste0(nms, "_", suffix))
    stats::setNames(df, new_nms)
  }
  a <- per_condition[[conditions[1]]] |>
    dplyr::select(dplyr::any_of(c(shared_keys, score_cols))) |>
    rename_with_suffix(conditions[1])
  b <- per_condition[[conditions[2]]] |>
    dplyr::select(dplyr::any_of(c(shared_keys, score_cols))) |>
    rename_with_suffix(conditions[2])
  combined <- dplyr::full_join(a, b, by = shared_keys)

  rank_col_a <- paste0("aggregate_rank_", conditions[1])
  rank_col_b <- paste0("aggregate_rank_", conditions[2])
  logfc_col_a <- paste0("logfc.logfc_comb_", conditions[1])
  logfc_col_b <- paste0("logfc.logfc_comb_", conditions[2])
  max_rank <- max(c(combined[[rank_col_a]], combined[[rank_col_b]]),
                  na.rm = TRUE) + 1

  combined$unique_to <- dplyr::case_when(
    is.na(combined[[rank_col_a]]) & is.na(combined[[rank_col_b]]) ~ NA_character_,
    is.na(combined[[rank_col_a]]) ~ conditions[2],
    is.na(combined[[rank_col_b]]) ~ conditions[1],
    TRUE                          ~ "both")
  combined[[paste0(rank_col_a, "_imp")]] <-
    ifelse(is.na(combined[[rank_col_a]]), max_rank, combined[[rank_col_a]])
  combined[[paste0(rank_col_b, "_imp")]] <-
    ifelse(is.na(combined[[rank_col_b]]), max_rank, combined[[rank_col_b]])
  combined$consensus_rank <-
    pmin(combined[[paste0(rank_col_a, "_imp")]],
         combined[[paste0(rank_col_b, "_imp")]])
  # Rank-based bias (legacy / for backwards-compat, but cohort-size confounded)
  combined$disease_bias_rank <-
    combined[[paste0(rank_col_b, "_imp")]] - combined[[paste0(rank_col_a, "_imp")]]
  # logFC-based bias (primary, sample-size invariant)
  if (logfc_col_a %in% colnames(combined) &&
      logfc_col_b %in% colnames(combined)) {
    a_lf <- combined[[logfc_col_a]]
    b_lf <- combined[[logfc_col_b]]
    # Where a side is missing, impute as zero (no signal) rather than NA so
    # the bias is interpretable for `unique_to` pairs as well.
    a_lf[is.na(a_lf)] <- 0
    b_lf[is.na(b_lf)] <- 0
    combined$disease_bias_logfc <- a_lf - b_lf  # NIU minus Viral
  } else {
    combined$disease_bias_logfc <- NA_real_
  }

  comb_path <- file.path(cc_paths$tables,
                         "liana_myeloid_to_tcell_combined.csv")
  utils::write.csv(combined, comb_path, row.names = FALSE)
  log_message("  liana combined: wrote ", nrow(combined), " LR rows to ",
              comb_path)

  # Family-level summary across all source -> target pairs in each arm.
  # Aggregates the rank of all LR pairs in a family to a per-arm signal.
  for (cond in conditions) {
    rcol <- paste0("aggregate_rank_", cond)
    if (!rcol %in% colnames(combined)) next
    fam_summary <- combined |>
      dplyr::filter(!is.na(.data[[rcol]]),
                    !is.na(.data$ligand_family),
                    .data$ligand_family != "Other") |>
      dplyr::group_by(.data$ligand_family) |>
      dplyr::summarise(
        n_pairs = dplyr::n(),
        mean_rank   = mean(.data[[rcol]], na.rm = TRUE),
        median_rank = stats::median(.data[[rcol]], na.rm = TRUE),
        n_top1pct   = sum(.data[[rcol]] <= 0.01, na.rm = TRUE),
        n_top5pct   = sum(.data[[rcol]] <= 0.05, na.rm = TRUE),
        .groups = "drop")
    fam_summary$condition <- cond
    fam_path <- file.path(cc_paths$tables,
                          paste0("liana_family_summary_", cond, ".csv"))
    utils::write.csv(fam_summary, fam_path, row.names = FALSE)
    log_message("  liana family[", cond, "]: wrote ", fam_path)
  }
  invisible(TRUE)
}
