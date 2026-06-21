# R/49_liana_bcell.R
# B cell-centered cell-cell communication via LIANA consensus rank.
# Mirrors R/47_liana_myeloid_tcell.R but iterates four directions:
#   myeloid -> bcell    (BAFF/APRIL, CXCL13, IL-6, etc.)
#   tcell   -> bcell    (Tfh-like help: CD40L, IL-21, ICOS-ICOSLG, CXCL13)
#   bcell   -> tcell    (MHC-II / costim / checkpoint feedback)
#   bcell   -> myeloid  (B as APC for myeloid sensing)
#
# Reuses .liana_run_one() and .liana_cell_count_audit() from R/47 (loaded by
# run_pipeline.R before this file). The ligand-family annotation is
# generalized so it can use myeloid_programs OR bcell_programs OR
# tcell_programs depending on which compartment is the SOURCE of the LR pair
# (which is where ligands live).
#
# Outputs under outputs/tables/cross_compartment/:
#   liana_<source>_to_<target>_NIU.csv
#   liana_<source>_to_<target>_Viral.csv
#   liana_<source>_to_<target>_combined.csv
#   liana_bcell_cluster_n_audit.csv
#   liana_bcell_family_summary_{NIU,Viral}.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Generalized ligand-family classifier. Same logic as
# .liana_classify_ligand_family() in R/47, but the programs list is supplied
# directly rather than pulled from cfg$myeloid_programs. Lets one call mix
# multiple compartment program panels (e.g., for direction == bcell_to_tcell
# the source is B, so we use bcell_programs as the ligand-side panel).
.liana_classify_with_programs <- function(ligands, programs_list) {
  if (length(programs_list) == 0L)
    return(rep(NA_character_, length(ligands)))
  fam_lookup <- character(0)
  for (fam in names(programs_list)) {
    g <- as.character(programs_list[[fam]])
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

# Pull the program panel that lives on the source compartment. Falls back
# to an empty list when the compartment's program block isn't configured.
.programs_for_compartment <- function(cfg, compartment) {
  switch(compartment,
         myeloid = cfg$myeloid_programs %||% list(),
         tcell   = cfg$tcell_programs   %||% list(),
         bcell   = cfg$bcell_programs   %||% list(),
         list())
}

# Build the merged Seurat object for one direction. Reads each source/target
# compartment object, stamps ccc_label = "<compartment>_<cluster_id>", merges
# them with compartment-prefixed cell IDs, joins layers, and returns the
# merged object plus the source/target ccc_label set so the LIANA filter can
# enforce direction.
.bcell_liana_build <- function(cfg, dir_cfg) {
  src_cmps <- dir_cfg$source_compartments
  tgt_cmps <- dir_cfg$target_compartments
  all_cmps <- unique(c(src_cmps, tgt_cmps))

  parts <- list()
  cell_id_tag <- character(0)
  for (cmp in all_cmps) {
    p <- get_target_paths(cfg, cmp)
    obj_path <- file.path(p$results_objects, "IntegratedSeuratObject.rds")
    if (!file.exists(obj_path)) {
      log_message("  liana_bcell[", dir_cfg$name, "]: ", cmp,
                  " object missing at ", obj_path, "; aborting direction.")
      return(NULL)
    }
    o <- readRDS(obj_path)
    o$ccc_label <- get_substate_key_vector(o, cmp)
    parts[[cmp]] <- o
    cell_id_tag <- c(cell_id_tag, cmp)
  }
  if (length(parts) == 1L) {
    merged <- parts[[1]]
  } else {
    first <- parts[[1]]
    rest  <- parts[-1]
    merged <- tryCatch(merge(first, y = rest,
                             add.cell.ids = cell_id_tag,
                             merge.data   = TRUE),
                       error = function(e) {
                         log_message("  liana_bcell[", dir_cfg$name,
                                     "]: merge failed: ",
                                     conditionMessage(e)); NULL })
    if (is.null(merged)) return(NULL)
  }
  rm(parts); gc(verbose = FALSE)

  Seurat::DefaultAssay(merged) <- "RNA"
  merged <- tryCatch(SeuratObject::JoinLayers(merged, assay = "RNA"),
                     error = function(e) {
                       log_message("  liana_bcell[", dir_cfg$name,
                                   "]: JoinLayers failed: ",
                                   conditionMessage(e)); merged })
  Seurat::Idents(merged) <- "ccc_label"
  merged
}

# Run one direction (source compartments -> target compartments) across both
# conditions. Returns a list keyed by condition with the per-condition LR
# table; also writes per-condition CSVs and combined+family summaries.
.bcell_liana_one_direction <- function(cfg, dir_cfg, cc_paths) {
  log_message("--- liana_bcell direction: ", dir_cfg$name, " ---")
  merged <- .bcell_liana_build(cfg, dir_cfg)
  if (is.null(merged)) return(invisible(NULL))

  lcfg     <- cfg$liana_bcell %||% list()
  methods  <- as.character(lcfg$methods    %||%
                           c("natmi","connectome","logfc","sca","cellphonedb"))
  resource <- as.character(lcfg$resource   %||% "Consensus")
  conds    <- as.character(lcfg$conditions %||% c("NIU","Viral"))
  min_cells <- as.integer(lcfg$min_cells_per_cluster %||% 10L)
  src_regex <- as.character(dir_cfg$source_regex)
  tgt_regex <- as.character(dir_cfg$target_regex)

  if (!"Phenotype_2" %in% colnames(merged[[]])) {
    log_message("  liana_bcell[", dir_cfg$name, "]: Phenotype_2 missing; ",
                "aborting direction.")
    return(invisible(NULL))
  }

  # Audit cell counts (per direction); writes once with a direction-suffixed name
  audit_path <- file.path(cc_paths$tables,
                          paste0("liana_bcell_cluster_n_audit_",
                                 dir_cfg$name, ".csv"))
  keep_labels <- .liana_cell_count_audit(merged[[]], "Phenotype_2",
                                         conds, min_cells, audit_path)
  if (length(keep_labels) < 2L) {
    log_message("  liana_bcell[", dir_cfg$name, "]: <2 labels pass min_cells=",
                min_cells, "; skipping direction.")
    return(invisible(NULL))
  }

  # Programs for ligand-family annotation: use the SOURCE compartment's panel.
  # Multi-source directions union the panels (not currently used, but the API
  # supports it).
  programs_merged <- do.call(c, lapply(dir_cfg$source_compartments,
                                       function(c) .programs_for_compartment(cfg, c)))

  per_condition <- list()
  for (cond in conds) {
    cells_keep <- colnames(merged)[merged$Phenotype_2 == cond &
                                   merged$ccc_label %in% keep_labels]
    if (length(cells_keep) < 200L) {
      log_message("  liana_bcell[", dir_cfg$name, "][", cond,
                  "]: <200 cells; skipping arm.")
      next
    }
    sub <- subset(merged, cells = cells_keep)
    Seurat::Idents(sub) <- "ccc_label"
    log_message("  liana_bcell[", dir_cfg$name, "][", cond,
                "]: running on ", length(unique(sub$ccc_label)),
                " labels x ", ncol(sub), " cells")

    df <- tryCatch(.liana_run_one(sub, methods, resource),
                   error = function(e) {
                     log_message("  liana_bcell[", dir_cfg$name, "][",
                                 cond, "] failed: ", conditionMessage(e))
                     NULL })
    if (is.null(df) || nrow(df) == 0L) next

    # Harmonize LIANA column names (same logic as R/47)
    df <- df |>
      dplyr::rename_with(~ "ligand_complex",
                         dplyr::any_of(c("ligand.complex", "ligand_complex"))) |>
      dplyr::rename_with(~ "receptor_complex",
                         dplyr::any_of(c("receptor.complex", "receptor_complex")))
    if (!"ligand_complex" %in% colnames(df) && "ligand" %in% colnames(df))
      df$ligand_complex <- df$ligand
    if (!"receptor_complex" %in% colnames(df) && "receptor" %in% colnames(df))
      df$receptor_complex <- df$receptor

    df <- df |>
      dplyr::filter(grepl(src_regex, .data$source),
                    grepl(tgt_regex, .data$target))
    if (!all(c("aggregate_rank","ligand_complex","receptor_complex") %in%
             colnames(df))) {
      log_message("  liana_bcell[", dir_cfg$name, "][", cond,
                  "]: required cols missing; skip write.")
      next
    }
    keep_cols <- c("source","target","ligand_complex","receptor_complex",
                   "aggregate_rank","mean_rank",
                   "logfc.logfc_comb","natmi.edge_specificity",
                   "sca.LRscore","cellphonedb.pvalue","connectome.weight_sc")
    df <- df |> dplyr::select(dplyr::any_of(keep_cols))
    df$ligand_family <- .liana_classify_with_programs(df$ligand_complex,
                                                     programs_merged)
    df$direction <- dir_cfg$name

    per_path <- file.path(cc_paths$tables,
                          paste0("liana_", dir_cfg$name, "_", cond, ".csv"))
    utils::write.csv(df, per_path, row.names = FALSE)
    log_message("  liana_bcell[", dir_cfg$name, "][", cond, "]: wrote ",
                nrow(df), " LR rows to ", per_path)
    per_condition[[cond]] <- df
  }
  if (length(per_condition) == 0L) return(invisible(NULL))

  # Per-condition top-N
  top_n <- as.integer(lcfg$top_n_for_dotplot %||% 25L)
  for (cond in names(per_condition)) {
    out <- per_condition[[cond]] |>
      dplyr::arrange(.data$aggregate_rank) |>
      dplyr::slice_head(n = top_n)
    utils::write.csv(out,
      file.path(cc_paths$tables,
                paste0("liana_", dir_cfg$name, "_top", top_n,
                       "_", cond, ".csv")),
      row.names = FALSE)
  }

  if (length(per_condition) == 2L) {
    # Combined with logFC-bias (mirrors R/47 lines 281+)
    shared_keys <- c("source","target","ligand_complex","receptor_complex",
                     "ligand_family","direction")
    rename_with_suffix <- function(df, sfx) {
      nms <- colnames(df)
      stats::setNames(df,
                      ifelse(nms %in% shared_keys, nms, paste0(nms, "_", sfx)))
    }
    score_cols <- intersect(
      c("aggregate_rank","mean_rank","logfc.logfc_comb",
        "natmi.edge_specificity","sca.LRscore","cellphonedb.pvalue",
        "connectome.weight_sc"),
      Reduce(intersect, lapply(per_condition, colnames)))

    a <- per_condition[[conds[1]]] |>
      dplyr::select(dplyr::any_of(c(shared_keys, score_cols))) |>
      rename_with_suffix(conds[1])
    b <- per_condition[[conds[2]]] |>
      dplyr::select(dplyr::any_of(c(shared_keys, score_cols))) |>
      rename_with_suffix(conds[2])
    combined <- dplyr::full_join(a, b, by = shared_keys)

    rcol_a <- paste0("aggregate_rank_", conds[1])
    rcol_b <- paste0("aggregate_rank_", conds[2])
    if (rcol_a %in% colnames(combined) && rcol_b %in% colnames(combined)) {
      max_rank <- max(c(combined[[rcol_a]], combined[[rcol_b]]),
                      na.rm = TRUE) + 1
      combined$unique_to <- dplyr::case_when(
        is.na(combined[[rcol_a]]) & is.na(combined[[rcol_b]]) ~ NA_character_,
        is.na(combined[[rcol_a]]) ~ conds[2],
        is.na(combined[[rcol_b]]) ~ conds[1],
        TRUE                      ~ "both")
      combined[[paste0(rcol_a, "_imp")]] <-
        ifelse(is.na(combined[[rcol_a]]), max_rank, combined[[rcol_a]])
      combined[[paste0(rcol_b, "_imp")]] <-
        ifelse(is.na(combined[[rcol_b]]), max_rank, combined[[rcol_b]])
      combined$consensus_rank <-
        pmin(combined[[paste0(rcol_a, "_imp")]],
             combined[[paste0(rcol_b, "_imp")]])
      lcol_a <- paste0("logfc.logfc_comb_", conds[1])
      lcol_b <- paste0("logfc.logfc_comb_", conds[2])
      if (lcol_a %in% colnames(combined) && lcol_b %in% colnames(combined)) {
        a_lf <- combined[[lcol_a]]; b_lf <- combined[[lcol_b]]
        a_lf[is.na(a_lf)] <- 0; b_lf[is.na(b_lf)] <- 0
        combined$disease_bias_logfc <- a_lf - b_lf   # cond1 minus cond2
      }
    }
    comb_path <- file.path(cc_paths$tables,
                           paste0("liana_", dir_cfg$name, "_combined.csv"))
    utils::write.csv(combined, comb_path, row.names = FALSE)
    log_message("  liana_bcell[", dir_cfg$name, "]: wrote combined to ",
                comb_path)
  }
  invisible(per_condition)
}

run_liana_bcell <- function(cfg) {
  lcfg <- cfg$liana_bcell %||% list()
  directions <- lcfg$directions %||% list()
  if (length(directions) == 0L) {
    log_message("liana_bcell: no directions configured; skipping.")
    return(invisible(FALSE))
  }
  cc_paths <- cfg$paths_cross_compartment %||% list(
    tables = "outputs/tables/cross_compartment/",
    viz    = "outputs/viz/cross_compartment/")
  ensure_dir(cc_paths$tables)

  log_message("=== liana_bcell: ", length(directions), " directions ===")

  all_per <- list()
  for (d in directions) {
    res <- tryCatch(.bcell_liana_one_direction(cfg, d, cc_paths),
                    error = function(e) {
                      log_message("  liana_bcell[", d$name, "]: ",
                                  conditionMessage(e))
                      NULL })
    if (!is.null(res)) {
      for (cond in names(res)) {
        all_per[[paste0(d$name, "_", cond)]] <- res[[cond]]
      }
    }
  }
  if (length(all_per) == 0L) {
    log_message("liana_bcell: no successful directions.")
    return(invisible(FALSE))
  }

  # Family-level summary per condition, pooled across directions. One row per
  # (direction, condition, ligand_family).
  conds <- as.character(lcfg$conditions %||% c("NIU","Viral"))
  for (cond in conds) {
    fam_rows <- list()
    for (key in names(all_per)) {
      if (!grepl(paste0("_", cond, "$"), key)) next
      df <- all_per[[key]]
      if (is.null(df) || nrow(df) == 0L) next
      dir_name <- sub(paste0("_", cond, "$"), "", key)
      tmp <- df |>
        dplyr::filter(!is.na(ligand_family), ligand_family != "Other") |>
        dplyr::group_by(ligand_family) |>
        dplyr::summarise(
          n_pairs      = dplyr::n(),
          mean_rank    = mean(aggregate_rank, na.rm = TRUE),
          median_rank  = stats::median(aggregate_rank, na.rm = TRUE),
          n_top1pct    = sum(aggregate_rank <= 0.01, na.rm = TRUE),
          n_top5pct    = sum(aggregate_rank <= 0.05, na.rm = TRUE),
          .groups = "drop")
      tmp$direction <- dir_name
      tmp$condition <- cond
      fam_rows[[key]] <- tmp
    }
    if (length(fam_rows) > 0L) {
      fam_summary <- dplyr::bind_rows(fam_rows)
      fam_path <- file.path(cc_paths$tables,
                            paste0("liana_bcell_family_summary_",
                                   cond, ".csv"))
      utils::write.csv(fam_summary, fam_path, row.names = FALSE)
      log_message("  liana_bcell family[", cond, "]: wrote ", fam_path)
    }
  }
  invisible(TRUE)
}

# run_liana_bcell iterates four cross-compartment directions
# (myeloid->bcell, tcell->bcell, bcell->tcell, bcell->myeloid) and for each
# merges the source + target compartment Seurat objects (each carrying its
# substate_key as ccc_label), filters per Phenotype_2 (NIU vs Viral), and
# calls the shared .liana_run_one() wrapper from R/47 to produce a rank-
# aggregated LR table restricted by source_regex/target_regex. Each direction
# writes per-condition CSVs, an NIU vs Viral combined table with logFC-based
# disease bias, and a family summary using the source compartment's program
# panel for ligand-family annotation. Outputs land under
# outputs/tables/cross_compartment/liana_<source>_to_<target>_*.csv. Reuses
# .liana_run_one(), .liana_cell_count_audit(), and the column-harmonization
# logic from R/47_liana_myeloid_tcell.R.
