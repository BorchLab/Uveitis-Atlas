# R/51_eye_blood_repertoire.R
# Per-subject eye<->blood overlap tables for TCR and BCR repertoires. Powers
# Figure 4 panel F (BCR overlap) and Figure 5 panel F/G (TCR overlap +
# alluvial substate-tracking). Operates on the full atlas; reuses the existing
# scRepertoire CTstrict assignments for TCR and the alakazam/dowser clone_id
# from outputs/tables/bcr_airr/*_airr.tsv for BCR (same definition the dowser
# tree-builder uses, more rigorous than scRepertoire defaults for SHM-laden
# BCRs).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Morisita-Horn similarity index between two count vectors (clones x freq).
# 1 = identical distributions, 0 = no overlap.
.morisita <- function(x, y) {
  shared <- intersect(names(x), names(y))
  if (length(shared) == 0) return(0)
  num <- 2 * sum(x[shared] * y[shared])
  den <- (sum(x ^ 2) / sum(x) + sum(y ^ 2) / sum(y)) * sum(x) * sum(y) / (sum(x) * sum(y))
  # Stable form: 2 * sum(xi * yi) / ((sum(xi^2)/X + sum(yi^2)/Y) * X * Y / (X*Y))
  # = 2 * sum(xi * yi) / ((sum(xi^2)/X + sum(yi^2)/Y) * 1)
  X <- sum(x); Y <- sum(y)
  num / ((sum(x ^ 2) / (X * X) + sum(y ^ 2) / (Y * Y)) * X * Y)
}

.jaccard <- function(x, y) {
  a <- names(x); b <- names(y)
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) 0 else inter / uni
}

# Load BCR clone IDs from per-sample AIRR tables. Returns a data.frame with
# cell_id_unique (barcode), subject_id, locus, c_call (isotype), clone_id.
.load_bcr_clones <- function(cfg) {
  bcr_dir <- file.path(cfg$paths$results_tables, "bcr_airr")
  if (!dir.exists(bcr_dir)) {
    log_message("  No BCR AIRR directory at ", bcr_dir, "; skipping BCR pass.")
    return(NULL)
  }
  files <- list.files(bcr_dir, pattern = "_airr\\.tsv$", full.names = TRUE)
  if (length(files) == 0) {
    log_message("  No AIRR tsv files; skipping BCR pass.")
    return(NULL)
  }
  log_message("  Reading ", length(files), " BCR AIRR tables...")
  cols_keep <- c("cell_id_unique", "subject_id", "locus", "c_call", "clone_id",
                 "v_call", "j_call", "junction_aa", "sample_id")
  bcr <- lapply(files, function(f) {
    df <- tryCatch(readr::read_tsv(f, col_types = readr::cols(.default = "c"),
                                   progress = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) return(NULL)
    cols_have <- intersect(cols_keep, colnames(df))
    df[, cols_have, drop = FALSE]
  })
  bcr <- bcr[!sapply(bcr, is.null)]
  if (length(bcr) == 0) return(NULL)
  bcr <- dplyr::bind_rows(bcr) |>
         dplyr::filter(!is.na(clone_id), !is.na(cell_id_unique))
  # Keep one row per (cell, clone) collapsing across heavy/light chains; prefer IGH (heavy)
  bcr <- bcr |>
    dplyr::group_by(cell_id_unique) |>
    dplyr::arrange(dplyr::desc(locus == "IGH")) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()
  bcr
}

# Join BCR clone IDs onto a Seurat-derived metadata frame by barcode. Returns
# the metadata frame augmented with a bcr_clone_id column.
.join_bcr_clones_meta <- function(meta, bcr) {
  meta$barcode_join <- if ("cell_id_unique" %in% colnames(meta)) {
    meta$cell_id_unique
  } else {
    rownames(meta)
  }
  meta <- dplyr::left_join(meta,
                           bcr |> dplyr::select(cell_id_unique, bcr_clone_id = clone_id, c_call),
                           by = c("barcode_join" = "cell_id_unique"))
  meta
}

# ---------------------------------------------------------------------------
# Per-modality computations
# ---------------------------------------------------------------------------

# Restrict meta to cells whose CTgene matches the modality locus prefix.
# Both TCR (scRepertoire) and BCR (R/03 AIRR pipeline) write to a single
# `CTstrict` column in the merged Seurat meta — R/03 places BCR first when
# deduplicating by barcode, so the BCR-derived "<subject>_<clone_id>" string
# wins for B cells. Without a locus check, those BCR clone IDs leak into the
# TCR pass (and would mirror the other way for any TCR cells annotated by
# both). The CTgene prefix (`^TR` vs `^IG`) is the same classifier R/03 uses
# internally; missing CTgene falls through unfiltered so old objects that
# pre-date the column don't get silently emptied.
.locus_filter <- function(meta, modality) {
  re <- switch(modality, "TCR" = "^TR", "BCR" = "^IG", NULL)
  if (is.null(re) || !"CTgene" %in% colnames(meta)) return(meta)
  dplyr::filter(meta, !is.na(.data$CTgene), grepl(re, .data$CTgene))
}

#' Per-subject eye<->blood overlap summary for a given clone column.
.run_overlap <- function(meta, modality, clone_col, cfg) {
  meta <- .locus_filter(meta, modality)
  out <- meta |>
    dplyr::filter(!is.na(.data[[clone_col]]), !is.na(Subject), !is.na(Tissue_1)) |>
    dplyr::mutate(clone = .data[[clone_col]])

  rows <- list()
  for (subj in unique(out$Subject)) {
    s <- out |> dplyr::filter(Subject == subj)
    eye_cells   <- s |> dplyr::filter(Tissue_1 == "Eye")
    blood_cells <- s |> dplyr::filter(Tissue_1 == "Blood")
    if (nrow(eye_cells) == 0 && nrow(blood_cells) == 0) next

    eye_counts   <- table(eye_cells$clone)
    blood_counts <- table(blood_cells$clone)

    n_shared <- length(intersect(names(eye_counts), names(blood_counts)))
    n_eye    <- length(eye_counts)
    n_blood  <- length(blood_counts)

    rows[[length(rows) + 1]] <- data.frame(
      subject     = subj,
      etiology    = unique(s$Etiology)[1],
      phenotype   = unique(s$Phenotype_2)[1],
      n_eye_clones    = n_eye,
      n_blood_clones  = n_blood,
      n_shared        = n_shared,
      n_eye_cells     = nrow(eye_cells),
      n_blood_cells   = nrow(blood_cells),
      frac_eye_shared_with_blood   = if (n_eye > 0)   n_shared / n_eye   else NA_real_,
      frac_blood_shared_with_eye   = if (n_blood > 0) n_shared / n_blood else NA_real_,
      jaccard       = .jaccard(eye_counts, blood_counts),
      morisita      = if (n_eye > 0 && n_blood > 0)
                        .morisita(eye_counts, blood_counts) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(NULL)
  out_df <- dplyr::bind_rows(rows)
  out_path <- file.path("outputs/tables/repertoire",
                        paste0(modality, "_eye_blood_overlap.csv"))
  write.csv(out_df, out_path, row.names = FALSE)
  log_message("  Saved: ", basename(out_path), " (", nrow(out_df), " subjects)")
  out_df
}

#' Top-N expanded eye clones per subject, with their detection in matched blood.
.run_top_expanded <- function(meta, modality, clone_col, cfg, top_n = 20) {
  meta <- .locus_filter(meta, modality)
  out <- meta |>
    dplyr::filter(!is.na(.data[[clone_col]]), !is.na(Subject), !is.na(Tissue_1)) |>
    dplyr::mutate(clone = .data[[clone_col]])

  rows <- list()
  for (subj in unique(out$Subject)) {
    s <- out |> dplyr::filter(Subject == subj)
    et <- unique(s$Etiology)[1]

    eye <- s |> dplyr::filter(Tissue_1 == "Eye") |>
                dplyr::count(clone, name = "n_eye") |>
                dplyr::arrange(dplyr::desc(n_eye)) |>
                dplyr::slice_head(n = top_n)
    if (nrow(eye) == 0) next

    blood <- s |> dplyr::filter(Tissue_1 == "Blood") |>
                  dplyr::count(clone, name = "n_blood")
    eye <- dplyr::left_join(eye, blood, by = "clone") |>
           dplyr::mutate(n_blood = tidyr::replace_na(n_blood, 0L),
                         found_in_blood = n_blood > 0,
                         freq_eye = n_eye / sum(eye$n_eye),
                         freq_blood = if (sum(blood$n_blood) > 0)
                                        n_blood / sum(blood$n_blood) else 0)

    rows[[length(rows) + 1]] <- eye |>
      dplyr::mutate(subject = subj, etiology = et, modality = modality) |>
      dplyr::select(subject, etiology, modality, clone_id = clone,
                    n_cells_eye = n_eye, n_cells_blood = n_blood,
                    found_in_blood, freq_eye, freq_blood)
  }
  if (length(rows) == 0) return(NULL)
  out_df <- dplyr::bind_rows(rows)
  out_path <- file.path("outputs/tables/repertoire",
                        paste0(modality, "_top_expanded_eye.csv"))
  write.csv(out_df, out_path, row.names = FALSE)
  log_message("  Saved: ", basename(out_path), " (", nrow(out_df), " clones)")
  out_df
}

#' Long-format clone x substate x tissue table for the alluvial / Sankey panel.
#' Reads each compartment's IntegratedSeuratObject.rds for substate -> barcode
#' mapping. Out-of-compartment cells get substate = "out_of_compartment".
.run_top_expanded_celltype <- function(meta, top_clones, modality, clone_col, cfg) {
  if (is.null(top_clones) || nrow(top_clones) == 0) return(NULL)
  meta <- .locus_filter(meta, modality)

  # Build a single barcode -> (compartment, substate) map across all compartments.
  cmp_map <- list()
  for (cmp in c("myeloid", "bcell", "tcell")) {
    p <- get_target_paths(cfg, cmp)
    obj_path <- file.path(p$results_objects, "IntegratedSeuratObject.rds")
    if (!file.exists(obj_path)) next
    obj_cmp <- readRDS(obj_path)
    cmp_map[[cmp]] <- data.frame(
      barcode  = colnames(obj_cmp),
      compartment = cmp,
      substate = paste0(cmp, "_", as.character(obj_cmp$knn.leiden.cluster)),
      stringsAsFactors = FALSE
    )
    rm(obj_cmp); invisible(gc(verbose = FALSE))
  }
  if (length(cmp_map) == 0) {
    log_message("  No compartment objects found; skipping per-celltype table.")
    return(NULL)
  }
  cmp_df <- dplyr::bind_rows(cmp_map)

  # Long format: each (clone x cell) becomes (clone x substate x tissue).
  # Capture barcodes from rownames BEFORE filtering, since dplyr drops them.
  meta$barcode_join <- if ("cell_id_unique" %in% colnames(meta)) {
    meta$cell_id_unique
  } else {
    rownames(meta)
  }
  expanded <- meta |>
    dplyr::filter(!is.na(.data[[clone_col]]), !is.na(Subject), !is.na(Tissue_1)) |>
    dplyr::mutate(clone = .data[[clone_col]])

  # Restrict to subjects/clones in top_clones
  keep <- top_clones |> dplyr::distinct(subject, clone_id)
  expanded <- expanded |>
    dplyr::semi_join(keep |>
                     dplyr::rename(Subject = subject, clone = clone_id),
                     by = c("Subject", "clone"))

  expanded <- dplyr::left_join(expanded,
                               cmp_df,
                               by = c("barcode_join" = "barcode"))
  expanded$substate[is.na(expanded$substate)] <- "out_of_compartment"

  out_df <- expanded |>
    dplyr::count(Subject, Etiology, clone, Tissue_1, substate,
                 name = "n_cells") |>
    dplyr::rename(subject = Subject, etiology = Etiology, clone_id = clone,
                  tissue = Tissue_1) |>
    dplyr::mutate(modality = modality)

  out_path <- file.path("outputs/tables/repertoire",
                        paste0(modality, "_top_expanded_eye_celltype.csv"))
  write.csv(out_df, out_path, row.names = FALSE)
  log_message("  Saved: ", basename(out_path), " (", nrow(out_df), " rows)")
  out_df
}

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

run_eye_blood_repertoire <- function(cfg) {
  paths_full <- get_target_paths(cfg, "all")
  obj_path   <- file.path(paths_full$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Full atlas IntegratedSeuratObject.rds not found at ", obj_path,
                ". Skipping eye-blood repertoire.")
    return(invisible(FALSE))
  }
  log_message("=== Phase 2b: eye<->blood repertoire summaries ===")
  ensure_dir("outputs/tables/repertoire")
  obj  <- readRDS(obj_path)
  meta <- obj@meta.data

  top_n <- cfg$repertoire_top_n_eye %||% 20

  # ---- TCR (CTstrict from scRepertoire) ----
  if ("CTstrict" %in% colnames(meta)) {
    log_message("TCR repertoire pass...")
    .run_overlap(meta, "TCR", "CTstrict", cfg)
    top_tcr <- .run_top_expanded(meta, "TCR", "CTstrict", cfg, top_n = top_n)
    .run_top_expanded_celltype(meta, top_tcr, "TCR", "CTstrict", cfg)
  } else {
    log_message("CTstrict column absent; skipping TCR pass.")
  }

  # ---- BCR (alakazam clone_id from AIRR tables) ----
  bcr <- .load_bcr_clones(cfg)
  if (!is.null(bcr)) {
    log_message("BCR repertoire pass (", nrow(bcr), " annotated cells)...")
    meta_b <- .join_bcr_clones_meta(meta, bcr)
    .run_overlap(meta_b, "BCR", "bcr_clone_id", cfg)
    top_bcr <- .run_top_expanded(meta_b, "BCR", "bcr_clone_id", cfg, top_n = top_n)
    .run_top_expanded_celltype(meta_b, top_bcr, "BCR", "bcr_clone_id", cfg)
  }

  log_message("=== Eye<->blood repertoire complete ===")
  invisible(TRUE)
}

# run_eye_blood_repertoire reads the full-atlas IntegratedSeuratObject and the
# per-sample BCR AIRR tables, then writes per-subject overlap and top-expanded
# clone tables for both TCR and BCR. The third table (clone x substate x
# tissue) requires the three compartment Seurat objects from
# 22_subset_compartments.R; if those are absent, the function still produces
# the first two tables and logs a warning. The function operates on the full
# atlas so eye<->blood paired clones can be detected; clone definitions are
# scRepertoire CTstrict for TCR and alakazam/dowser clone_id for BCR (same
# definition used by 52_bcr_lineage.R for tree-building).
