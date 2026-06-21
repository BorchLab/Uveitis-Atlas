# R/56_bcell_lineage_trees.R
# Top-N BCR lineage trees per etiology (NIU, Viral) + per-etiology
# cluster-sharing circos plots. See design doc:
#   docs/plans/2026-05-27-bcell-lineage-trees-circos-design.md
#
# Outputs:
#   outputs/viz/eye/bcell/10_lineage_arch/lineage_trees/bcell_lineage_<etiology>_top<rank>_<subject>_clone<id>.pdf
#   outputs/viz/eye/bcell/08_repertoire/bcell_clonal_sharing_circos_{NIU,Viral}.{pdf,png}
#   outputs/tables/eye/bcell/bcell_top_lineage_summary.csv
#   outputs/tables/eye/bcell/bcell_clonal_sharing_matrix_{NIU,Viral}.csv
#   outputs/objects/ibex/bcell_lineage_top_trees.rds  (tree cache)

suppressPackageStartupMessages({
  library(Seurat)
  library(dowser)
  library(alakazam)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggtree)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Load bcell compartment meta with phenotype + eye-filtered AIRR rows.
# Shared by run_bcell_lineage_trees and run_bcell_clonal_sharing_circos.
#
# Phenotype mapping routes through cfg$etiology_groups$niu/$viral rather
# than the .drop_healthy() helper (R/74) because the bcell compartment
# metadata reliably carries Etiology subtypes, whereas Phenotype_2 may
# carry "Healthy" or NA for some upstream pipeline versions. Keeps R/56
# robust to either convention.
#
# Returns NULL if the bcell object is missing, or a list with:
#   bm        : bcell meta with cell_id_unique + phenotype columns
#   bcr_eye   : AIRR rows for productive IGH heavy chains of eye cells
#               (NULL if AIRR dir is empty)
#   paths     : get_target_paths(cfg, "bcell")
#   paths_all : get_target_paths(cfg, "all")
.bcell_load_eye_airr <- function(cfg) {
  paths     <- get_target_paths(cfg, "bcell")
  paths_all <- get_target_paths(cfg, "all")
  bcell_obj_path <- file.path(paths$results_objects,
                              "IntegratedSeuratObject.rds")
  if (!file.exists(bcell_obj_path)) {
    log_message("  bcell object missing: ", bcell_obj_path)
    return(NULL)
  }
  obj <- readRDS(bcell_obj_path)
  bm <- obj@meta.data
  bm$cell_id_unique <- rownames(bm)
  niu_set   <- as.character(cfg$etiology_groups$niu   %||% character(0))
  viral_set <- as.character(cfg$etiology_groups$viral %||% character(0))
  bm$phenotype <- dplyr::case_when(
    bm$Etiology %in% niu_set   ~ "NIU",
    bm$Etiology %in% viral_set ~ "Viral",
    TRUE                       ~ NA_character_)
  bm <- bm |> dplyr::filter(!is.na(phenotype))
  if (nrow(bm) == 0L) {
    log_message("  no NIU/Viral B cells in compartment.")
    return(list(bm = bm, bcr_eye = NULL,
                paths = paths, paths_all = paths_all))
  }

  airr_dir <- file.path(paths_all$results_tables, "bcr_airr")
  airr_files <- list.files(airr_dir, pattern = "_airr\\.tsv$",
                           full.names = TRUE)
  if (length(airr_files) == 0L) {
    log_message("  no AIRR tables in ", airr_dir)
    return(list(bm = bm, bcr_eye = NULL,
                paths = paths, paths_all = paths_all))
  }
  bcr_db <- dplyr::bind_rows(lapply(airr_files,
                                    utils::read.delim,
                                    stringsAsFactors = FALSE))
  bcr_eye <- bcr_db |>
    dplyr::filter(productive == TRUE, locus == "IGH",
                  !is.na(clone_id),
                  cell_id_unique %in% bm$cell_id_unique)
  list(bm = bm, bcr_eye = bcr_eye,
       paths = paths, paths_all = paths_all)
}

# Sum heavy-chain CDR + FWR R+S SHM frequency. Mirrors .shm_total_col in
# R/86_viz_bcell.R but reads from an arbitrary df, not obj@meta.data.
.bcell_lineage_shm_total <- function(df) {
  cols <- intersect(c("mu_freq_cdr_r_heavy", "mu_freq_cdr_s_heavy",
                      "mu_freq_fwr_r_heavy", "mu_freq_fwr_s_heavy"),
                    colnames(df))
  if (length(cols) == 0L) return(rep(NA_real_, nrow(df)))
  rowSums(as.matrix(df[, cols, drop = FALSE]), na.rm = TRUE)
}

# Build dowser PML trees for the requested (subject, clone_id) candidates.
# bcr_eye must be pre-filtered to productive IGH rows for eye cells only.
# Uses formatClones() before getTrees() (required by dowser >=2.0). Caches
# to disk; subsequent calls re-use the cache and only build missing clones.
.bcell_lineage_build_trees <- function(candidates, bcr_eye, cache_path) {
  cached <- if (file.exists(cache_path)) {
    log_message("  reading tree cache: ", cache_path)
    readRDS(cache_path)
  } else list()

  for (i in seq_len(nrow(candidates))) {
    subject <- candidates$subject[i]
    clone   <- as.character(candidates$clone_id[i])
    key     <- paste0(subject, "|", clone)
    if (!is.null(cached[[key]])) {
      log_message("  [cache hit] ", key)
      next
    }
    log_message("  [build] ", key, " (",
                candidates$phenotype[i], ", n_eye=",
                candidates$n_cells_eye[i], ")")
    clone_seqs <- bcr_eye |>
      dplyr::filter(subject_id == subject,
                    as.character(clone_id) == clone)
    if (nrow(clone_seqs) < 3L) {
      log_message("    skip: only ", nrow(clone_seqs),
                  " heavy-chain seqs")
      next
    }
    if (!"germline_alignment_d_mask" %in% colnames(clone_seqs) ||
        all(is.na(clone_seqs$germline_alignment_d_mask))) {
      log_message("    skip: no germline_alignment_d_mask column",
                  " (re-run R/03 to add)")
      next
    }
    fc <- tryCatch(
      dowser::formatClones(clone_seqs,
                           traits = c("cell_id_unique", "c_call"),
                           text_fields = "c_call"),
      error = function(e) {
        log_message("    formatClones failed: ", conditionMessage(e))
        NULL })
    if (is.null(fc) || nrow(fc) == 0L) next
    trees <- tryCatch(
      dowser::getTrees(fc, build = "pml", nproc = 1),
      error = function(e) {
        log_message("    getTrees failed: ", conditionMessage(e))
        NULL })
    if (!is.null(trees) && nrow(trees) > 0L) {
      cached[[key]] <- trees
      log_message("    built (",
                  length(trees$trees[[1]]$tip.label), " tips)")
    }
  }

  ensure_dir(dirname(cache_path))
  saveRDS(cached, cache_path)
  log_message("  wrote tree cache (", length(cached), " trees): ",
              cache_path)
  cached
}

# Render one dowser tree object as a ggtree PDF. tree_row is one row of
# the tibble returned by dowser::getTrees(): trees[[1]] is the ape::phylo,
# data[[1]] is an airrClone whose @data carries (sequence_id, cell_id_unique).
# bcell_meta has substate + isotype + SHM joined. Returns summary stats.
.bcell_lineage_plot_one <- function(tree_row, candidate_row, bcell_meta,
                                    out_path, w = 7, h = 7) {
  ph <- tree_row$trees[[1]]  # ape::phylo
  if (!inherits(ph, "phylo")) return(NULL)

  # Map tip labels (sequence_id) -> cell_id_unique via airrClone@data, then
  # join cell_id_unique -> bcell_meta for substate / SHM. Isotype comes
  # straight from the AIRR c_call (bcell meta lacks c_call_heavy). The
  # germline root carries NAs and renders in the na.value color/shape.
  airr <- tree_row$data[[1]]@data
  tip_df <- data.frame(label = ph$tip.label, stringsAsFactors = FALSE)
  airr_idx <- match(tip_df$label, airr$sequence_id)
  tip_df$cell_id_unique <- airr$cell_id_unique[airr_idx]
  iso_raw <- if ("c_call" %in% colnames(airr)) airr$c_call[airr_idx]
             else rep(NA_character_, length(airr_idx))
  iso_clean <- sub("\\*.*$", "", as.character(iso_raw))
  tip_df$isotype <- dplyr::case_when(
    grepl("^IGHM", iso_clean) ~ "IGHM",
    grepl("^IGHD", iso_clean) ~ "IGHD",
    grepl("^IGHG", iso_clean) ~ "IGHG",
    grepl("^IGHA", iso_clean) ~ "IGHA",
    grepl("^IGHE", iso_clean) ~ "IGHE",
    is.na(iso_clean)          ~ NA_character_,
    TRUE                      ~ "Other")
  meta_lookup <- bcell_meta[match(tip_df$cell_id_unique,
                                  bcell_meta$cell_id_unique), ]
  tip_df$substate <- meta_lookup$substate
  tip_df$shm      <- meta_lookup$SHM_total
  tip_df$tissue   <- meta_lookup$Tissue_1

  mean_shm <- mean(tip_df$shm, na.rm = TRUE)
  max_shm  <- max(tip_df$shm,  na.rm = TRUE)
  tree_depth <- tryCatch(max(ape::node.depth.edgelength(ph)),
                         error = function(e) NA_real_)

  title <- sprintf("%s | clone %s | n_eye=%d | mean SHM=%.3f",
                   candidate_row$subject, candidate_row$clone_id,
                   candidate_row$n_cells_eye, mean_shm)

  p <- ggtree::ggtree(ph) %<+% tip_df +
    ggtree::geom_tippoint(aes(color = substate, shape = isotype),
                          size = 4, na.rm = FALSE) +
    scale_color_viridis_d(option = "viridis", na.value = "grey80",
                          name = "Substate") +
    scale_shape_manual(values = c(IGHM = 16, IGHD = 17, IGHG = 15,
                                  IGHA = 18, IGHE = 8, Other = 4),
                       na.value = 4, name = "Isotype") +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold", size = 10),
          legend.position = "right")

  ensure_dir(dirname(out_path))
  ggsave(out_path, p, width = w, height = h)
  log_message("  saved: ", basename(out_path))
  list(
    n_seqs_in_tree = length(ph$tip.label),
    mean_shm = mean_shm,
    max_shm = max_shm,
    tree_depth = tree_depth,
    n_substates_in_tree = length(unique(stats::na.omit(tip_df$substate))),
    isotype_set = paste(sort(unique(stats::na.omit(tip_df$isotype))),
                        collapse = ",")
  )
}

# Build a per-cell lookup with substate, collapsed isotype, SHM, tissue.
.bcell_lineage_meta <- function(cfg) {
  paths <- get_target_paths(cfg, "bcell")
  obj <- readRDS(file.path(paths$results_objects,
                           "IntegratedSeuratObject.rds"))
  m <- obj@meta.data
  m$cell_id_unique <- m$cell_id_unique %||% rownames(m)
  m$substate <- substate_labels(cfg, "bcell", m$knn.leiden.cluster)
  iso_col <- intersect(c("c_call_heavy", "c_call", "isotype"),
                       colnames(m))[1]
  m$isotype_collapsed <- if (!is.na(iso_col)) {
    cc <- sub("\\*.*$", "", as.character(m[[iso_col]]))
    dplyr::case_when(
      grepl("^IGHM", cc) ~ "IGHM",
      grepl("^IGHD", cc) ~ "IGHD",
      grepl("^IGHG", cc) ~ "IGHG",
      grepl("^IGHA", cc) ~ "IGHA",
      grepl("^IGHE", cc) ~ "IGHE",
      TRUE              ~ "Other")
  } else NA_character_
  m$SHM_total <- .bcell_lineage_shm_total(m)
  m[, c("cell_id_unique", "substate", "isotype_collapsed",
        "SHM_total", "Tissue_1")]
}

run_bcell_lineage_trees <- function(cfg, top_n = 10L) {
  log_message("=== bcell lineage trees (top ", top_n,
              " per etiology) ===")
  loaded <- .bcell_load_eye_airr(cfg)
  if (is.null(loaded) || is.null(loaded$bcr_eye) ||
      nrow(loaded$bcr_eye) == 0L) {
    log_message("  load failed or no eye IGH AIRR rows; skipping.")
    return(invisible(FALSE))
  }
  bm        <- loaded$bm
  bcr_eye   <- loaded$bcr_eye
  paths     <- loaded$paths
  paths_all <- loaded$paths_all
  log_message("  eye IGH AIRR rows: ", nrow(bcr_eye))

  # ---- Per-clone eye-cell counts; assign phenotype via subject ----
  subj_phen <- bm |>
    dplyr::distinct(Subject, phenotype) |>
    dplyr::rename(subject_id = Subject)
  clone_counts <- bcr_eye |>
    dplyr::distinct(subject_id, clone_id, cell_id_unique) |>
    dplyr::count(subject_id, clone_id, name = "n_cells_eye") |>
    dplyr::left_join(subj_phen, by = "subject_id") |>
    dplyr::filter(!is.na(phenotype))

  # ---- Top top_n + buffer per etiology, by n_cells_eye ----
  # Buffer is generous because dowser's formatClones/getTrees often drops
  # small clones whose sequences are identical (especially NIU clones,
  # which carry less SHM diversity than viral GC-driven clones).
  buffer <- 25L
  candidates <- clone_counts |>
    dplyr::group_by(phenotype) |>
    dplyr::arrange(dplyr::desc(n_cells_eye), .by_group = TRUE) |>
    dplyr::slice_head(n = top_n + buffer) |>
    dplyr::mutate(rank_within_etiology = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::rename(subject = subject_id)
  log_message("  candidate clones: NIU=",
              sum(candidates$phenotype == "NIU"),
              " Viral=", sum(candidates$phenotype == "Viral"))

  # ---- Build trees (cached) ----
  cache_path <- file.path(paths_all$results_objects,
                          "ibex/bcell_lineage_top_trees.rds")
  trees_cached <- .bcell_lineage_build_trees(candidates, bcr_eye,
                                             cache_path)
  log_message("  trees available: ", length(trees_cached))

  # ---- Trim candidates to those with a built tree; take top_n per etiology
  candidates$key <- paste0(candidates$subject, "|", candidates$clone_id)
  candidates$tree_built <- candidates$key %in% names(trees_cached)
  selected <- candidates |>
    dplyr::filter(tree_built) |>
    dplyr::group_by(phenotype) |>
    dplyr::arrange(dplyr::desc(n_cells_eye), .by_group = TRUE) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()
  log_message("  selected for rendering: NIU=",
              sum(selected$phenotype == "NIU"),
              " Viral=", sum(selected$phenotype == "Viral"))

  # ---- Render one PDF per selected clone ----
  bcell_meta <- .bcell_lineage_meta(cfg)
  out_dir <- file.path(viz_subdir(paths, "lineage_arch"), "lineage_trees")
  ensure_dir(out_dir)

  summary_rows <- list()
  for (i in seq_len(nrow(selected))) {
    cand <- selected[i, , drop = FALSE]
    tree_row <- trees_cached[[cand$key]]
    if (is.null(tree_row)) next
    fname <- sprintf("bcell_lineage_%s_top%02d_%s_clone%s.pdf",
                     cand$phenotype, cand$rank_within_etiology,
                     gsub("[^A-Za-z0-9]", "_", cand$subject),
                     gsub("[^A-Za-z0-9]", "_", cand$clone_id))
    stats <- tryCatch(
      .bcell_lineage_plot_one(tree_row, cand, bcell_meta,
                              file.path(out_dir, fname)),
      error = function(e) {
        log_message("    render failed for ", cand$key, ": ",
                    conditionMessage(e)); NULL })
    if (!is.null(stats)) {
      summary_rows[[i]] <- cbind(
        cand[, c("subject", "clone_id", "phenotype",
                 "rank_within_etiology", "n_cells_eye")],
        as.data.frame(stats, stringsAsFactors = FALSE))
    }
  }

  out_table <- file.path("outputs/tables/eye/bcell",
                         "bcell_top_lineage_summary.csv")
  ensure_dir(dirname(out_table))
  if (length(summary_rows) > 0L) {
    summ <- dplyr::bind_rows(summary_rows)
    utils::write.csv(summ, out_table, row.names = FALSE)
    log_message("  wrote summary table (", nrow(summ), " rows): ",
                out_table)
  } else {
    log_message("  no trees rendered; summary table not written.")
  }
  invisible(TRUE)
}

.bcell_circos_save <- function(M, etiology, out_base, w = 6, h = 6) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    log_message("  circlize not installed; skipping ", etiology, " chord.")
    return(invisible())
  }
  if (sum(M) == 0L) {
    log_message("  ", etiology, ": no shared clones; skipping chord plot.")
    return(invisible())
  }
  substates <- rownames(M)
  pal <- viridis::viridis(length(substates), option = "viridis", end = 0.9)
  names(pal) <- substates
  .draw <- function() {
    circlize::circos.clear()
    circlize::circos.par(start.degree = 90, gap.degree = 4)
    circlize::chordDiagram(
      M, grid.col = pal, transparency = 0.4,
      annotationTrack = c("name", "grid"),
      annotationTrackHeight = c(0.05, 0.05))
    title(sprintf("%s: BCR clonal sharing across B/plasma substates",
                  etiology))
    circlize::circos.clear()
  }
  ensure_dir(dirname(out_base))
  grDevices::pdf(paste0(out_base, ".pdf"), width = w, height = h); .draw()
  grDevices::dev.off()
  grDevices::png(paste0(out_base, ".png"), width = w * 100,
                 height = h * 100, res = 100); .draw()
  grDevices::dev.off()
  log_message("  saved: ", basename(out_base), ".{pdf,png}")
}

run_bcell_clonal_sharing_circos <- function(cfg) {
  log_message("=== bcell clonal sharing circos ===")
  loaded <- .bcell_load_eye_airr(cfg)
  if (is.null(loaded) || is.null(loaded$bcr_eye) ||
      nrow(loaded$bcr_eye) == 0L) {
    log_message("  load failed or no eye IGH AIRR rows; skipping.")
    return(invisible(FALSE))
  }
  bm      <- loaded$bm
  bcr_eye <- loaded$bcr_eye
  paths   <- loaded$paths
  bm$substate <- substate_labels(cfg, "bcell", bm$knn.leiden.cluster)
  bm <- bm |> dplyr::filter(!is.na(substate))
  if (nrow(bm) == 0L) {
    log_message("  no cells with substate; skipping.")
    return(invisible(FALSE))
  }

  # Build per-cell (subject, clone_id, cell_id_unique, substate, phenotype)
  cell_clone <- bcr_eye |>
    dplyr::distinct(subject_id, clone_id, cell_id_unique) |>
    dplyr::left_join(bm[, c("cell_id_unique", "substate", "phenotype")],
                     by = "cell_id_unique") |>
    dplyr::filter(!is.na(substate), !is.na(phenotype)) |>
    dplyr::mutate(clone_key = paste0(subject_id, "|", clone_id))

  out_tab_dir <- "outputs/tables/eye/bcell"
  ensure_dir(out_tab_dir)
  matrices <- list()
  for (eti in c("NIU", "Viral")) {
    df <- cell_clone |> dplyr::filter(phenotype == eti)
    if (nrow(df) == 0L) {
      log_message("  ", eti, ": no cells; skipping.")
      next
    }
    clone_substates <- df |>
      dplyr::distinct(clone_key, substate)
    substates <- sort(unique(clone_substates$substate))
    K <- length(substates)
    M <- matrix(0L, nrow = K, ncol = K,
                dimnames = list(substates, substates))
    for (cid in unique(clone_substates$clone_key)) {
      s <- clone_substates$substate[clone_substates$clone_key == cid]
      if (length(s) < 2L) next
      pairs <- utils::combn(s, 2)
      for (j in seq_len(ncol(pairs))) {
        a <- pairs[1, j]; b <- pairs[2, j]
        M[a, b] <- M[a, b] + 1L
        M[b, a] <- M[b, a] + 1L
      }
    }
    long <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
    colnames(long) <- c("from", "to", "n_shared_clones")
    long <- long[long$from != long$to, ]
    out_csv <- file.path(out_tab_dir,
                         paste0("bcell_clonal_sharing_matrix_",
                                eti, ".csv"))
    utils::write.csv(long, out_csv, row.names = FALSE)
    log_message("  ", eti, ": ", sum(long$n_shared_clones) / 2L,
                " inter-substate clone pairs (", nrow(M), " substates);",
                " wrote ", out_csv)
    matrices[[eti]] <- M

    out_base <- file.path(viz_subdir(paths, "repertoire"),
                          paste0("bcell_clonal_sharing_circos_", eti))
    .bcell_circos_save(M, eti, out_base)
  }
  invisible(matrices)
}
