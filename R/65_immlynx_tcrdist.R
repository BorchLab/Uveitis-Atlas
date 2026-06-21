# R/65_immlynx_tcrdist.R
# tcrdist + ESM-2 embedding via immLynx (BorchLab/immLynx).
# Builds a SingleCellExperiment from prepare_intraocular_trb() (R/64) and runs:
#   - runEmbeddings()      -> ESM-2 protein language model embedding
#   - scater::runUMAP()    -> 2-D viz coordinates on the embedding
#   - runTCRdist()         -> pairwise tcrdist3 distance matrix (list output)
#   - runClustTCR()        -> CDR3 motif clusters (MCL)
# Metaclonotypist is intentionally skipped because it wraps tcrdist3 on the
# same data with extra edit-distance prefilter; runTCRdist + runClustTCR
# already cover that ground.
#
# Each immLynx step is wrapped in tryCatch so a downstream failure does not
# prevent partial results from being persisted to ImmLynxTcrdistResults.rds.

`%||%` <- function(x, y) if (is.null(x)) y else x

# tcrdist3 deduplicates by (v_b_gene, j_b_gene, cdr3_b_aa) before computing
# distances, so immLynx::runTCRdist returns a clone-level matrix smaller than
# length(td$barcodes) whenever the input has any repeated clonotypes. Expand
# it back to a cell-level matrix indexed by barcode by mapping each cell to
# its first-occurrence clone index. Verified empirically (2026-05-21) that
# tcrdist3 preserves first-occurrence order in its TCRrep dedup, so
# match(clone_key, clone_key[!duplicated(clone_key)]) gives the correct
# row index in the returned matrix. Downstream consumers (R/67, R/68)
# require cell-barcode rownames.
.tcrdist_to_cell_matrix <- function(td) {
  dmat <- td$distances$pw_beta %||% td$distances$pw_cdr3_b_aa
  if (is.null(dmat) || is.null(td$barcodes)) return(NULL)
  dmat <- as.matrix(dmat)
  n_cells <- length(td$barcodes)
  if (nrow(dmat) == n_cells) {
    rownames(dmat) <- td$barcodes
    colnames(dmat) <- td$barcodes
    return(dmat)
  }
  if (is.null(td$tcr_data)) {
    stop("tcrdist matrix has ", nrow(dmat), " rows but td$barcodes has ",
         n_cells, " entries, and td$tcr_data is NULL — cannot expand.")
  }
  ckey <- paste(td$tcr_data$v_b_gene, td$tcr_data$j_b_gene,
                td$tcr_data$cdr3_b_aa, sep = "|")
  ukey <- ckey[!duplicated(ckey)]
  if (length(ukey) != nrow(dmat)) {
    stop("tcrdist matrix dim (", nrow(dmat), ") does not match unique ",
         "(v,j,cdr3) tuples (", length(ukey),
         "); tcrdist3 dedup ordering assumption violated. Inspect ",
         "td$tcr_data and tcrdist3 version.")
  }
  idx <- match(ckey, ukey)
  out <- dmat[idx, idx, drop = FALSE]
  rownames(out) <- td$barcodes
  colnames(out) <- td$barcodes
  out
}

run_immlynx_tcrdist <- function(cfg) {
  if (!isTRUE(cfg$tcrdist$enable) || !isTRUE(cfg$steps$tcrdist)) {
    log_message("tcrdist (immLynx) disabled. Skipping.")
    return(invisible(TRUE))
  }
  if (!requireNamespace("immLynx", quietly = TRUE)) {
    log_message("immLynx not installed; skipping tcrdist step.")
    return(invisible(FALSE))
  }

  log_message("Starting immLynx tcrdist + embedding (Viral vs NIU intraocular)...")

  obj_path <- file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Integrated object not found. Skipping tcrdist.")
    return(invisible(TRUE))
  }
  obj <- readRDS(obj_path)

  out_tables <- file.path(cfg$paths$results_tables, "repertoire")
  out_obj    <- cfg$paths$results_objects
  ensure_dir(out_tables); ensure_dir(out_obj)

  # prepare_intraocular_trb() + trb_frame_to_sce() come from R/64_immgliph.R.
  trb <- prepare_intraocular_trb(obj, cfg)
  log_message("  Intraocular TRB rows (per cell): ", nrow(trb))
  if (nrow(trb) < 50) {
    log_message("  Too few TRB sequences. Aborting tcrdist.")
    return(invisible(FALSE))
  }

  # Deduplicate to unique (TRBV, TRBJ, CDR3b) tuples before embedding/UMAP.
  # Without this dedup, the per-cell SCE feeds ~70k rows into ESM-2 +
  # scater::runUMAP where many rows are identical CDR3b sequences; UMAP's
  # local repulsion then spreads each repeated sequence into a small
  # circular cluster of dots at near-identical coordinates (visible as
  # "rings" in the tcrdist UMAP). The tcrdist3 backend already
  # internally dedupes via match(v_b_gene, j_b_gene, cdr3_b_aa), so the
  # downstream cell-level distance matrix is recovered losslessly by
  # .tcrdist_to_cell_matrix() — we just move that dedup forward so the
  # embedding and UMAP are computed once per unique clonotype rather than
  # once per cell.
  clone_key <- paste(trb$TRBV %||% NA, trb$TRBJ %||% NA, trb$CDR3b,
                     sep = "|")
  first_per_clone <- !duplicated(clone_key)
  trb_unique <- trb[first_per_clone, , drop = FALSE]
  # Carry the per-cell -> unique-clone index for downstream consumers.
  cell_to_clone_idx <- match(clone_key, clone_key[first_per_clone])
  log_message("  Unique (TRBV, TRBJ, CDR3b) clonotypes: ", nrow(trb_unique),
              " (", round(100 * nrow(trb_unique) / nrow(trb), 1),
              "% of cells; ", nrow(trb) - nrow(trb_unique),
              " duplicate cells collapsed)")

  sce <- trb_frame_to_sce(trb_unique)

  # ---- Embedding (ESM-2 via basilisk) -------------------------------------
  model_name <- cfg$tcrdist$model_name %||% "facebook/esm2_t12_35M_UR50D"
  log_message("  Running runEmbeddings (model_name=", model_name, ")...")
  sce <- tryCatch(
    immLynx::runEmbeddings(sce, chains = "TRB",
                           model_name = model_name,
                           return_object = TRUE),
    error = function(e) {
      log_message("  runEmbeddings failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(sce)) return(invisible(FALSE))

  reds <- SingleCellExperiment::reducedDimNames(sce)
  log_message("  reducedDimNames after runEmbeddings: ",
              paste(reds, collapse = ", "))
  red_name <- intersect(c("tcr_esm", "ESM_TRB", "tcr_esm_TRB"), reds)[1]
  if (is.na(red_name) && length(reds) > 0) red_name <- reds[1]

  if (!is.na(red_name)) {
    log_message("  Running UMAP on '", red_name, "' embedding...")
    sce <- tryCatch(
      scater::runUMAP(sce, dimred = red_name, name = "UMAP_TCR"),
      error = function(e) {
        log_message("  UMAP failed: ", conditionMessage(e)); sce
      }
    )
  } else {
    log_message("  No reduction found post-runEmbeddings; UMAP skipped.")
  }
  umap_coords <- tryCatch(
    SingleCellExperiment::reducedDim(sce, "UMAP_TCR"),
    error = function(e) NULL
  )
  if (!is.null(umap_coords))
    log_message("  UMAP_TCR coords (per unique clonotype): ",
                nrow(umap_coords), " x ", ncol(umap_coords))

  # Expand per-clone UMAP back to per-cell so downstream consumers (R/67,
  # R/77 Panel F via GliphTcrdistJoint) can still join on cell barcode.
  # Cells from the same clonotype land at identical coords — visually they
  # overlap (no UMAP-induced rings from running on duplicate inputs).
  umap_per_clone <- umap_coords
  if (!is.null(umap_coords)) {
    umap_per_cell <- umap_coords[cell_to_clone_idx, , drop = FALSE]
    rownames(umap_per_cell) <- trb$barcode
    umap_coords <- umap_per_cell
    log_message("  UMAP_TCR coords expanded to per-cell: ",
                nrow(umap_coords), " rows (one per cell, ",
                "identical coords within clonotype).")
  }

  # ---- tcrdist distance matrix -------------------------------------------
  # Output: list(distances = list(pw_alpha, pw_beta, pw_cdr3_a_aa,
  # pw_cdr3_b_aa), barcodes, tcr_data).
  log_message("  Running runTCRdist (chains=beta)...")
  td <- tryCatch(
    immLynx::runTCRdist(sce, chains = "beta",
                        compute_distances = TRUE,
                        add_to_object     = FALSE),
    error = function(e) {
      log_message("  runTCRdist failed: ", conditionMessage(e)); NULL
    }
  )

  neigh_df    <- NULL
  per_subject <- NULL
  glm_res     <- NULL

  if (!is.null(td)) {
    log_message("  runTCRdist slots: ", paste(names(td), collapse = ", "))
    if (!is.null(td$distances))
      log_message("  distances slots: ",
                  paste(names(td$distances), collapse = ", "))

    td_barcodes <- td$barcodes
    dmat_clone <- tryCatch(.tcrdist_to_cell_matrix(td),
                           error = function(e) {
                             log_message("  tcrdist matrix expansion failed: ",
                                         conditionMessage(e)); NULL
                           })

    if (is.null(dmat_clone) || is.null(td_barcodes)) {
      log_message("  tcrdist returned no usable matrix (pw_beta missing).")
    } else {
      # dmat_clone is now per-unique-clonotype (because SCE was deduped above).
      # Expand to per-cell for downstream R/67-68 (which expect cell-barcode
      # rownames). Cells from the same clonotype get identical distance rows
      # to all other cells — correct, since tcrdist is defined at the clone
      # level. clone_bc_lookup maps each unique-clone barcode to its index.
      clone_bc_lookup <- match(td_barcodes, trb_unique$barcode)
      # cell_to_unique_row[i] = which row of dmat_clone corresponds to cell i.
      cell_to_unique_row <- match(cell_to_clone_idx, clone_bc_lookup)
      dmat <- dmat_clone[cell_to_unique_row, cell_to_unique_row, drop = FALSE]
      rownames(dmat) <- trb$barcode
      colnames(dmat) <- trb$barcode
      saveRDS(dmat, file.path(out_obj, "tcrdist_pw_beta.rds"))
      log_message("  Saved: tcrdist_pw_beta.rds (",
                  nrow(dmat), " x ", ncol(dmat),
                  " per-cell; expanded from ", nrow(dmat_clone),
                  " unique-clonotype rows)")

      K <- cfg$tcrdist$knn %||% 25
      K <- min(K, nrow(dmat) - 1)
      log_message("  KNN density (K=", K, ") on ", nrow(dmat), " cells...")

      trb_by_bc <- trb[match(rownames(dmat), trb$barcode), , drop = FALSE]
      phen_vec  <- trb_by_bc$Phenotype_2

      nn_density <- t(apply(dmat, 1, function(row) {
        ord <- order(row)
        nn  <- ord[2:(K + 1)]
        c(viral_nn = sum(phen_vec[nn] == "Viral", na.rm = TRUE),
          niu_nn   = sum(phen_vec[nn] == "NIU",   na.rm = TRUE))
      }))

      neigh_df <- data.frame(
        barcode            = rownames(dmat),
        Subject            = trb_by_bc$Subject,
        Etiology           = trb_by_bc$Etiology,
        Phenotype_2        = phen_vec,
        CDR3b              = trb_by_bc$CDR3b,
        HLA_B27_pathogenic = trb_by_bc$HLA_B27_pathogenic,
        viral_nn           = nn_density[, "viral_nn"],
        niu_nn             = nn_density[, "niu_nn"],
        K                  = K,
        viral_frac         = nn_density[, "viral_nn"] / K,
        niu_frac           = nn_density[, "niu_nn"]   / K,
        same_group_frac    = ifelse(
          phen_vec == "Viral",
          nn_density[, "viral_nn"] / K,
          nn_density[, "niu_nn"]   / K
        ),
        stringsAsFactors = FALSE
      )
      write.csv(neigh_df,
                file.path(out_tables, "tcrdist_neighborhoods.csv"),
                row.names = FALSE)
      log_message("  Saved: tcrdist_neighborhoods.csv")

      per_subject <- neigh_df |>
        dplyr::group_by(Subject, Etiology, Phenotype_2) |>
        dplyr::summarise(n_cells              = dplyr::n(),
                         mean_same_group_frac = mean(same_group_frac),
                         med_same_group_frac  = stats::median(same_group_frac),
                         .groups = "drop")
      write.csv(per_subject,
                file.path(out_tables, "tcrdist_per_subject.csv"),
                row.names = FALSE)
      log_message("  Saved: tcrdist_per_subject.csv")

      glm_res <- tryCatch({
        fit <- stats::lm(mean_same_group_frac ~ Phenotype_2 + log10(n_cells),
                         data = per_subject)
        cs  <- summary(fit)$coefficients
        data.frame(metric   = "mean_same_group_frac",
                   term     = rownames(cs),
                   estimate = cs[, "Estimate"],
                   std_err  = cs[, "Std. Error"],
                   t        = cs[, "t value"],
                   p        = cs[, "Pr(>|t|)"])
      }, error = function(e) {
        log_message("  GLM failed: ", conditionMessage(e)); NULL
      })
      if (!is.null(glm_res)) {
        write.csv(glm_res,
                  file.path(out_tables, "tcrdist_glm_results.csv"),
                  row.names = FALSE)
        log_message("  Saved: tcrdist_glm_results.csv")
      }
    }
  }

  # ---- clusTCR motif clusters --------------------------------------------
  log_message("  Running runClustTCR (mcl)...")
  cluster_df <- tryCatch({
    res_sce <- immLynx::runClustTCR(sce, chains = "TRB", method = "mcl",
                                    return_object = TRUE)
    if (inherits(res_sce, "SingleCellExperiment")) {
      cd <- SummarizedExperiment::colData(res_sce)
      cluster_col <- intersect(
        c("clustcr_TRB_cluster", "clustcr_cluster", "clustcr"),
        colnames(cd))[1]
      if (is.na(cluster_col)) {
        m <- grep("^clustcr", colnames(cd), value = TRUE)
        cluster_col <- if (length(m) > 0) m[1] else NA_character_
      }
      if (!is.na(cluster_col)) {
        data.frame(barcode = colnames(res_sce),
                   cluster = cd[[cluster_col]],
                   stringsAsFactors = FALSE)
      } else NULL
    } else as.data.frame(res_sce)
  }, error = function(e) {
    log_message("  runClustTCR failed: ", conditionMessage(e)); NULL
  })
  if (!is.null(cluster_df)) {
    # cluster_df is per-unique-clonotype (SCE was deduped above). Expand to
    # per-cell so downstream R/67 can join on cell barcode; cells from the
    # same clonotype inherit the same clusTCR cluster label.
    clone_to_cluster <- setNames(cluster_df$cluster, cluster_df$barcode)
    rep_bc <- trb_unique$barcode[cell_to_clone_idx]
    cluster_df <- data.frame(barcode = trb$barcode,
                             cluster = unname(clone_to_cluster[rep_bc]),
                             stringsAsFactors = FALSE)
    write.csv(cluster_df,
              file.path(out_tables, "clustcr_clusters.csv"),
              row.names = FALSE)
    log_message("  Saved: clustcr_clusters.csv (", nrow(cluster_df),
                " per-cell rows; expanded from ",
                length(unique(rep_bc)), " unique clonotypes)")
  }

  # ---- Persist for viz ---------------------------------------------------
  # umap = per-cell coords (backward-compatible with R/67, R/77 Panel F).
  # umap_per_clone = clone-level coords for analyses that want one dot per
  # unique TCR (e.g., when overlay-rendering motifs at the clonotype level).
  # trb_unique is the per-clonotype frame; cell_to_clone_idx is the
  # per-cell -> clonotype index for downstream consumers that need it.
  saveRDS(list(trb              = trb,
               trb_unique       = trb_unique,
               cell_to_clone_idx = cell_to_clone_idx,
               sce              = sce,
               umap             = umap_coords,
               umap_per_clone   = umap_per_clone,
               neigh            = neigh_df,
               clusters         = cluster_df,
               per_subject      = per_subject,
               glm              = glm_res),
          file.path(out_obj, "ImmLynxTcrdistResults.rds"))
  log_message("  Saved: ImmLynxTcrdistResults.rds (umap=",
              !is.null(umap_coords), ", neigh=", !is.null(neigh_df),
              ", clusters=", !is.null(cluster_df), ")")

  log_message("immLynx tcrdist analysis complete.")
  invisible(TRUE)
}
