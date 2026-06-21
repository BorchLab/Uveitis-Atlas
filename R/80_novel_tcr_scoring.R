# R/80_novel_tcr_scoring.R
# Score every GLIPH motif cluster as a candidate novel autoimmune-pathogenic
# motif by integrating Pgen, GLIPH cluster properties + enrichment, network
# centrality, public/expanded clone status, dominant T cell cluster + its
# MiloR DA logFC, the Figure-4 LR module engagement count, the HLA-B27 known-
# motif overlap, and a per-cell receptor Z-score over the autoimmune-biased
# receptors.
#
# Inputs (all already on disk after the upstream Phase 1d / 3d / 3e steps):
#   outputs/objects/ImmGLIPHResults.rds       (g$clusters, g$enrich)
#   outputs/tables/repertoire/gliph_cluster_properties.csv
#   outputs/tables/repertoire/gliph_enrichment_viral_vs_niu.csv
#   outputs/tables/repertoire/gliph_network_node_attrs.csv  (per-CDR3 centrality)
#   outputs/tables/repertoire/olga_pgen_per_clone.csv       (per-cell Pgen)
#   outputs/tables/repertoire/hla_b27_pathogenic_clones.csv
#   outputs/tables/eye/tcell/milo_da_cluster_calls.csv
#   outputs/tables/cross_compartment/liana_myeloid_to_tcell_combined.csv
#   outputs/objects/eye/tcell/IntegratedSeuratObject.rds    (for receptor Z + Leiden join)
#
# Outputs:
#   outputs/tables/repertoire/novel_tcr_candidates_ranked.csv
#   outputs/objects/NovelTcrCandidates.rds
#
# Entry: run_novel_tcr_scoring(cfg). Helpers prefixed .ntsc_*.

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------------------
# Load every input table once. Each is wrapped so a missing file degrades to
# NULL rather than killing the whole step — the scorer then NA-fills the
# affected feature and logs a warning.
# ---------------------------------------------------------------------------
.ntsc_load_inputs <- function(cfg) {
  rep_tbl <- file.path(get_target_paths(cfg, "all")$results_tables, "repertoire")
  rep_obj <- get_target_paths(cfg, "all")$results_objects
  t_tbl   <- get_target_paths(cfg, "tcell")$results_tables
  cc_tbl  <- file.path(cfg$paths$results_tables, "cross_compartment")

  read_csv_safe <- function(path) {
    if (!file.exists(path)) {
      log_message("  .ntsc_load_inputs: missing ", basename(path))
      return(NULL)
    }
    utils::read.csv(path, stringsAsFactors = FALSE)
  }

  inputs <- list(
    gliph_rds_path = file.path(rep_obj, "ImmGLIPHResults.rds"),
    gliph_clusters = read_csv_safe(file.path(rep_tbl, "gliph_clusters.csv")),
    gliph_props    = read_csv_safe(file.path(rep_tbl,
                                             "gliph_cluster_properties.csv")),
    gliph_enrich   = read_csv_safe(file.path(rep_tbl,
                                             "gliph_enrichment_viral_vs_niu.csv")),
    gliph_node     = read_csv_safe(file.path(rep_tbl,
                                             "gliph_network_node_attrs.csv")),
    pgen           = read_csv_safe(file.path(rep_tbl,
                                             "olga_pgen_per_clone.csv")),
    b27            = read_csv_safe(file.path(rep_tbl,
                                             "hla_b27_pathogenic_clones.csv")),
    milo           = read_csv_safe(file.path(t_tbl,
                                             "milo_da_cluster_calls.csv")),
    liana          = read_csv_safe(file.path(cc_tbl,
                                             "liana_myeloid_to_tcell_combined.csv")),
    tcell_obj_path = file.path(get_target_paths(cfg, "tcell")$results_objects,
                               "IntegratedSeuratObject.rds")
  )

  required <- c("gliph_clusters", "pgen")
  missing <- required[vapply(required, function(k) is.null(inputs[[k]]),
                             logical(1))]
  if (length(missing)) {
    stop(".ntsc_load_inputs: required inputs missing: ",
         paste(missing, collapse = ", "))
  }
  inputs
}

# ---------------------------------------------------------------------------
# Build the (cluster_id x CDR3b x cell) membership table by joining GLIPH
# memberships to the per-cell Pgen rows (the Pgen table is the per-cell
# universe — one row per barcode in the intraocular TRB filter). This is
# the single source of truth every downstream feature joins onto.
# ---------------------------------------------------------------------------
.ntsc_build_membership <- function(inputs) {
  cls  <- inputs$gliph_clusters
  pgen <- inputs$pgen

  # Drop empty CDR3 rows up front.
  cls  <- cls[!is.na(cls$CDR3b)  & nchar(cls$CDR3b) > 0, , drop = FALSE]
  pgen <- pgen[!is.na(pgen$CDR3b) & nchar(pgen$CDR3b) > 0, , drop = FALSE]

  # A CDR3b can belong to multiple GLIPH clusters (multi-motif hits). Keep
  # the (cluster_id, CDR3b) edges as-is and many-to-many join to Pgen rows.
  pairs <- dplyr::distinct(cls[, c("cluster_id", "CDR3b")])

  mem <- pgen |>
    dplyr::inner_join(pairs, by = "CDR3b",
                      relationship = "many-to-many")

  # Per-CDR3 cell count within each (cluster_id, CDR3b) bin.
  mem <- mem |>
    dplyr::group_by(.data$cluster_id, .data$CDR3b) |>
    dplyr::mutate(n_cells_in_clone = dplyr::n()) |>
    dplyr::ungroup()

  log_message(sprintf("  .ntsc_build_membership: %d cells across %d (cluster x CDR3b) pairs in %d clusters.",
                      nrow(mem), nrow(pairs),
                      dplyr::n_distinct(mem$cluster_id)))
  mem
}

# ---------------------------------------------------------------------------
# Per-cluster features that do NOT need the Seurat object. Each feature is
# computed defensively — if its source input is NULL the column is NA and
# the scorer treats it as a neutral (median) rank.
# ---------------------------------------------------------------------------
.ntsc_per_cluster_features <- function(mem, inputs, cfg) {
  expanded_cutoff <- cfg$tcr_advanced$fig5$novel_tcr_discovery$expanded_clone_cutoff %||% 3L
  min_cluster_size <- cfg$tcr_advanced$fig5$novel_tcr_discovery$min_cluster_size %||% 2L

  # ---- Membership-derived counts -----------------------------------------
  base <- mem |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::summarise(
      n_cdr3     = dplyr::n_distinct(.data$CDR3b),
      n_subjects = dplyr::n_distinct(.data$Subject),
      n_cells    = dplyr::n(),
      n_expanded_clones = dplyr::n_distinct(.data$CDR3b[
        .data$n_cells_in_clone >= expanded_cutoff]),
      frac_in_NIU = mean(.data$Phenotype_2 == "NIU", na.rm = TRUE),
      subjects_csv = paste(sort(unique(.data$Subject)), collapse = ","),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$n_cdr3 >= min_cluster_size)

  if (!nrow(base)) {
    log_message("  .ntsc_per_cluster_features: no clusters pass min_cluster_size; aborting.")
    return(base)
  }

  # ---- NIU vs Viral Fisher per cluster -----------------------------------
  pgen <- inputs$pgen
  pheno <- pgen$Phenotype_2
  N_NIU_total <- sum(pheno == "NIU",   na.rm = TRUE)
  N_VIR_total <- sum(pheno == "Viral", na.rm = TRUE)
  fisher_one <- function(cid) {
    mm <- mem[mem$cluster_id == cid, , drop = FALSE]
    n_in_NIU  <- sum(mm$Phenotype_2 == "NIU",   na.rm = TRUE)
    n_in_VIR  <- sum(mm$Phenotype_2 == "Viral", na.rm = TRUE)
    n_out_NIU <- max(0, N_NIU_total - n_in_NIU)
    n_out_VIR <- max(0, N_VIR_total - n_in_VIR)
    mat <- matrix(c(n_in_NIU, n_in_VIR, n_out_NIU, n_out_VIR), nrow = 2)
    p <- tryCatch(stats::fisher.test(mat)$p.value,
                  error = function(e) NA_real_)
    p
  }
  base$niu_vs_viral_fisher_p <- vapply(base$cluster_id, fisher_one,
                                       numeric(1))
  base$niu_vs_viral_fisher_FDR <- stats::p.adjust(base$niu_vs_viral_fisher_p,
                                                  method = "BH")

  # ---- Pgen per-cluster median/min ---------------------------------------
  pgen_safe <- mem |>
    dplyr::filter(is.finite(.data$log10_Pgen)) |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::summarise(
      median_log10_Pgen = stats::median(.data$log10_Pgen),
      min_log10_Pgen    = min(.data$log10_Pgen),
      .groups = "drop"
    )
  base <- dplyr::left_join(base, pgen_safe, by = "cluster_id")
  # Floor handling: replace -Inf (numerically 0 Pgen) with min(finite)-0.5.
  pgen_finite <- inputs$pgen$log10_Pgen[is.finite(inputs$pgen$log10_Pgen)]
  pgen_floor <- if (length(pgen_finite)) min(pgen_finite) - 0.5 else -20
  base$pgen_floored <- !is.finite(base$median_log10_Pgen) |
                       !is.finite(base$min_log10_Pgen)
  base$median_log10_Pgen[!is.finite(base$median_log10_Pgen)] <- pgen_floor
  base$min_log10_Pgen[!is.finite(base$min_log10_Pgen)]       <- pgen_floor

  # ---- GLIPH cluster properties join ------------------------------------
  if (!is.null(inputs$gliph_props)) {
    gp <- inputs$gliph_props[, c("cluster_id", "fisher.score", "OvE",
                                 "total.score")]
    gp$gliph_fisher_score <- gp$fisher.score
    gp$gliph_OvE          <- gp$OvE
    gp$gliph_total_score  <- gp$total.score
    gp <- gp[, c("cluster_id", "gliph_fisher_score", "gliph_OvE",
                 "gliph_total_score")]
    base <- dplyr::left_join(base, gp, by = "cluster_id")
  } else {
    base$gliph_fisher_score <- NA_real_
    base$gliph_OvE          <- NA_real_
    base$gliph_total_score  <- NA_real_
  }

  # ---- NIU enrichment FDR (from gliph_enrichment_viral_vs_niu) ----------
  if (!is.null(inputs$gliph_enrich)) {
    en <- inputs$gliph_enrich
    en_niu <- en[en$direction == "NIU_enriched", c("cluster_id", "FDR")]
    en_niu$niu_enrichment_FDR <- en_niu$FDR
    en_niu$FDR <- NULL
    base <- dplyr::left_join(base, en_niu, by = "cluster_id")
  } else {
    base$niu_enrichment_FDR <- NA_real_
  }

  # ---- Network centrality (per-CDR3 attrs aggregated to cluster) --------
  if (!is.null(inputs$gliph_node)) {
    nd <- inputs$gliph_node |>
      dplyr::group_by(.data$cluster_id) |>
      dplyr::summarise(
        network_eigen_centrality = mean(.data$eigen_centrality,
                                        na.rm = TRUE),
        network_degree           = mean(.data$degree, na.rm = TRUE),
        network_betweenness      = mean(.data$betweenness, na.rm = TRUE),
        .groups = "drop"
      )
    base <- dplyr::left_join(base, nd, by = "cluster_id")
  } else {
    base$network_eigen_centrality <- NA_real_
    base$network_degree           <- NA_real_
    base$network_betweenness      <- NA_real_
  }

  # ---- B27 known-motif overlap -------------------------------------------
  if (!is.null(inputs$b27)) {
    b27_cdr3 <- unique(inputs$b27$CDR3b)
    b27_by_cluster <- mem |>
      dplyr::group_by(.data$cluster_id) |>
      dplyr::summarise(
        frac_b27_pathogenic = mean(unique(.data$CDR3b) %in% b27_cdr3),
        .groups = "drop"
      )
    base <- dplyr::left_join(base, b27_by_cluster, by = "cluster_id")
  } else {
    base$frac_b27_pathogenic <- 0
  }
  base$is_b27_known <- !is.na(base$frac_b27_pathogenic) &
                       base$frac_b27_pathogenic >= 0.5

  base
}

# ---------------------------------------------------------------------------
# Dominant T cell cluster + Milo DA logFC + LIANA LR engagement. Needs a
# barcode -> knn.leiden.cluster map, which we pull lazily from the tcell
# Seurat object metadata. The object is loaded once and released.
# ---------------------------------------------------------------------------
.ntsc_attach_cluster_context <- function(features, mem, inputs, cfg) {
  if (!file.exists(inputs$tcell_obj_path)) {
    log_message("  .ntsc_attach_cluster_context: tcell object missing; cluster context NA.")
    features$dominant_tcell_cluster <- NA_character_
    features$dominant_cluster_frac  <- NA_real_
    features$dominant_cluster_DA_logFC <- NA_real_
    features$n_LR_autoimmune_targeting_dominant_cluster <- NA_integer_
    return(features)
  }
  bc_to_clu <- local({
    obj <- readRDS(inputs$tcell_obj_path)
    md  <- obj@meta.data
    clu_col <- if ("knn.leiden.cluster" %in% colnames(md))
                 "knn.leiden.cluster" else NULL
    if (is.null(clu_col)) {
      log_message("  .ntsc_attach_cluster_context: knn.leiden.cluster column missing.")
      return(NULL)
    }
    out <- setNames(as.character(md[[clu_col]]), colnames(obj))
    rm(obj); invisible(gc())
    out
  })
  if (is.null(bc_to_clu)) {
    features$dominant_tcell_cluster <- NA_character_
    features$dominant_cluster_frac  <- NA_real_
    features$dominant_cluster_DA_logFC <- NA_real_
    features$n_LR_autoimmune_targeting_dominant_cluster <- NA_integer_
    return(features)
  }

  mem$knn_cluster <- unname(bc_to_clu[mem$barcode])
  mem_has <- mem[!is.na(mem$knn_cluster), , drop = FALSE]

  dom <- mem_has |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::summarise(
      dominant_tcell_cluster = {
        tbl <- sort(table(.data$knn_cluster), decreasing = TRUE)
        names(tbl)[1]
      },
      dominant_cluster_frac = {
        tbl <- table(.data$knn_cluster)
        max(tbl) / sum(tbl)
      },
      .groups = "drop"
    )
  features <- dplyr::left_join(features, dom, by = "cluster_id")
  # Mixed-cluster flag: treat low-fraction dominants as "mixed".
  mixed <- !is.na(features$dominant_cluster_frac) &
           features$dominant_cluster_frac < 0.4
  features$dominant_tcell_cluster[mixed] <- "mixed"

  # ---- MiloR DA logFC join -----------------------------------------------
  if (!is.null(inputs$milo)) {
    milo <- inputs$milo
    milo$cluster <- as.character(milo$cluster)
    # Signed logFC: positive == NIU-up. milo$direction == "NIU_up" already
    # uses that convention for mean_logFC.
    milo$signed_logFC <- ifelse(milo$direction == "NIU_up",
                                abs(milo$mean_logFC),
                                ifelse(milo$direction == "Viral_up",
                                       -abs(milo$mean_logFC),
                                       milo$mean_logFC))
    milo_join <- milo[, c("cluster", "signed_logFC")]
    names(milo_join) <- c("dominant_tcell_cluster", "dominant_cluster_DA_logFC")
    features <- dplyr::left_join(features, milo_join,
                                 by = "dominant_tcell_cluster")
  } else {
    features$dominant_cluster_DA_logFC <- NA_real_
  }

  # ---- LIANA autoimmune LR engagement count ------------------------------
  if (!is.null(inputs$liana)) {
    lia <- inputs$liana
    if (!"aggregate_rank_NIU" %in% colnames(lia)) {
      log_message("  .ntsc_attach_cluster_context: LIANA aggregate_rank_NIU missing; skipping LR feature.")
      features$n_LR_autoimmune_targeting_dominant_cluster <- NA_integer_
    } else {
      # Top quartile of NIU rank (lower rank = stronger pair) AND disease_bias
      # favors NIU (positive logfc).
      rk <- lia$aggregate_rank_NIU
      q1 <- stats::quantile(rk, 0.25, na.rm = TRUE)
      keep <- !is.na(rk) & rk <= q1 &
              !is.na(lia$disease_bias_logfc) & lia$disease_bias_logfc > 0
      lr_keep <- lia[keep, , drop = FALSE]
      lr_count <- table(lr_keep$target)
      lr_df <- data.frame(
        target  = names(lr_count),
        n_LR    = as.integer(lr_count),
        stringsAsFactors = FALSE
      )
      features$liana_target_key <- ifelse(is.na(features$dominant_tcell_cluster) |
                                          features$dominant_tcell_cluster == "mixed",
                                          NA_character_,
                                          paste0("tcell_", features$dominant_tcell_cluster))
      features <- dplyr::left_join(features,
                                   dplyr::rename(lr_df,
                                                 liana_target_key = "target",
                                                 n_LR_autoimmune_targeting_dominant_cluster = "n_LR"),
                                   by = "liana_target_key")
      features$liana_target_key <- NULL
      features$n_LR_autoimmune_targeting_dominant_cluster[is.na(
        features$n_LR_autoimmune_targeting_dominant_cluster)] <- 0L
      features$n_LR_autoimmune_targeting_dominant_cluster[mixed] <- NA_integer_
    }
  } else {
    features$n_LR_autoimmune_targeting_dominant_cluster <- NA_integer_
  }

  features
}

# ---------------------------------------------------------------------------
# Per-cell receptor Z-score (logcounts) averaged to per-CDR3b -> per-cluster.
# Loads the tcell Seurat object once, restricts to membership barcodes,
# Z-scales each receptor across that subset, row-means to per-cell, then
# aggregates to cluster. Skip cleanly if the object or any of the receptors
# are missing.
# ---------------------------------------------------------------------------
.ntsc_receptor_z <- function(features, mem, inputs, cfg) {
  if (!file.exists(inputs$tcell_obj_path)) {
    log_message("  .ntsc_receptor_z: tcell object missing; receptor_Z_score NA.")
    features$receptor_Z_score <- NA_real_
    return(features)
  }
  genes <- cfg$tcr_advanced$fig5$novel_tcr_discovery$receptor_z_genes %||%
           c("CTLA4","CXCR3","IL2RA","CD44","ITGA4","CD69","PTPRC")
  bcs <- unique(mem$barcode)
  if (!length(bcs)) {
    features$receptor_Z_score <- NA_real_
    return(features)
  }
  z_by_bc <- local({
    obj <- readRDS(inputs$tcell_obj_path)
    bcs_have <- intersect(bcs, colnames(obj))
    if (!length(bcs_have)) {
      log_message("  .ntsc_receptor_z: no membership barcodes match the tcell object.")
      rm(obj); invisible(gc())
      return(NULL)
    }
    obj_sub <- obj[, bcs_have]
    rm(obj); invisible(gc())
    have_genes <- intersect(genes, rownames(obj_sub))
    missing_g <- setdiff(genes, have_genes)
    if (length(missing_g)) {
      log_message("  .ntsc_receptor_z: skipping missing receptors: ",
                  paste(missing_g, collapse = ", "))
    }
    if (!length(have_genes)) {
      rm(obj_sub); invisible(gc())
      return(NULL)
    }
    mat <- Seurat::FetchData(obj_sub, vars = have_genes,
                             layer = "data")
    # Z per gene across membership subset.
    mat_z <- scale(mat)
    mat_z[!is.finite(mat_z)] <- 0
    per_cell <- rowMeans(mat_z, na.rm = TRUE)
    rm(obj_sub, mat, mat_z); invisible(gc())
    setNames(per_cell, bcs_have)
  })

  if (is.null(z_by_bc)) {
    features$receptor_Z_score <- NA_real_
    return(features)
  }
  mem$receptor_Z <- unname(z_by_bc[mem$barcode])
  zc <- mem |>
    dplyr::filter(!is.na(.data$receptor_Z)) |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::summarise(receptor_Z_score = mean(.data$receptor_Z),
                     .groups = "drop")
  features <- dplyr::left_join(features, zc, by = "cluster_id")
  features
}

# ---------------------------------------------------------------------------
# Composite score: per feature, rank within the candidate set with the
# correct direction, min-max scale to [0,1], weighted sum. NAs collapse to
# the neutral midpoint (0.5).
# ---------------------------------------------------------------------------
.ntsc_composite <- function(features, cfg) {
  weights <- cfg$tcr_advanced$fig5$novel_tcr_discovery$feature_weights %||% list()

  # (feature column, direction, weight key) tuples. direction = +1 means
  # higher-is-better; -1 means lower-is-better.
  #
  # Deliberately excluded: gliph_fisher_score and gliph_total_score. Those
  # are per-cluster internal GLIPH metrics that inflate for very small
  # clusters with perfect convergence and pull the composite away from
  # the biologically meaningful signal (cell count + disease bias).
  feature_spec <- list(
    list(col = "median_log10_Pgen", dir = -1, w = "log_pgen"),
    list(col = "min_log10_Pgen",    dir = -1, w = "log_pgen"),
    list(col = "niu_enrichment_FDR",dir = -1, w = "niu_enrichment"),
    list(col = "frac_in_NIU",       dir = +1, w = "frac_in_NIU"),
    list(col = "niu_vs_viral_fisher_FDR", dir = -1, w = "fisher_p"),
    list(col = "n_subjects",        dir = +1, w = "n_subjects"),
    list(col = "n_expanded_clones", dir = +1, w = "n_subjects"),
    list(col = "n_cells",           dir = +1, w = "n_cells"),
    list(col = "network_eigen_centrality", dir = +1, w = "network_centrality"),
    list(col = "dominant_cluster_DA_logFC", dir = +1, w = "dominant_DA_logFC"),
    list(col = "n_LR_autoimmune_targeting_dominant_cluster", dir = +1, w = "LR_targeting"),
    list(col = "receptor_Z_score",  dir = +1, w = "receptor_Z")
  )

  scaled <- vapply(feature_spec, function(sp) {
    v <- features[[sp$col]]
    if (is.null(v)) return(rep(0.5, nrow(features)))
    if (sp$dir == -1) v <- -v
    r <- rank(v, ties.method = "average", na.last = "keep")
    r_min <- min(r, na.rm = TRUE); r_max <- max(r, na.rm = TRUE)
    if (!is.finite(r_min) || !is.finite(r_max) || r_max == r_min) {
      out <- rep(0.5, length(r))
    } else {
      out <- (r - r_min) / (r_max - r_min)
    }
    out[is.na(out)] <- 0.5
    out
  }, numeric(nrow(features)))

  w_vec <- vapply(feature_spec, function(sp) {
    as.numeric(weights[[sp$w]] %||% 1)
  }, numeric(1))

  features$composite_score <- as.numeric(scaled %*% w_vec) / sum(w_vec)
  features$composite_rank  <- rank(-features$composite_score,
                                   ties.method = "min")
  features <- features[order(features$composite_rank), , drop = FALSE]
  features
}

# ---------------------------------------------------------------------------
# Top-level entry. Returns the candidate table invisibly and writes both
# CSV and RDS.
# ---------------------------------------------------------------------------
run_novel_tcr_scoring <- function(cfg) {
  log_message("Starting novel TCR candidate scoring...")
  inputs <- .ntsc_load_inputs(cfg)
  mem    <- .ntsc_build_membership(inputs)
  features <- .ntsc_per_cluster_features(mem, inputs, cfg)
  if (!nrow(features)) {
    log_message("  run_novel_tcr_scoring: no candidates after min_cluster_size filter; nothing to write.")
    return(invisible(NULL))
  }
  features <- .ntsc_attach_cluster_context(features, mem, inputs, cfg)
  features <- .ntsc_receptor_z(features, mem, inputs, cfg)
  features <- .ntsc_composite(features, cfg)

  # Carry the motif string for reporting.
  if (!is.null(inputs$gliph_props)) {
    motif_lk <- setNames(inputs$gliph_props$motif,
                         inputs$gliph_props$cluster_id)
    features$motif <- unname(motif_lk[features$cluster_id])
  }

  # ---- Outputs -----------------------------------------------------------
  out_tbl <- file.path(get_target_paths(cfg, "all")$results_tables,
                       "repertoire", "novel_tcr_candidates_ranked.csv")
  ensure_dir(dirname(out_tbl))
  utils::write.csv(features, out_tbl, row.names = FALSE)
  log_message("  Wrote ", basename(out_tbl), " (", nrow(features), " rows).")

  out_rds <- file.path(get_target_paths(cfg, "all")$results_objects,
                       "NovelTcrCandidates.rds")
  ensure_dir(dirname(out_rds))
  saveRDS(features, out_rds)
  log_message("  Wrote ", basename(out_rds), ".")

  # ---- Sanity log --------------------------------------------------------
  log_message(sprintf("  Subjects represented: %d unique.",
                      length(unique(unlist(strsplit(features$subjects_csv, ","))))))
  log_message(sprintf("  B27-known motifs flagged: %d.", sum(features$is_b27_known)))
  if (any(features$is_b27_known)) {
    b27_top <- min(features$composite_rank[features$is_b27_known])
    log_message(sprintf("  Top B27-known motif sits at composite_rank = %d.", b27_top))
  }

  invisible(features)
}
