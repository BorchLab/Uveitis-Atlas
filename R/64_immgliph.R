# R/64_immgliph.R
# GLIPH2 convergence-group discovery on intraocular TRB repertoires.
# Compares Viral vs NIU using immGLIPH (BorchLab/immGLIPH) with bootstrap
# resampling to control for cell-number disparity (35-7,213 TCR+ per sample).
# Outputs: gliph_clusters.csv, gliph_enrichment_viral_vs_niu.csv,
# gliph_motif_logos.rds (group-level position matrices for ggseqlogo).

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Shared preprocessor (sourced by 65 + 66 too) ---------------------------
# Builds a clean intraocular TRB data frame from the integrated Seurat object.
# Filters: T cells, eye tissue, viral or NIU, in-frame productive TRB,
# single chain (no semicolons), length 8 - 60 aa.
prepare_intraocular_trb <- function(obj, cfg) {
  meta <- obj@meta.data

  cb_col <- resolve_celltype_broad(meta)
  if (is.null(cb_col)) stop("No celltype_broad column found in metadata.")

  trb_full  <- stringr::str_split(meta$CTaa,  "_", simplify = TRUE)[, 2]
  gene_full <- stringr::str_split(meta$CTgene,"_", simplify = TRUE)[, 2]
  trbv <- stringr::str_extract(gene_full, "TRBV[0-9A-Z\\-]+")
  trbj <- stringr::str_extract(gene_full, "TRBJ[0-9A-Z\\-]+")

  min_len <- cfg$gliph$min_cdr3_length %||% 8

  keep <- meta[[cb_col]] == "T cell" &
          meta$Tissue_1 == "Eye" &
          !is.na(meta$CTaa) &
          meta$Phenotype_2 %in% c("Viral", "NIU") &
          !is.na(trb_full) &
          trb_full != "" &
          !is.na(trb_full) &
          !grepl(";", trb_full) &
          nchar(trb_full) >= min_len &
          nchar(trb_full) <= 60

  keep[is.na(keep)] <- FALSE

  hla_b27_flag <- flag_hla_b27_pathogenic_tcr(meta$CTgene[keep],
                                              meta$CTaa[keep])

  data.frame(
    barcode             = colnames(obj)[keep],
    Subject             = meta$Subject[keep],
    Sample              = meta$Sample[keep] %||% NA,
    Etiology            = meta$Etiology[keep],
    Phenotype_2         = meta$Phenotype_2[keep],
    TRBV                = trbv[keep],
    CDR3b               = trb_full[keep],
    TRBJ                = trbj[keep],
    clone_id            = meta$CTstrict[keep],
    CTaa                = meta$CTaa[keep],
    CTgene              = meta$CTgene[keep],
    HLA_B27_pathogenic  = hla_b27_flag,
    stringsAsFactors    = FALSE
  )
}

# ---- Shared SCE builder (used by R/65 + R/66) -----------------------------
# immLynx functions consume SingleCellExperiment objects with scRepertoire
# colData (CTaa, CTgene, CTstrict). We carry the *original* CTaa/CTgene
# (not reconstructed stubs) so immApex::getIR() can parse both chains.
trb_frame_to_sce <- function(trb) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("SingleCellExperiment required.")

  cd <- S4Vectors::DataFrame(
    barcode            = trb$barcode,
    Subject            = trb$Subject,
    Etiology           = trb$Etiology,
    Phenotype_2        = trb$Phenotype_2,
    TRBV               = trb$TRBV,
    TRBJ               = trb$TRBJ,
    CDR3b              = trb$CDR3b,
    CTaa               = trb$CTaa,
    CTgene             = trb$CTgene,
    CTstrict           = trb$clone_id,
    clone_id           = trb$clone_id,
    HLA_B27_pathogenic = trb$HLA_B27_pathogenic
  )
  rownames(cd) <- trb$barcode

  m <- Matrix::Matrix(0, nrow = 1, ncol = nrow(trb), sparse = TRUE)
  colnames(m) <- trb$barcode
  rownames(m) <- "placeholder"

  SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = m),
    colData = cd
  )
}

# ---- Bootstrap helper (Fisher OR per GLIPH cluster across B subsamples) ----
.bootstrap_gliph_enrichment <- function(trb, gliph_clusters,
                                        n_iter = 200, seed = 42) {
  set.seed(seed)
  # gliph_clusters is a long table of (cluster_id, CDR3b) memberships.
  # We resample CDR3 rows from each Phenotype_2 group with replacement to a
  # common size (the smaller of the two), then recompute Fisher's OR per cluster.
  n_v <- sum(trb$Phenotype_2 == "Viral")
  n_n <- sum(trb$Phenotype_2 == "NIU")
  n_sub <- min(n_v, n_n)
  if (n_sub < 100) {
    log_message("  WARN: bootstrap pool < 100 (n=", n_sub,
                "). Enrichment will be unstable.")
  }

  clusters <- unique(gliph_clusters$cluster_id)
  cluster_members <- split(gliph_clusters$CDR3b, gliph_clusters$cluster_id)

  results <- vector("list", n_iter)
  for (b in seq_len(n_iter)) {
    iv <- sample(which(trb$Phenotype_2 == "Viral"), n_sub, replace = TRUE)
    in_ <- sample(which(trb$Phenotype_2 == "NIU"),   n_sub, replace = TRUE)
    cdr3_v <- trb$CDR3b[iv]
    cdr3_n <- trb$CDR3b[in_]
    rows <- lapply(clusters, function(cl) {
      mem <- cluster_members[[as.character(cl)]]
      a <- sum(cdr3_v %in% mem)
      c <- sum(cdr3_n %in% mem)
      b_ <- n_sub - a
      d  <- n_sub - c
      if (a + c < 3) return(NULL)
      ft <- suppressWarnings(fisher.test(matrix(c(a, c, b_, d), nrow = 2)))
      data.frame(cluster_id = cl,
                 iter       = b,
                 n_viral    = a,
                 n_niu      = c,
                 OR         = unname(ft$estimate),
                 p          = ft$p.value)
    })
    results[[b]] <- do.call(rbind, rows)
  }
  do.call(rbind, results)
}

# ---- Main entry --------------------------------------------------------------
run_immgliph <- function(cfg) {
  if (!isTRUE(cfg$gliph$enable) || !isTRUE(cfg$steps$gliph)) {
    log_message("GLIPH analysis disabled. Skipping.")
    return(invisible(TRUE))
  }
  if (!requireNamespace("immGLIPH", quietly = TRUE)) {
    log_message("immGLIPH not installed; cannot run GLIPH step.")
    return(invisible(FALSE))
  }

  log_message("Starting immGLIPH (Viral vs NIU intraocular TRB)...")

  obj_path <- file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Integrated object not found. Skipping GLIPH.")
    return(invisible(TRUE))
  }
  obj <- readRDS(obj_path)

  out_tables <- file.path(cfg$paths$results_tables, "repertoire")
  out_obj    <- cfg$paths$results_objects
  ensure_dir(out_tables); ensure_dir(out_obj)

  trb <- prepare_intraocular_trb(obj, cfg)
  log_message("  Intraocular TRB rows: ", nrow(trb),
              " (Viral=", sum(trb$Phenotype_2 == "Viral"),
              ", NIU=",   sum(trb$Phenotype_2 == "NIU"), ")")

  # HLA-B27 pathogenic signature (TRAV21 + CDR3b [YF]S[TS] flank motif)
  n_b27 <- sum(trb$HLA_B27_pathogenic, na.rm = TRUE)
  log_message("  HLA-B27 pathogenic candidates: ", n_b27,
              " cells (", length(unique(trb$CDR3b[trb$HLA_B27_pathogenic])),
              " unique CDR3-beta, ",
              length(unique(trb$Subject[trb$HLA_B27_pathogenic])),
              " subjects)")
  if (n_b27 > 0) {
    b27_df <- trb[trb$HLA_B27_pathogenic, , drop = FALSE]
    utils::write.csv(b27_df,
                     file.path(out_tables, "hla_b27_pathogenic_clones.csv"),
                     row.names = FALSE)
    log_message("  Saved: hla_b27_pathogenic_clones.csv")
  }

  if (nrow(trb) < 200) {
    log_message("  Too few TRB sequences (<200). GLIPH motif discovery aborted.")
    return(invisible(FALSE))
  }

  # immGLIPH::runGLIPH recognizes columns CDR3b, TRBV, patient, HLA, counts.
  # 1) Drop rows with NA / non-canonical TRBV.
  # 2) Restrict CDR3-beta to standard amino acids only (no *, _, X, lowercase)
  #    -- some scRepertoire outputs include stop codons or "Multi" placeholders
  #    that crash GLIPH's internal padding step with a names() length error.
  # 3) Aggregate to unique (subject, CDR3b, TRBV) with `counts` so GLIPH gets
  #    ~10-15k rows instead of ~60k -- much faster, less prone to chunk-level
  #    errors, and keeps publicity (per-subject) tracking intact.
  AA_RE <- "^[ACDEFGHIKLMNPQRSTVWY]+$"
  gliph_input <- trb |>
    dplyr::filter(!is.na(TRBV), nzchar(TRBV), grepl("^TRBV", TRBV),
                  !is.na(CDR3b), grepl(AA_RE, CDR3b)) |>
    dplyr::group_by(Subject, CDR3b, TRBV) |>
    dplyr::summarise(counts = dplyr::n(), .groups = "drop") |>
    dplyr::rename(patient = Subject) |>
    as.data.frame()

  log_message("  GLIPH input after sanitization: ", nrow(gliph_input),
              " unique (subject, CDR3b, TRBV) rows ",
              "(dropped ", nrow(trb) - sum(gliph_input$counts),
              " malformed rows).")
  if (nrow(gliph_input) < 50) {
    log_message("  Too few clean rows for GLIPH. Aborting.")
    return(invisible(FALSE))
  }

  method   <- cfg$gliph$method     %||% "gliph2"
  sim_dep  <- cfg$gliph$n_sims     %||% 1000
  refdb    <- cfg$gliph$refdb_beta %||% "human_v2.0_CD48"
  n_cores  <- cfg$gliph$n_cores    %||% 0
  if (!is.numeric(n_cores) || n_cores < 1)
    n_cores <- max(1, parallel::detectCores() - 1)

  .gliph_call <- function(input_df, cores) {
    immGLIPH::runGLIPH(
      input_df,
      method     = method,
      chains     = "TRB",
      refdb_beta = refdb,
      sim_depth  = sim_dep,
      n_cores    = cores,
      verbose    = FALSE
    )
  }

  log_message("  Running immGLIPH::runGLIPH (method=", method,
              ", sim_depth=", sim_dep, ", refdb_beta=", refdb,
              ", n_cores=", n_cores, ")...")

  gliph_res <- tryCatch(
    .gliph_call(gliph_input, n_cores),
    error = function(e) {
      log_message("  immGLIPH parallel call failed: ", conditionMessage(e))
      if (n_cores > 1) {
        log_message("  Retrying serially (n_cores=1) to surface the underlying error...")
        tryCatch(.gliph_call(gliph_input, 1L),
                 error = function(e2) {
                   log_message("  Serial retry also failed: ",
                               conditionMessage(e2))
                   NULL
                 })
      } else NULL
    }
  )
  if (is.null(gliph_res)) return(invisible(FALSE))

  # immGLIPH return slots (per documentation):
  #   cluster_properties : data.frame, one row per convergence group, with
  #     columns: type, tag (== cluster_id), cluster_size,
  #     unique_cdr3_sample, unique_cdr3_ref, OvE, fisher.score,
  #     total.score, network.size.score, cdr3.length.score, vgene.score,
  #     members (space-separated CDR3 string)
  #   cluster_list : named list (names == tag) of data.frames with columns
  #     seq_ID, CDR3b, TRBV, patient, ultCDR3b
  #
  # We flatten cluster_list to long format (one row per CDR3 in each group)
  # for the per-CDR3 enrichment downstream, and keep cluster_properties
  # for cluster-level motif + score reporting.
  cluster_list  <- gliph_res$cluster_list      %||% list()
  cluster_props <- gliph_res$cluster_properties %||% data.frame()

  if (length(cluster_list) == 0 || nrow(cluster_props) == 0) {
    log_message("  immGLIPH returned no convergence groups ",
                "(cluster_list empty). ",
                "Try a higher sim_depth (>= 1000) or check input size.")
    return(invisible(FALSE))
  }

  cluster_df <- do.call(rbind, lapply(names(cluster_list), function(tag) {
    d <- cluster_list[[tag]]
    if (is.null(d) || nrow(d) == 0) return(NULL)
    data.frame(cluster_id = tag,
               CDR3b      = d$CDR3b,
               TRBV       = d$TRBV,
               patient    = d$patient,
               stringsAsFactors = FALSE)
  }))

  # Cluster-level motif: tag format is "<motif>_<int>_<int>" for local
  # patterns and "global-<...>" or similar for GLIPH2 global hits.
  cluster_props$cluster_id <- cluster_props$tag
  cluster_props$motif      <- sub("_.*", "", cluster_props$tag)

  min_grp <- cfg$gliph$min_group_size %||% 3
  grp_n   <- table(cluster_df$cluster_id)
  keep_g  <- names(grp_n)[grp_n >= min_grp]
  cluster_df    <- cluster_df[cluster_df$cluster_id %in% keep_g, ,
                              drop = FALSE]
  cluster_props <- cluster_props[cluster_props$cluster_id %in% keep_g, ,
                                 drop = FALSE]
  log_message("  Retained ", length(keep_g),
              " convergence groups (>= ", min_grp, " members; ",
              nrow(cluster_df), " CDR3-beta memberships).")

  # Carry motif into per-CDR3 table so downstream viz that joins on
  # cluster_id has the motif label available.
  cluster_df <- dplyr::left_join(
    cluster_df,
    cluster_props[, c("cluster_id", "motif", "type")],
    by = "cluster_id"
  )

  write.csv(cluster_df,
            file.path(out_tables, "gliph_clusters.csv"),
            row.names = FALSE)
  write.csv(cluster_props,
            file.path(out_tables, "gliph_cluster_properties.csv"),
            row.names = FALSE)
  log_message("  Saved: gliph_clusters.csv, gliph_cluster_properties.csv")

  # ---- Bootstrap Viral vs NIU enrichment -----------------------------------
  n_iter <- cfg$gliph$bootstrap_iterations %||% 200
  log_message("  Bootstrap enrichment (B=", n_iter, ")...")
  boot <- .bootstrap_gliph_enrichment(trb,
                                      cluster_df[, c("cluster_id", "CDR3b")],
                                      n_iter = n_iter,
                                      seed   = cfg$seed %||% 42)
  if (is.null(boot) || nrow(boot) == 0) {
    log_message("  Bootstrap produced no rows; skipping enrichment table.")
    return(invisible(FALSE))
  }

  enrich <- boot |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarise(
      n_iter        = dplyr::n(),
      median_OR     = stats::median(OR, na.rm = TRUE),
      OR_q025       = stats::quantile(OR, 0.025, na.rm = TRUE),
      OR_q975       = stats::quantile(OR, 0.975, na.rm = TRUE),
      median_p      = stats::median(p,  na.rm = TRUE),
      n_viral_mean  = mean(n_viral),
      n_niu_mean    = mean(n_niu),
      .groups       = "drop"
    ) |>
    dplyr::mutate(
      direction = ifelse(median_OR > 1, "Viral_enriched", "NIU_enriched"),
      FDR       = stats::p.adjust(median_p, method = "BH")
    ) |>
    dplyr::arrange(FDR)

  # Attach motif text and member counts for downstream viz.
  motif_lookup <- cluster_df |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarise(motif      = dplyr::first(motif),
                     n_members  = dplyr::n(),
                     example_CDR3 = paste(utils::head(CDR3b, 3), collapse = ";"),
                     .groups = "drop")
  enrich <- dplyr::left_join(enrich, motif_lookup, by = "cluster_id")

  write.csv(enrich,
            file.path(out_tables, "gliph_enrichment_viral_vs_niu.csv"),
            row.names = FALSE)
  log_message("  Saved: gliph_enrichment_viral_vs_niu.csv (",
              nrow(enrich), " groups)")

  # ---- Persist for viz -----------------------------------------------------
  saveRDS(list(trb            = trb,
               clusters       = cluster_df,
               cluster_props  = cluster_props,
               enrich         = enrich,
               raw            = gliph_res),
          file.path(out_obj, "ImmGLIPHResults.rds"))
  log_message("  Saved: ImmGLIPHResults.rds")

  log_message("immGLIPH analysis complete.")
  invisible(TRUE)
}
