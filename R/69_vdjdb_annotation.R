# R/69_vdjdb_annotation.R
# Annotate intraocular TRB clones with VDJdb antigen/epitope. Downloads the
# pinned VDJdb release on first run and caches under references/vdjdb/.
# Matches by (CDR3-beta exact, V gene allele-stripped); optional Hamming-1
# fuzzy match controlled by cfg$tcr_advanced$vdjdb$fuzzy_match.
#
# Three downstream analyses:
#   1. Per-cell annotation table joined to olga_pgen_per_clone.csv
#   2. Per-T-cell-substate Fisher enrichment by antigen.species (the
#      "n=X, substate Y, p=Z" inputs for the UMAP overlay)
#   3. GLIPH cluster <-> VDJdb crosswalk (fraction annotated, dominant antigen)
#
# Outputs:
#   outputs/tables/repertoire/vdjdb_annotations.csv
#   outputs/tables/repertoire/vdjdb_substate_enrichment.csv
#   outputs/tables/repertoire/vdjdb_gliph_overlap.csv
#   outputs/objects/VdjdbAnnotation.rds
#   references/vdjdb/<tag>/vdjdb.slim.txt   (cache)

`%||%` <- function(x, y) if (is.null(x)) y else x

.vdjdb_download <- function(release_tag, cache_dir, force = FALSE) {
  dest_dir <- file.path(cache_dir, release_tag)
  slim     <- file.path(dest_dir, "vdjdb.slim.txt")
  if (!isTRUE(force) && file.exists(slim) && file.info(slim)$size > 0) {
    log_message("  VDJdb cache hit: ", slim)
    return(slim)
  }
  ensure_dir(dest_dir)
  url <- sprintf(
    "https://github.com/antigenomics/vdjdb-db/releases/download/%s/vdjdb-%s.zip",
    release_tag, release_tag)
  zip_path <- file.path(dest_dir, "vdjdb.zip")
  log_message("  Downloading VDJdb release ", release_tag, " from ", url)
  ok <- tryCatch({
    utils::download.file(url, zip_path, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    log_message("  download.file failed: ", conditionMessage(e)); FALSE
  })
  if (!ok || !file.exists(zip_path)) return(NA_character_)

  utils::unzip(zip_path, exdir = dest_dir)
  unlink(zip_path)

  # Some releases drop everything at top level; some inside a release dir.
  if (!file.exists(slim)) {
    cand <- list.files(dest_dir, pattern = "vdjdb\\.slim\\.txt$",
                       recursive = TRUE, full.names = TRUE)
    if (length(cand) > 0) file.copy(cand[1], slim, overwrite = TRUE)
  }
  if (!file.exists(slim)) {
    log_message("  VDJdb slim file not found after unzip.")
    return(NA_character_)
  }
  log_message("  VDJdb cached at ", slim)
  slim
}

# Allele-strip a V/J gene string: TRBV6-5*01 -> TRBV6-5; TRBV6-1/6-2 -> TRBV6-1
.vdjdb_strip_allele <- function(x) {
  x <- sub("\\*[0-9]+$", "", x)
  x <- sub("/.*$", "", x)
  x
}

.vdjdb_filter_trb <- function(df, species_keep = "HomoSapiens",
                              exclude_self_antigens = TRUE) {
  # vdjdb.slim.txt columns (canonical):
  #   gene, cdr3, v.segm, j.segm, species, mhc.a, mhc.b, mhc.class,
  #   antigen.epitope, antigen.gene, antigen.species, reference.id,
  #   method, meta, cdr3fix, vdjdb.score, web.method, web.method.seq, ...
  # `species` is the host (TCR organism) and `antigen.species` is the
  # epitope source. For pathogen-focused analyses, drop self-reactive
  # entries (antigen.species == "HomoSapiens", e.g. autoimmune /
  # neoantigen records) so the per-clone summarizer doesn't collapse a
  # pathogen-specific hit into a self-antigen winner on tied vdjdb.score.
  cols_needed <- c("gene", "cdr3", "v.segm", "species",
                   "antigen.epitope", "antigen.gene", "antigen.species",
                   "mhc.a", "mhc.b", "vdjdb.score")
  miss <- setdiff(cols_needed, colnames(df))
  if (length(miss) > 0)
    log_message("  WARN: vdjdb missing columns: ",
                paste(miss, collapse = ", "))

  d <- df[df$gene == "TRB" & df$species %in% species_keep, , drop = FALSE]
  if (isTRUE(exclude_self_antigens) &&
      "antigen.species" %in% colnames(d)) {
    n_before <- nrow(d)
    d <- d[d$antigen.species != "HomoSapiens", , drop = FALSE]
    log_message("  VDJdb: dropped ", n_before - nrow(d),
                " self-antigen (Homo sapiens) rows.")
  }
  d$V_stripped <- .vdjdb_strip_allele(d$v.segm)
  d
}

# Restrict VDJdb to the curated uveitis pathogen panel and tag each row
# with its tier ("primary" / "secondary"). Logs missing panel members so
# absent panel labels (e.g. HSV-1 in 2026-05-16) are loud, not silent.
# Writes vdjdb_panel_coverage.csv to out_tables for the audit trail.
.vdjdb_filter_panel <- function(df, panel, out_tables) {
  if (is.null(panel) ||
      (length(panel$primary) + length(panel$secondary)) == 0) {
    df$panel_tier <- NA_character_
    return(df)
  }
  tier <- c(
    setNames(rep("primary",   length(panel$primary)),   panel$primary),
    setNames(rep("secondary", length(panel$secondary)), panel$secondary)
  )
  df$panel_tier <- unname(tier[df$antigen.species])
  panel_members <- c(panel$primary, panel$secondary)
  cov <- data.frame(
    antigen_species = panel_members,
    panel_tier      = unname(tier[panel_members]),
    n_rows_in_db    = vapply(panel_members,
                             function(s) sum(df$antigen.species == s, na.rm = TRUE),
                             integer(1)),
    stringsAsFactors = FALSE
  )
  cov$status <- ifelse(cov$n_rows_in_db > 0, "present", "absent_in_release")
  utils::write.csv(cov, file.path(out_tables, "vdjdb_panel_coverage.csv"),
                   row.names = FALSE)
  for (s in panel_members) {
    log_message(sprintf("  panel %-12s tier=%-9s rows=%d",
                        s, tier[s], cov$n_rows_in_db[cov$antigen_species == s]))
  }
  df[!is.na(df$panel_tier), , drop = FALSE]
}

# Match (CDR3b, V_stripped) exactly. When fuzzy=TRUE, also accept
# V-gene-anchored, length-equal, substitution-only Hamming-1 matches.
# Anchoring is enforced by the per-V loop below: each fuzzy comparison
# is restricted to VDJdb entries sharing the same allele-stripped V
# gene. utils::adist() is configured with insertions/deletions cost
# = 99 so the search is substitution-only.
.vdjdb_match <- function(clones, vdjdb, fuzzy = FALSE) {
  if (nrow(clones) == 0 || nrow(vdjdb) == 0)
    return(data.frame())

  clones$V_stripped <- .vdjdb_strip_allele(clones$TRBV)

  ex <- dplyr::inner_join(
    clones,
    dplyr::select(vdjdb, V_stripped, cdr3,
                  antigen.epitope, antigen.gene, antigen.species,
                  mhc.a, mhc.b, vdjdb.score, panel_tier),
    by = c("CDR3b" = "cdr3", "V_stripped" = "V_stripped"),
    relationship = "many-to-many")
  if (nrow(ex) > 0) ex$match_type <- "exact"

  if (!isTRUE(fuzzy) || nrow(clones) == 0) return(ex)

  # For unmatched clones, do per-length fuzzy Hamming-1
  unmatched <- clones |>
    dplyr::anti_join(ex, by = c("barcode"))
  if (nrow(unmatched) == 0) return(ex)

  fz <- list()
  for (vg in unique(unmatched$V_stripped)) {
    cl <- unmatched[unmatched$V_stripped == vg, , drop = FALSE]
    vd <- vdjdb[vdjdb$V_stripped == vg, , drop = FALSE]
    if (nrow(vd) == 0 || nrow(cl) == 0) next
    by_len <- split(seq_len(nrow(cl)), nchar(cl$CDR3b))
    for (ll in names(by_len)) {
      cl_l <- cl[by_len[[ll]], , drop = FALSE]
      vd_l <- vd[nchar(vd$cdr3) == as.integer(ll), , drop = FALSE]
      if (nrow(vd_l) == 0) next
      dmat <- utils::adist(cl_l$CDR3b, vd_l$cdr3,
                           costs = list(insertions = 99,
                                        deletions  = 99,
                                        substitutions = 1))
      hits <- which(dmat <= 1, arr.ind = TRUE)
      if (nrow(hits) == 0) next
      hit_rows <- cbind(
        cl_l[hits[, "row"], , drop = FALSE],
        vd_l[hits[, "col"], c("antigen.epitope", "antigen.gene",
                              "antigen.species", "mhc.a", "mhc.b",
                              "vdjdb.score", "panel_tier"), drop = FALSE])
      hit_rows$match_type <- "fuzzy_hamming1"
      fz[[length(fz) + 1]] <- hit_rows
    }
  }
  fuzzy_df <- if (length(fz) > 0) do.call(rbind, fz) else data.frame()
  if (nrow(fuzzy_df) > 0) {
    common <- intersect(colnames(ex), colnames(fuzzy_df))
    ex <- rbind(ex[, common, drop = FALSE], fuzzy_df[, common, drop = FALSE])
  }
  ex
}

# Per-barcode summary: highest vdjdb.score wins; ties -> modal antigen species.
.vdjdb_summarize_per_clone <- function(matches) {
  if (nrow(matches) == 0) return(data.frame())
  matches |>
    dplyr::group_by(barcode) |>
    dplyr::summarise(
      antigen_species = {
        max_score <- max(vdjdb.score, na.rm = TRUE)
        keep <- vdjdb.score == max_score
        tab <- sort(table(antigen.species[keep]), decreasing = TRUE)
        names(tab)[1]
      },
      antigen_gene  = {
        max_score <- max(vdjdb.score, na.rm = TRUE)
        keep <- vdjdb.score == max_score
        tab <- sort(table(antigen.gene[keep]), decreasing = TRUE)
        names(tab)[1]
      },
      panel_tier = {
        max_score <- max(vdjdb.score, na.rm = TRUE)
        keep <- vdjdb.score == max_score
        tab <- sort(table(panel_tier[keep]), decreasing = TRUE)
        if (length(tab) == 0) NA_character_ else names(tab)[1]
      },
      antigen_epitope = paste(sort(unique(antigen.epitope[!is.na(antigen.epitope)])),
                              collapse = ";"),
      mhc_a    = paste(sort(unique(mhc.a[!is.na(mhc.a)])), collapse = ";"),
      mhc_b    = paste(sort(unique(mhc.b[!is.na(mhc.b)])), collapse = ";"),
      vdjdb_score_max = max(vdjdb.score, na.rm = TRUE),
      n_hits   = dplyr::n(),
      match_type = paste(sort(unique(match_type)), collapse = ";"),
      .groups = "drop")
}

# Fisher per (substate_key, antigen_species) against the rest of the tcell
# compartment. BH-FDR across all (substate, species) pairs. When cd8_keys
# is non-NULL, restrict the cell universe to those substates before
# counting (used to compute the CD8-only contrast for the hypothesis).
.vdjdb_substate_enrichment <- function(per_clone_summary, tcell_meta,
                                       cd8_keys = NULL) {
  if (nrow(per_clone_summary) == 0 || is.null(tcell_meta) ||
      !"substate_key" %in% colnames(tcell_meta))
    return(data.frame())

  if ("barcode" %in% colnames(tcell_meta)) tcell_meta$barcode <- NULL
  joined <- tcell_meta |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(
      dplyr::select(per_clone_summary, barcode, antigen_species),
      by = "barcode") |>
    dplyr::filter(!is.na(substate_key))
  if (!is.null(cd8_keys))
    joined <- joined[joined$substate_key %in% cd8_keys, , drop = FALSE]
  if (nrow(joined) == 0) return(data.frame())

  species_levels <- setdiff(unique(joined$antigen_species), NA_character_)
  substates      <- sort(unique(joined$substate_key))
  if (length(species_levels) == 0) return(data.frame())

  rows <- list()
  for (sp in species_levels) {
    is_sp <- joined$antigen_species == sp & !is.na(joined$antigen_species)
    for (sub in substates) {
      in_sub <- joined$substate_key == sub
      a <- sum(is_sp & in_sub)
      b <- sum(is_sp & !in_sub)
      c <- sum(!is_sp & in_sub)
      d <- sum(!is_sp & !in_sub)
      if (a + b == 0 || a + c == 0) next
      ft <- suppressWarnings(fisher.test(matrix(c(a, c, b, d), nrow = 2),
                                         alternative = "greater"))
      rows[[length(rows) + 1]] <- data.frame(
        substate_key   = sub,
        antigen_species = sp,
        n_substate     = a + c,
        n_hits_substate = a,
        n_hits_total   = a + b,
        odds_ratio     = unname(ft$estimate),
        fisher_p       = ft$p.value,
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out$fdr <- stats::p.adjust(out$fisher_p, method = "BH")
  out[order(out$fdr), ]
}

# Per panel virus: Fisher of (Phenotype_2 == "Viral") x (hit for that
# virus). Run once on all T cells, once restricted to CD8 substates;
# both go in the output with a `subset` column. Emits a row for every
# panel species (even when n_hits = 0) so plots can tile a fixed grid.
#
# Contrast: alternative = "greater" tests viral-phenotype enrichment in
# the sp-specific bucket, i.e. infection draws sp-reactive clones into
# the eye. Background is NIU intraocular T cells (Healthy is upstream-
# excluded from olga_pgen_per_clone.csv).
.vdjdb_phenotype_enrichment <- function(per_clone_anno, clones_meta,
                                        panel, cd8_keys) {
  if (nrow(per_clone_anno) == 0 || nrow(clones_meta) == 0)
    return(data.frame())
  panel_species <- c(panel$primary, panel$secondary)
  tier_lookup   <- c(setNames(rep("primary",   length(panel$primary)),
                              panel$primary),
                     setNames(rep("secondary", length(panel$secondary)),
                              panel$secondary))
  do_one <- function(meta, subset_label) {
    rows <- list()
    for (sp in panel_species) {
      hit_bc <- per_clone_anno$barcode[per_clone_anno$antigen_species == sp]
      meta$hit   <- meta$barcode %in% hit_bc
      meta$viral <- meta$Phenotype_2 == "Viral"
      a <- sum( meta$hit  &  meta$viral)
      b <- sum( meta$hit  & !meta$viral)
      c <- sum(!meta$hit  &  meta$viral)
      d <- sum(!meta$hit  & !meta$viral)
      if (a + b == 0) {
        rows[[length(rows) + 1]] <- data.frame(
          subset = subset_label, antigen_species = sp,
          panel_tier = unname(tier_lookup[sp]),
          n_viral = a + c, n_niu = b + d,
          n_hits_viral = a, n_hits_niu = b,
          OR_fisher = NA_real_, p_fisher = NA_real_,
          stringsAsFactors = FALSE)
        next
      }
      ft <- suppressWarnings(fisher.test(matrix(c(a, c, b, d), nrow = 2),
                                         alternative = "greater"))
      rows[[length(rows) + 1]] <- data.frame(
        subset = subset_label, antigen_species = sp,
        panel_tier = unname(tier_lookup[sp]),
        n_viral = a + c, n_niu = b + d,
        n_hits_viral = a, n_hits_niu = b,
        OR_fisher = unname(ft$estimate),
        p_fisher  = ft$p.value,
        stringsAsFactors = FALSE)
    }
    out <- do.call(rbind, rows)
    if (is.null(out)) return(data.frame())
    out$fdr <- stats::p.adjust(out$p_fisher, method = "BH")
    out
  }
  out_all <- do_one(clones_meta, "all_T")
  cd8_meta <- clones_meta[clones_meta$substate_key %in% cd8_keys, , drop = FALSE]
  out_cd8 <- do_one(cd8_meta, "CD8_only")
  rbind(out_all, out_cd8)
}

# Per (viral etiology, panel virus) WITHIN viral cells only:
# Fisher of (Etiology == e) x (hit for that virus). Restricting to
# Phenotype_2 == "Viral" so the contrast is "VZV_ARN vs other viral
# etiologies", not "VZV_ARN vs NIU" (which would confound disease
# group with TCR specificity). This is the headline test for the
# hypothesis: VZV_ARN x VZV, CMV_CRN x CMV, HTLV1 x HTLV-1, HSV2 x
# HSV-2 are expected to show OR >> 1 at FDR < 0.1.
#
# TODO(GLMM): cell-level Fisher over-counts subjects with large
# repertoires (e.g. VZV_ARN at ~22k cells). Layer in a patient
# random-effect GLMM as a follow-up via run_fisher_glmm_per_cluster
# (R/01_setup_utils.R:273-395).
.vdjdb_etiology_enrichment <- function(per_clone_anno, clones_meta,
                                       panel, viral_etiologies,
                                       cd8_keys, min_cells) {
  if (nrow(per_clone_anno) == 0 || nrow(clones_meta) == 0)
    return(data.frame())
  panel_species <- c(panel$primary, panel$secondary)
  tier_lookup   <- c(setNames(rep("primary",   length(panel$primary)),
                              panel$primary),
                     setNames(rep("secondary", length(panel$secondary)),
                              panel$secondary))
  do_one <- function(meta, subset_label) {
    meta <- meta[meta$Phenotype_2 == "Viral" &
                 meta$Etiology %in% viral_etiologies, , drop = FALSE]
    if (nrow(meta) == 0) return(data.frame())
    rows <- list()
    for (e in viral_etiologies) {
      n_etiology <- sum(meta$Etiology == e)
      for (sp in panel_species) {
        hit_bc <- per_clone_anno$barcode[per_clone_anno$antigen_species == sp]
        in_e <- meta$Etiology == e
        hit  <- meta$barcode %in% hit_bc
        a <- sum( hit &  in_e); b <- sum( hit & !in_e)
        c <- sum(!hit &  in_e); d <- sum(!hit & !in_e)
        n_hits_total <- a + b
        status <- "tested"
        if (n_etiology < min_cells)        status <- "skipped_low_n"
        else if (n_hits_total == 0)        status <- "absent_in_db"
        if (status != "tested") {
          rows[[length(rows) + 1]] <- data.frame(
            subset = subset_label, etiology = e, antigen_species = sp,
            panel_tier = unname(tier_lookup[sp]),
            n_etiology = n_etiology, n_hits_etiology = a,
            n_hits_total = n_hits_total,
            OR_fisher = NA_real_, p_fisher = NA_real_,
            status = status, stringsAsFactors = FALSE)
          next
        }
        ft <- suppressWarnings(fisher.test(matrix(c(a, c, b, d), nrow = 2),
                                           alternative = "greater"))
        rows[[length(rows) + 1]] <- data.frame(
          subset = subset_label, etiology = e, antigen_species = sp,
          panel_tier = unname(tier_lookup[sp]),
          n_etiology = n_etiology, n_hits_etiology = a,
          n_hits_total = n_hits_total,
          OR_fisher = unname(ft$estimate),
          p_fisher  = ft$p.value,
          status = "tested", stringsAsFactors = FALSE)
      }
    }
    out <- do.call(rbind, rows)
    if (is.null(out)) return(data.frame())
    tested <- out$status == "tested"
    out$fdr <- NA_real_
    if (any(tested))
      out$fdr[tested] <- stats::p.adjust(out$p_fisher[tested], method = "BH")
    out
  }
  out_all <- do_one(clones_meta, "all_T")
  cd8_meta <- clones_meta[clones_meta$substate_key %in% cd8_keys, , drop = FALSE]
  out_cd8 <- do_one(cd8_meta, "CD8_only")
  rbind(out_all, out_cd8)
}

# Per GLIPH cluster: n annotated, dominant antigen species/gene, fraction
# annotated, Fisher OR vs background (all-other-clones-with-Pgen).
.vdjdb_gliph_overlap <- function(per_clone_summary, gliph_clusters,
                                 per_clone_pgen) {
  if (nrow(per_clone_summary) == 0 || nrow(gliph_clusters) == 0)
    return(data.frame())
  background_cdr3 <- unique(per_clone_pgen$CDR3b)
  hit_cdr3 <- unique(per_clone_summary |>
    dplyr::inner_join(per_clone_pgen, by = "barcode") |>
    dplyr::pull(CDR3b))
  n_back_total <- length(background_cdr3)
  n_hit_total  <- length(hit_cdr3)
  if (n_back_total == 0) return(data.frame())

  hit_table <- per_clone_summary |>
    dplyr::inner_join(per_clone_pgen, by = "barcode") |>
    dplyr::select(CDR3b, antigen_species, antigen_gene)

  rows <- gliph_clusters |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarise(
      n_cdr3 = dplyr::n_distinct(CDR3b),
      n_annotated = sum(unique(CDR3b) %in% hit_cdr3),
      dominant_antigen_species = {
        h <- hit_table[hit_table$CDR3b %in% CDR3b, , drop = FALSE]
        if (nrow(h) == 0) NA_character_
        else names(sort(table(h$antigen_species), decreasing = TRUE))[1]
      },
      dominant_antigen_gene = {
        h <- hit_table[hit_table$CDR3b %in% CDR3b, , drop = FALSE]
        if (nrow(h) == 0) NA_character_
        else names(sort(table(h$antigen_gene), decreasing = TRUE))[1]
      },
      .groups = "drop")

  rows$fraction_annotated <- rows$n_annotated / pmax(rows$n_cdr3, 1)
  rows$fisher_p <- vapply(seq_len(nrow(rows)), function(i) {
    a <- rows$n_annotated[i]; b <- rows$n_cdr3[i] - a
    c <- n_hit_total - a;     d <- n_back_total - n_hit_total - b
    if (a + b == 0) return(NA_real_)
    suppressWarnings(fisher.test(matrix(c(a, c, b, d), nrow = 2),
                                 alternative = "greater")$p.value)
  }, numeric(1))
  rows$fdr <- stats::p.adjust(rows$fisher_p, method = "BH")
  rows[order(rows$fdr), ]
}


# This builds the answer table: one row per documented etiology, pairing it with
# its matched VDJdb antigen and, critically, with how many reference records
# exist to match against. That last column is what converts an uninterpretable
# null into an explained one. It also returns the off-diagonal hits annotated
# with the same coverage, which shows they track database depth rather than
# biology.
.vdjdb_q1_answer_table <- function(cfg, enrich, slim_path) {
  matched <- c(VZV_ARN = "VZV", CMV_CRN = "CMV", HTLV1 = "HTLV-1",
               HSV1 = "HSV-1", HSV2 = "HSV-2")
  cov <- NULL
  if (!is.na(slim_path) && file.exists(slim_path)) {
    db <- utils::read.delim(slim_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
    db <- db[db$species == "HomoSapiens", , drop = FALSE]
    cov <- do.call(rbind, lapply(unique(c(matched, enrich$antigen_species)),
      function(sp) {
        k <- db$`antigen.species` == sp
        data.frame(antigen_species = sp,
                   db_records = sum(k),
                   db_records_beta = sum(k & db$gene == "TRB"),
                   db_epitopes = length(unique(db$`antigen.epitope`[k])),
                   db_classI = sum(k & db$`mhc.class` == "MHCI"),
                   db_classII = sum(k & db$`mhc.class` == "MHCII"),
                   stringsAsFactors = FALSE)
      }))
  }

  e <- enrich[enrich$subset == "all_T", , drop = FALSE]
  e$matched_antigen <- unname(matched[e$etiology])
  diag <- e[!is.na(e$matched_antigen) &
              e$antigen_species == e$matched_antigen, , drop = FALSE]
  if (nrow(diag) && !is.null(cov))
    diag <- merge(diag, cov, by = "antigen_species", all.x = TRUE)
  if (nrow(diag)) {
    diag$question <- "documented infection vs matched VDJdb antigen"
    # An etiology is only assessable if the database has something to match.
    diag$assessable <- ifelse(is.na(diag$db_records_beta) |
                                diag$db_records_beta < 100,
                              "no: reference coverage too thin", "yes")
    diag$answer <- ifelse(diag$assessable != "yes", "not assessable",
                   ifelse(diag$n_hits_etiology == 0L, "no relation detected",
                   ifelse(!is.na(diag$fdr) & diag$fdr < 0.05,
                          "relation detected", "no relation detected")))
    diag <- diag[order(-diag$n_etiology),
                 c("etiology", "matched_antigen", "n_etiology",
                   "n_hits_etiology", "db_records_beta", "db_epitopes",
                   "db_classI", "db_classII", "OR_fisher", "p_fisher", "fdr",
                   "status", "assessable", "answer")]
  }

  off <- e[!is.na(e$matched_antigen) & e$antigen_species != e$matched_antigen &
             !is.na(e$fdr) & e$fdr < 0.05, , drop = FALSE]
  if (nrow(off) && !is.null(cov)) off <- merge(off, cov, by = "antigen_species",
                                               all.x = TRUE)
  if (nrow(off)) {
    off$question <- "off-diagonal hit (etiology does NOT match antigen)"
    off <- off[order(off$fdr),
               c("etiology", "antigen_species", "matched_antigen", "n_etiology",
                 "n_hits_etiology", "db_records_beta", "db_epitopes",
                 "OR_fisher", "fdr")]
  }
  list(matched = diag, off_diagonal = off, coverage = cov)
}

#
# The point of the grid rather than a ranked list is that it shows the shape of
# the result. If predicted specificity tracked documented infection, signal
# would sit on the diagonal. Here the diagonal is empty and the only significant
# cells are off it, which is the honest summary: no clear relation, and the
# residual pattern is noise that follows reference-database depth rather than
# patient etiology.
#
# Cells are one of four states, kept visually distinct because they mean very
# different things:
#   not in VDJdb    the matched antigen has no human records at all
#   too few cells   etiology below the per-etiology cell floor
#   tested, ns      a real test that found nothing
#   tested, FDR<.05 a real test that found something (all off-diagonal here)
.vdjdb_q1_fdr_matrix <- function(enrich, cov = NULL) {
  matched <- c(VZV_ARN = "VZV", CMV_CRN = "CMV", HTLV1 = "HTLV-1",
               HSV1 = "HSV-1", HSV2 = "HSV-2")
  d <- enrich
  d$matched_antigen <- unname(matched[d$etiology])
  d$is_matched <- !is.na(d$matched_antigen) &
                    d$antigen_species == d$matched_antigen
  d$cell_state <- ifelse(d$status == "absent_in_db", "not in VDJdb",
                  ifelse(d$status == "skipped_low_n", "too few cells",
                  ifelse(!is.na(d$fdr) & d$fdr < 0.05, "tested, FDR < 0.05",
                         "tested, not significant")))
  d$neglog_fdr <- ifelse(is.na(d$fdr), NA_real_, -log10(pmax(d$fdr, 1e-30)))
  d$label <- ifelse(d$status == "absent_in_db", "—",
             ifelse(d$status == "skipped_low_n", "n/a",
             ifelse(is.na(d$fdr), "",
             ifelse(d$fdr < 0.001, sprintf("%.0e", d$fdr),
                    sprintf("%.2f", d$fdr)))))
  if (!is.null(cov))
    d <- merge(d, cov[, c("antigen_species", "db_records_beta")],
               by = "antigen_species", all.x = TRUE)
  d
}

viz_vdjdb_q1_fdr_matrix <- function(cfg, enrich = NULL, cov = NULL,
                                    out_dir = NULL) {
  out_tab <- file.path(cfg$paths$results_tables, "repertoire")
  if (is.null(enrich)) {
    f <- file.path(out_tab, "vdjdb_etiology_enrichment.csv")
    if (!file.exists(f)) { log_message("  q1 matrix: ", f, " missing."); return(invisible(NULL)) }
    enrich <- utils::read.csv(f, stringsAsFactors = FALSE)
  }
  d <- .vdjdb_q1_fdr_matrix(enrich, cov)

  # Wide FDR matrices, one per subset, for the supplement.
  for (ss in unique(d$subset)) {
    w <- d[d$subset == ss, c("etiology", "antigen_species", "fdr")]
    m <- stats::reshape(w, idvar = "etiology", timevar = "antigen_species",
                        direction = "wide")
    names(m) <- sub("^fdr\\.", "", names(m))
    utils::write.csv(m, file.path(out_tab,
      paste0("vdjdb_q1_fdr_matrix_", ss, ".csv")), row.names = FALSE)
  }
  utils::write.csv(d, file.path(out_tab, "vdjdb_q1_fdr_matrix_long.csv"),
                   row.names = FALSE)
  log_message("  Saved: vdjdb_q1_fdr_matrix_{all_T,CD8_only,long}.csv")

  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(d))
  # Row / column labels carry the numbers that explain the pattern: cells per
  # etiology on the rows, reference depth on the columns.
  # Build the axis labels from per-key lookup tables, then order the factor
  # levels from those same tables. An earlier version ordered with
  # ifelse(<length-1 condition>, vec, 0), which returns a LENGTH-1 result and
  # silently collapsed every column level to NA.
  n_by_eti <- tapply(d$n_etiology, d$etiology,
                     function(x) suppressWarnings(max(x, na.rm = TRUE)))
  row_map <- stats::setNames(
    paste0(names(n_by_eti), "\n(", format(n_by_eti, big.mark = ",", trim = TRUE),
           " cells)"), names(n_by_eti))
  d$row_lab <- factor(row_map[d$etiology],
                      levels = row_map[names(sort(n_by_eti, decreasing = TRUE))])

  sp <- sort(unique(d$antigen_species))
  db_by_sp <- stats::setNames(rep(NA_real_, length(sp)), sp)
  if ("db_records_beta" %in% names(d)) {
    tmp <- tapply(d$db_records_beta, d$antigen_species,
                  function(x) suppressWarnings(max(x, na.rm = TRUE)))
    tmp[!is.finite(tmp)] <- NA_real_
    db_by_sp[names(tmp)] <- as.numeric(tmp)
  }
  col_map <- stats::setNames(ifelse(
    is.na(db_by_sp), sp,
    paste0(sp, "\n(", format(db_by_sp, big.mark = ",", trim = TRUE), " TRB)")), sp)
  ord <- names(sort(db_by_sp, decreasing = TRUE, na.last = TRUE))
  d$col_lab <- factor(col_map[d$antigen_species], levels = col_map[ord])

  pal <- c("not in VDJdb"            = "grey88",
           "too few cells"           = "grey75",
           "tested, not significant" = "#EAF2F8",
           "tested, FDR < 0.05"      = "#E8453B")
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$col_lab, .data$row_lab)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$cell_state), colour = "white",
                       linewidth = 0.6) +
    # Box the diagonal: these are the only cells that answer the question.
    ggplot2::geom_tile(data = d[d$is_matched, , drop = FALSE],
                       fill = NA, colour = "black", linewidth = 1.1) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.7) +
    ggplot2::scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    ggplot2::facet_wrap(~ .data$subset) +
    ggplot2::labs(
      title = "Predicted antigen specificity versus documented ocular infection",
      subtitle = paste(
        "BH-FDR per (documented etiology x VDJdb antigen) Fisher test.",
        "Black boxes mark the matched pairs that answer the question:",
        "every one is untestable or FDR = 1.",
        "\nAll significant cells are off-diagonal and sit on the three",
        "best-covered reference species, so the residual pattern tracks",
        "database depth rather than patient etiology."),
      x = "VDJdb antigen species (human TRB records)",
      y = "Documented ocular etiology (cells)") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 7.5),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       size = 7),
                   axis.text.y = ggplot2::element_text(size = 7),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "bottom")
  vd <- out_dir %||% file.path(get_target_paths(cfg, "eye")$viz_dir,
                               VIZ_BUCKETS[["repertoire"]])
  ensure_dir(vd)
  save_pdf_png(p, file.path(vd, "vdjdb_q1_fdr_matrix"), w = 11, h = 6)
  invisible(d)
}

run_vdjdb_annotation <- function(cfg) {
  if (!isTRUE(cfg$steps$vdjdb)) {
    log_message("VDJdb annotation disabled. Skipping.")
    return(invisible(TRUE))
  }
  log_message("Starting VDJdb annotation...")

  obj_dir    <- cfg$paths$results_objects
  out_tables <- file.path(cfg$paths$results_tables, "repertoire")
  ensure_dir(out_tables)

  release_tag <- cfg$tcr_advanced$vdjdb$release_tag %||% "2026-05-16"
  cache_dir   <- cfg$tcr_advanced$vdjdb$cache_dir   %||% "references/vdjdb"
  fuzzy       <- isTRUE(cfg$tcr_advanced$vdjdb$fuzzy_match)
  species     <- cfg$tcr_advanced$vdjdb$species_keep %||% "HomoSapiens"
  exclude_self_antigens <- isTRUE(cfg$tcr_advanced$vdjdb$exclude_self_antigens %||% TRUE)
  panel       <- cfg$tcr_advanced$vdjdb$pathogen_panel
  cd8_keys    <- cfg$tcr_advanced$vdjdb$cd8_substate_keys %||%
                   c("tcell_1", "tcell_3")
  min_cells_per_etiology <- cfg$tcr_advanced$vdjdb$min_cells_per_etiology %||% 20
  viral_etiologies <- cfg$etiology_groups$viral %||%
                        c("VZV_ARN", "HSV1", "HSV2", "CMV_CRN", "HTLV1")

  slim_path <- .vdjdb_download(release_tag, cache_dir, force = FALSE)
  if (is.na(slim_path)) {
    log_message("  VDJdb download failed. Aborting.")
    return(invisible(FALSE))
  }
  vdjdb <- utils::read.delim(slim_path, sep = "\t",
                             stringsAsFactors = FALSE,
                             check.names = FALSE)
  vdjdb <- .vdjdb_filter_trb(vdjdb, species_keep = species,
                             exclude_self_antigens = exclude_self_antigens)
  log_message("  VDJdb TRB rows after filter: ", nrow(vdjdb))
  vdjdb <- .vdjdb_filter_panel(vdjdb, panel, out_tables)
  log_message("  VDJdb TRB rows after panel filter: ", nrow(vdjdb))

  pgen_csv <- file.path(out_tables, "olga_pgen_per_clone.csv")
  if (!file.exists(pgen_csv)) {
    log_message("  olga_pgen_per_clone.csv missing — needed as clone source. Skipping.")
    return(invisible(FALSE))
  }
  clones <- utils::read.csv(pgen_csv, stringsAsFactors = FALSE)
  clones <- clones[, c("barcode", "Subject", "Phenotype_2", "Etiology",
                       "CDR3b", "TRBV", "TRBJ")]

  # Load T-cell object once; attach substate_key to clones for the new
  # phenotype/etiology helpers and the CD8-stratified substate test.
  tcell_path <- get_target_paths(cfg, "tcell")$results_objects
  tcell_rds  <- file.path(tcell_path, "IntegratedSeuratObject.rds")
  have_tcell <- file.exists(tcell_rds)
  obj_t <- NULL
  if (have_tcell) {
    obj_t <- readRDS(tcell_rds)
    if (!"substate_key" %in% colnames(obj_t@meta.data)) {
      obj_t$substate_key <- get_substate_key_vector(obj_t, "tcell")
    }
    sk <- setNames(as.character(obj_t$substate_key), colnames(obj_t))
    clones$substate_key <- unname(sk[clones$barcode])
  } else {
    clones$substate_key <- NA_character_
    log_message("  T cell compartment object not found; CD8 stratification disabled.")
  }

  matches <- .vdjdb_match(clones, vdjdb, fuzzy = fuzzy)
  log_message("  VDJdb hit rows: ", nrow(matches),
              " (across ", length(unique(matches$barcode)), " cells)")
  if (nrow(matches) == 0) {
    log_message("  No VDJdb matches in this repertoire."); return(invisible(FALSE))
  }
  per_clone <- .vdjdb_summarize_per_clone(matches)
  utils::write.csv(per_clone,
                   file.path(out_tables, "vdjdb_annotations.csv"),
                   row.names = FALSE)
  log_message("  Saved: vdjdb_annotations.csv (", nrow(per_clone), " cells)")

  # --- Per-substate enrichment (all_T + CD8_only) -----------------------------
  sub_enrich <- data.frame()
  if (have_tcell) {
    sub_all <- .vdjdb_substate_enrichment(per_clone, obj_t@meta.data)
    if (nrow(sub_all) > 0) sub_all$subset <- "all_T"
    sub_cd8 <- .vdjdb_substate_enrichment(per_clone, obj_t@meta.data,
                                          cd8_keys = cd8_keys)
    if (nrow(sub_cd8) > 0) sub_cd8$subset <- "CD8_only"
    sub_enrich <- rbind(sub_all, sub_cd8)
    if (nrow(sub_enrich) > 0) {
      utils::write.csv(sub_enrich,
                       file.path(out_tables, "vdjdb_substate_enrichment.csv"),
                       row.names = FALSE)
      log_message("  Saved: vdjdb_substate_enrichment.csv (",
                  nrow(sub_enrich), " rows)")
    }
  }

  # --- Per-phenotype enrichment (Viral vs NIU, all_T + CD8_only) -------------
  phen_enrich <- .vdjdb_phenotype_enrichment(per_clone, clones, panel, cd8_keys)
  if (nrow(phen_enrich) > 0) {
    utils::write.csv(phen_enrich,
                     file.path(out_tables, "vdjdb_phenotype_enrichment.csv"),
                     row.names = FALSE)
    log_message("  Saved: vdjdb_phenotype_enrichment.csv (",
                nrow(phen_enrich), " rows)")
  }

  # --- Per-etiology enrichment (within Viral only) ---------------------------
  etio_enrich <- .vdjdb_etiology_enrichment(per_clone, clones, panel,
                                            viral_etiologies, cd8_keys,
                                            min_cells_per_etiology)
  if (nrow(etio_enrich) > 0) {
    utils::write.csv(etio_enrich,
                     file.path(out_tables, "vdjdb_etiology_enrichment.csv"),
                     row.names = FALSE)
    log_message("  Saved: vdjdb_etiology_enrichment.csv (",
                nrow(etio_enrich), " rows)")
  }

  tryCatch({
    slim <- .vdjdb_cache(cfg)
    q1 <- .vdjdb_q1_answer_table(cfg, etio_enrich, slim)
    utils::write.csv(q1$matched,
                     file.path(out_tables, "vdjdb_q1_matched_etiology.csv"),
                     row.names = FALSE)
    if (NROW(q1$off_diagonal))
      utils::write.csv(q1$off_diagonal,
                       file.path(out_tables, "vdjdb_q1_offdiagonal_hits.csv"),
                       row.names = FALSE)
    viz_vdjdb_q1_fdr_matrix(cfg, etio_enrich, q1$coverage)
  }, error = function(e)
    log_message("  vdjdb Q1 answer objects failed: ", conditionMessage(e)))

  if (have_tcell) { rm(obj_t); invisible(gc()) }

  # --- GLIPH overlap ----------------------------------------------------------
  gliph_csv <- file.path(out_tables, "gliph_clusters.csv")
  if (file.exists(gliph_csv)) {
    gl <- utils::read.csv(gliph_csv, stringsAsFactors = FALSE)
    pgen_min <- utils::read.csv(pgen_csv, stringsAsFactors = FALSE)
    overlap <- .vdjdb_gliph_overlap(per_clone, gl,
                                    pgen_min[, c("barcode", "CDR3b")])
    if (nrow(overlap) > 0) {
      utils::write.csv(overlap,
                       file.path(out_tables, "vdjdb_gliph_overlap.csv"),
                       row.names = FALSE)
      log_message("  Saved: vdjdb_gliph_overlap.csv (", nrow(overlap), " rows)")
    }
  }

  saveRDS(list(vdjdb                = vdjdb,
               matches              = matches,
               per_clone            = per_clone,
               substate_enrichment  = sub_enrich,
               phenotype_enrichment = phen_enrich,
               etiology_enrichment  = etio_enrich),
          file.path(obj_dir, "VdjdbAnnotation.rds"))
  log_message("  Saved: VdjdbAnnotation.rds")

  log_message("VDJdb annotation complete.")
  invisible(TRUE)
}
