# R/67_gliph_tcrdist_correlation.R
# Concordance between GLIPH convergence groups and tcrdist / clusTCR
# clusterings on the intraocular TRB repertoire.
#
# Two complementary tests:
#   1. Assignment-level: ARI / AMI / NMI between GLIPH primary cluster and
#      clusTCR cluster (per-cell, exploded multi-membership = highest fisher
#      score).
#   2. Distance-level: per GLIPH cluster, mean tcrdist within the cluster vs
#      a Subject-stratified background pool. 1000-permutation null + BH-FDR.
# Plus: bipartite Fisher of (GLIPH cluster x clusTCR cluster) contingency
# counts → log2OR + BH q per cell.
#
# Inputs (all written by 64/65 + saved pw_beta from 65):
#   outputs/objects/ImmGLIPHResults.rds        — gliph cluster_df + props
#   outputs/objects/ImmLynxTcrdistResults.rds  — sce + barcodes
#   outputs/objects/tcrdist_pw_beta.rds        — barcode-indexed distance matrix
#   outputs/tables/repertoire/gliph_clusters.csv
#   outputs/tables/repertoire/clustcr_clusters.csv
# Outputs:
#   outputs/tables/repertoire/gliph_tcrdist_concordance.csv
#   outputs/tables/repertoire/gliph_clustcr_overlap.csv
#   outputs/objects/GliphTcrdistJoint.rds

`%||%` <- function(x, y) if (is.null(x)) y else x

# Build a per-barcode joint table linking GLIPH (multi-membership preserved
# as a list-col) to clusTCR (single membership). Picks a single "primary"
# GLIPH cluster per CDR3 (highest fisher.score) for ARI/AMI computations.
.gtc_build_join <- function(gliph_clusters, gliph_props, clustcr_df,
                            trb_frame) {
  # gliph_clusters is keyed by (cluster_id, CDR3b, TRBV). The same (CDR3b,
  # TRBV) can belong to >1 cluster. Pick the highest-fisher.score cluster
  # per CDR3 as the primary, and keep the full list as audit info.
  score_lookup <- setNames(as.numeric(gliph_props$fisher.score),
                           gliph_props$cluster_id)
  gliph_clusters$.score <- score_lookup[gliph_clusters$cluster_id]
  gliph_clusters$.score[is.na(gliph_clusters$.score)] <- 0

  by_cdr3 <- gliph_clusters |>
    dplyr::group_by(CDR3b) |>
    dplyr::summarise(
      gliph_primary_cluster = cluster_id[which.max(.score)][1],
      gliph_cluster_ids     = paste(sort(unique(cluster_id)), collapse = ";"),
      n_gliph_memberships   = dplyr::n_distinct(cluster_id),
      .groups = "drop"
    )

  # trb_frame has one row per cell (barcode, Subject, Phenotype_2, CDR3b)
  joint <- trb_frame |>
    dplyr::select(barcode, Subject, Phenotype_2, CDR3b) |>
    dplyr::left_join(by_cdr3, by = "CDR3b") |>
    dplyr::left_join(clustcr_df |> dplyr::rename(clustcr_cluster = cluster),
                     by = "barcode")
  joint
}

# Assignment-level metrics across cells where BOTH a GLIPH and a clusTCR
# label exist. NMI computed inline so we don't pull in aricode unless we
# really need it.
.gtc_nmi <- function(a, b) {
  ta <- table(a)
  tb <- table(b)
  N  <- sum(ta)
  if (N < 2) return(NA_real_)
  ent <- function(p) {
    p <- p / sum(p); p <- p[p > 0]
    -sum(p * log(p))
  }
  Ha <- ent(ta); Hb <- ent(tb)
  joint <- table(a, b)
  Hab <- ent(as.numeric(joint))
  mi  <- Ha + Hb - Hab
  if (Ha + Hb == 0) return(NA_real_)
  2 * mi / (Ha + Hb)  # normalized MI, range 0..1
}

.gtc_assignment_metrics <- function(joint) {
  keep <- !is.na(joint$gliph_primary_cluster) & !is.na(joint$clustcr_cluster)
  if (sum(keep) < 10) {
    return(data.frame(metric = c("ari", "nmi"),
                      value  = c(NA_real_, NA_real_),
                      n_cells = sum(keep)))
  }
  a <- as.character(joint$gliph_primary_cluster[keep])
  b <- as.character(joint$clustcr_cluster[keep])
  ari <- tryCatch(mclust::adjustedRandIndex(a, b), error = function(e) NA_real_)
  nmi <- tryCatch(.gtc_nmi(a, b), error = function(e) NA_real_)
  data.frame(metric  = c("ari", "nmi"),
             value   = c(ari, nmi),
             n_cells = sum(keep))
}

# Subject-stratified within vs background permutation for a single GLIPH
# cluster. Drawn pairs are restricted to subjects represented in the cluster
# to preserve the publicity structure (clones from the same subject have
# correlated TRBV usage which would inflate background distances).
.gtc_distance_perm <- function(member_bcs, pw_beta, bc_to_subj,
                               n_perm = 1000, seed = 42L) {
  if (length(member_bcs) < 2) {
    return(c(mean_within = NA_real_, mean_between = NA_real_,
             z = NA_real_, perm_p = NA_real_))
  }
  member_bcs <- intersect(member_bcs, rownames(pw_beta))
  if (length(member_bcs) < 2) {
    return(c(mean_within = NA_real_, mean_between = NA_real_,
             z = NA_real_, perm_p = NA_real_))
  }
  sub_dmat <- pw_beta[member_bcs, member_bcs, drop = FALSE]
  within_vals <- sub_dmat[upper.tri(sub_dmat)]
  mean_within <- mean(within_vals, na.rm = TRUE)

  cluster_subjects <- unique(bc_to_subj[member_bcs])
  pool <- names(bc_to_subj)[bc_to_subj %in% cluster_subjects]
  pool <- intersect(pool, rownames(pw_beta))
  if (length(pool) < length(member_bcs) * 2) {
    return(c(mean_within = mean_within, mean_between = NA_real_,
             z = NA_real_, perm_p = NA_real_))
  }

  set.seed(seed)
  k <- length(member_bcs)
  bg <- vapply(seq_len(n_perm), function(i) {
    samp <- sample(pool, size = k, replace = FALSE)
    d <- pw_beta[samp, samp, drop = FALSE]
    mean(d[upper.tri(d)], na.rm = TRUE)
  }, numeric(1))
  bg <- bg[is.finite(bg)]
  mean_between <- mean(bg)
  sd_between   <- stats::sd(bg)
  z <- if (is.finite(sd_between) && sd_between > 0)
         (mean_within - mean_between) / sd_between else NA_real_
  perm_p <- (1 + sum(bg <= mean_within)) / (1 + length(bg))
  c(mean_within = mean_within, mean_between = mean_between,
    z = z, perm_p = perm_p)
}

# Bipartite Fisher per GLIPH x clusTCR cell — large tables are sparse so we
# loop only over observed (gliph, clustcr) co-occurrences.
.gtc_overlap_fisher <- function(joint) {
  d <- joint[!is.na(joint$gliph_primary_cluster) &
             !is.na(joint$clustcr_cluster), , drop = FALSE]
  if (nrow(d) == 0) return(data.frame())
  tab <- as.data.frame(table(gliph = d$gliph_primary_cluster,
                             clustcr = d$clustcr_cluster),
                       stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  n_total <- nrow(d)
  gliph_totals   <- tapply(rep(1, nrow(d)), d$gliph_primary_cluster, sum)
  clustcr_totals <- tapply(rep(1, nrow(d)), d$clustcr_cluster, sum)
  rows <- lapply(seq_len(nrow(tab)), function(i) {
    g  <- tab$gliph[i]; cc <- tab$clustcr[i]
    a  <- tab$Freq[i]
    b  <- gliph_totals[[as.character(g)]] - a
    cv <- clustcr_totals[[as.character(cc)]] - a
    d_ <- n_total - a - b - cv
    ft <- suppressWarnings(fisher.test(matrix(c(a, cv, b, d_), nrow = 2),
                                       alternative = "greater"))
    data.frame(gliph_cluster_id = g,
               clustcr_cluster  = cc,
               n_shared         = a,
               log2OR           = log2(unname(ft$estimate) + 1e-6),
               fisher_p         = ft$p.value,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$fisher_p, method = "BH")
  out[order(out$fdr), ]
}

run_gliph_tcrdist_correlation <- function(cfg) {
  if (!isTRUE(cfg$steps$gliph_tcrdist)) {
    log_message("GLIPH x tcrdist correlation disabled. Skipping.")
    return(invisible(TRUE))
  }

  log_message("Starting GLIPH x tcrdist correlation...")

  obj_dir    <- cfg$paths$results_objects
  out_tables <- file.path(cfg$paths$results_tables, "repertoire")
  ensure_dir(out_tables)

  gliph_rds  <- file.path(obj_dir, "ImmGLIPHResults.rds")
  tcr_rds    <- file.path(obj_dir, "ImmLynxTcrdistResults.rds")
  pwb_rds    <- file.path(obj_dir, "tcrdist_pw_beta.rds")
  clustcr_csv <- file.path(out_tables, "clustcr_clusters.csv")

  if (!file.exists(gliph_rds) || !file.exists(tcr_rds)) {
    log_message("  Required RDS missing (gliph or tcrdist). Skipping.")
    return(invisible(FALSE))
  }

  gliph <- readRDS(gliph_rds)
  tcr   <- readRDS(tcr_rds)
  clustcr_df <- if (file.exists(clustcr_csv))
                  utils::read.csv(clustcr_csv, stringsAsFactors = FALSE)
                else data.frame(barcode = character(0), cluster = integer(0))

  joint <- .gtc_build_join(gliph$clusters, gliph$cluster_props, clustcr_df,
                           gliph$trb)
  saveRDS(joint, file.path(obj_dir, "GliphTcrdistJoint.rds"))
  log_message("  Joint table: ", nrow(joint), " cells, ",
              sum(!is.na(joint$gliph_primary_cluster)), " with GLIPH, ",
              sum(!is.na(joint$clustcr_cluster)), " with clusTCR.")

  # --- Assignment-level metrics ----------------------------------------------
  assign_m <- .gtc_assignment_metrics(joint)
  log_message("  ARI = ", round(assign_m$value[assign_m$metric == "ari"], 3),
              " ; NMI = ", round(assign_m$value[assign_m$metric == "nmi"], 3))

  # --- Per-cluster distance permutation --------------------------------------
  pw_beta <- if (file.exists(pwb_rds)) readRDS(pwb_rds) else NULL
  cluster_rows <- data.frame()
  if (is.null(pw_beta)) {
    log_message("  WARN: tcrdist_pw_beta.rds not found; skipping distance ",
                "permutation. Re-run R/65 with the updated saver to populate it.")
  } else {
    # Map CDR3 -> barcodes (for member lookup) and barcode -> Subject.
    cdr3_to_bcs <- split(gliph$trb$barcode, gliph$trb$CDR3b)
    bc_to_subj  <- setNames(gliph$trb$Subject, gliph$trb$barcode)

    cluster_ids <- unique(gliph$clusters$cluster_id)
    n_perm <- 1000L
    seed   <- cfg$seed %||% 42L
    rows <- lapply(cluster_ids, function(cl) {
      cdr3s <- unique(gliph$clusters$CDR3b[gliph$clusters$cluster_id == cl])
      bcs   <- unlist(cdr3_to_bcs[cdr3s], use.names = FALSE)
      if (length(bcs) < 2) return(NULL)
      bcs   <- intersect(bcs, rownames(pw_beta))
      m     <- .gtc_distance_perm(bcs, pw_beta, bc_to_subj,
                                  n_perm = n_perm, seed = seed)

      cc_overlap <- joint |>
        dplyr::filter(gliph_primary_cluster == cl,
                      !is.na(clustcr_cluster))
      dom_cc <- NA_character_
      jac    <- NA_real_
      if (nrow(cc_overlap) > 0) {
        tab <- sort(table(cc_overlap$clustcr_cluster), decreasing = TRUE)
        dom_cc <- names(tab)[1]
        a <- as.integer(tab[1])
        b <- sum(joint$clustcr_cluster == dom_cc &
                 !is.na(joint$clustcr_cluster) &
                 (is.na(joint$gliph_primary_cluster) |
                  joint$gliph_primary_cluster != cl))
        c <- nrow(cc_overlap) - a
        jac <- a / (a + b + c)
      }

      data.frame(
        cluster_id           = cl,
        n_cdr3               = length(cdr3s),
        n_cells              = length(bcs),
        mean_within_tcrdist  = unname(m["mean_within"]),
        mean_between_tcrdist = unname(m["mean_between"]),
        z                    = unname(m["z"]),
        perm_p               = unname(m["perm_p"]),
        dominant_clustcr_cluster = dom_cc,
        jaccard_clustcr      = jac,
        stringsAsFactors     = FALSE
      )
    })
    cluster_rows <- do.call(rbind, rows)
    if (!is.null(cluster_rows) && nrow(cluster_rows) > 0)
      cluster_rows$fdr <- stats::p.adjust(cluster_rows$perm_p, method = "BH")
  }

  # Append assignment-level ARI/NMI as global rows (cluster_id = "GLOBAL_ari" etc).
  global_rows <- data.frame(
    cluster_id           = paste0("GLOBAL_", assign_m$metric),
    n_cdr3               = NA_integer_,
    n_cells              = assign_m$n_cells,
    mean_within_tcrdist  = NA_real_,
    mean_between_tcrdist = NA_real_,
    z                    = assign_m$value,
    perm_p               = NA_real_,
    dominant_clustcr_cluster = NA_character_,
    jaccard_clustcr      = NA_real_,
    fdr                  = NA_real_,
    stringsAsFactors     = FALSE
  )
  if (nrow(cluster_rows) == 0) cluster_rows <- global_rows[0, , drop = FALSE]
  concordance <- rbind(cluster_rows, global_rows)
  utils::write.csv(concordance,
                   file.path(out_tables, "gliph_tcrdist_concordance.csv"),
                   row.names = FALSE)
  log_message("  Saved: gliph_tcrdist_concordance.csv (",
              nrow(concordance), " rows)")

  # --- Bipartite Fisher --------------------------------------------------------
  overlap <- .gtc_overlap_fisher(joint)
  if (nrow(overlap) == 0) {
    log_message("  No GLIPH-clusTCR co-occurrences; skipping overlap table.")
  } else {
    utils::write.csv(overlap,
                     file.path(out_tables, "gliph_clustcr_overlap.csv"),
                     row.names = FALSE)
    log_message("  Saved: gliph_clustcr_overlap.csv (", nrow(overlap), " rows)")
  }

  log_message("GLIPH x tcrdist correlation complete.")
  invisible(TRUE)
}
