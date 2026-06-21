# R/53_public_bcr_clones.R
# Cross-subject "public" IGH clone detection via fuzzy CDR3 matching.
#
# Approach: SHM-laden BCRs almost never share identical heavy-chain CDR3 across
# subjects, so exact-match public-clone detection (which works for TCR) misses
# convergent antibody responses. We bucket on V gene + J gene + CDR3 amino-acid
# length, then within each bucket compute pairwise Hamming distances on the
# junction_aa column and single-link cluster on a threshold of
# floor(max_hamming_pct * L). A cluster is "public" when it contains CDR3s from
# at least min_subjects_per_cluster distinct Subjects.
#
# Reads:  outputs/tables/bcr_airr/*_airr.tsv  (productive IGH only)
#         inputs/data/metadata.csv             (subject -> etiology, phenotype)
# Writes: outputs/tables/repertoire/BCR_public_clones.csv
#         outputs/tables/repertoire/BCR_public_cluster_summary.csv
#         outputs/tables/repertoire/BCR_public_clones_permnull.csv  (optional)
#
# Permutation null: shuffles subject labels perm_null_iters times and recomputes
# n_public_clusters. Establishes that the observed convergence rate exceeds
# what bucketing alone produces by chance.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Strip allele suffix (e.g. "IGHV3-43*02" -> "IGHV3-43"); collapse comma-listed
# multi-gene calls to the first hit, which is alakazam's convention.
.strip_allele <- function(x) {
  x <- sub(",.*$", "", as.character(x))
  sub("\\*.*$", "", x)
}

# Load productive IGH rows from all per-sample AIRR tables. Keeps only the
# columns we need for matching + downstream annotation.
.load_igh_airr <- function(cfg) {
  bcr_dir <- file.path(cfg$paths$results_tables, "bcr_airr")
  if (!dir.exists(bcr_dir)) {
    log_message("  No BCR AIRR directory at ", bcr_dir, "; skipping public-clone pass.")
    return(NULL)
  }
  files <- list.files(bcr_dir, pattern = "_airr\\.tsv$", full.names = TRUE)
  if (length(files) == 0L) return(NULL)
  log_message("  Reading ", length(files), " BCR AIRR tables for IGH fuzzy matching...")
  cols_keep <- c("cell_id_unique", "sequence_id", "subject_id", "sample_id",
                 "locus", "productive", "v_call", "j_call", "c_call",
                 "junction_aa", "mu_freq", "mu_freq_heavy", "clone_id")
  rows <- lapply(files, function(f) {
    df <- tryCatch(readr::read_tsv(f, col_types = readr::cols(.default = "c"),
                                   progress = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) return(NULL)
    have <- intersect(cols_keep, colnames(df))
    df <- df[, have, drop = FALSE]
    # Standardize productive (some pipelines write "T"/"F", others "TRUE"/"FALSE")
    if ("productive" %in% colnames(df)) {
      df$productive <- toupper(substr(df$productive, 1, 1)) %in% c("T")
    } else {
      df$productive <- TRUE
    }
    df[df$locus == "IGH" & df$productive & !is.na(df$junction_aa) &
       nchar(df$junction_aa) > 0L, , drop = FALSE]
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(NULL)
  out <- dplyr::bind_rows(rows)
  out$v_gene <- .strip_allele(out$v_call)
  out$j_gene <- .strip_allele(out$j_call)
  out$cdr3_len <- nchar(out$junction_aa)
  out
}

# Within one (v, j, length) bucket, cluster CDR3s by Hamming distance.
# Returns an integer vector of cluster ids aligned to rows. Single-linkage via
# igraph::components on the threshold graph. Diagonal not included.
.cluster_one_bucket <- function(seqs, max_hamming) {
  n <- length(seqs)
  if (n <= 1L) return(rep(1L, n))
  if (!requireNamespace("stringdist", quietly = TRUE))
    stop("Install 'stringdist' for fuzzy IGH public-clone matching.")
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Install 'igraph' for fuzzy IGH public-clone matching.")
  d <- stringdist::stringdistmatrix(seqs, seqs, method = "hamming")
  diag(d) <- NA
  hit <- which(d <= max_hamming, arr.ind = TRUE)
  if (nrow(hit) == 0L) return(seq_len(n))
  hit <- hit[hit[, 1] < hit[, 2], , drop = FALSE]   # undirected, upper triangle
  if (nrow(hit) == 0L) return(seq_len(n))
  g <- igraph::make_undirected_graph(c(t(hit)), n = n)
  comp <- igraph::components(g)$membership
  as.integer(comp)
}

# Compute public-clone clusters within a stratum (etiology or "all"). Returns
# a list of two data.frames: cell-level and cluster-level summary.
.run_public_one_stratum <- function(igh, stratum_label, cfg, perm_iters = 0L) {
  pcfg <- cfg$bcr_public
  pct  <- as.numeric(pcfg$max_hamming_pct %||% 0.10)
  min_len <- as.integer(pcfg$min_cdr3_length %||% 10L)
  min_subj <- as.integer(pcfg$min_subjects_per_cluster %||% 2L)

  d <- igh |>
    dplyr::filter(cdr3_len >= min_len,
                  !is.na(v_gene), !is.na(j_gene),
                  !is.na(subject_id), nzchar(subject_id))
  if (nrow(d) == 0L) return(NULL)

  log_message("  [", stratum_label, "] bucketing on (v_gene, j_gene, length): ",
              nrow(d), " IGH sequences across ", length(unique(d$subject_id)),
              " subjects")

  # Bucket and cluster per bucket
  d$bucket <- paste(d$v_gene, d$j_gene, d$cdr3_len, sep = "|")
  buckets <- split(seq_len(nrow(d)), d$bucket)
  cluster_assign <- integer(nrow(d))
  next_id <- 1L
  for (b in names(buckets)) {
    idx <- buckets[[b]]
    seqs <- d$junction_aa[idx]
    L <- d$cdr3_len[idx[1]]
    max_h <- floor(pct * L)
    if (max_h < 1L) max_h <- 1L
    cl <- .cluster_one_bucket(seqs, max_hamming = max_h)
    cluster_assign[idx] <- cl + (next_id - 1L)
    next_id <- next_id + max(cl)
  }
  d$cluster_id <- cluster_assign

  # Per cluster: how many distinct subjects, etc.
  per_cluster <- d |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarise(
      n_cells          = dplyr::n(),
      n_subjects       = dplyr::n_distinct(subject_id),
      v_gene           = dplyr::first(v_gene),
      j_gene           = dplyr::first(j_gene),
      cdr3_len         = dplyr::first(cdr3_len),
      cdr3_consensus   = names(sort(table(junction_aa), decreasing = TRUE))[1],
      isotypes         = paste(sort(unique(stats::na.omit(.strip_allele(c_call)))),
                               collapse = ","),
      subjects         = paste(sort(unique(subject_id)), collapse = ","),
      .groups = "drop"
    )

  pub_clusters <- per_cluster |>
    dplyr::filter(n_subjects >= min_subj)
  log_message("  [", stratum_label, "] public clusters (>= ", min_subj,
              " subjects): ", nrow(pub_clusters), " of ", nrow(per_cluster),
              " total clusters")

  cell_rows <- d |>
    dplyr::semi_join(pub_clusters, by = "cluster_id") |>
    dplyr::select(dplyr::any_of(c("cell_id_unique","sequence_id","subject_id",
                                  "sample_id","v_gene","j_gene","cdr3_len",
                                  "junction_aa","c_call","mu_freq",
                                  "mu_freq_heavy","cluster_id"))) |>
    dplyr::mutate(stratum = stratum_label)
  pub_clusters$stratum <- stratum_label

  # Permutation null
  null_out <- NULL
  if (perm_iters > 0L && nrow(d) > 0L) {
    set.seed(cfg$seed %||% 42L)
    null_counts <- integer(perm_iters)
    for (k in seq_len(perm_iters)) {
      shuf_subjects <- sample(d$subject_id)   # permute subject labels
      tmp <- data.frame(cluster_id = d$cluster_id, subject_id = shuf_subjects)
      tmp_pc <- tmp |>
        dplyr::group_by(cluster_id) |>
        dplyr::summarise(n_subjects = dplyr::n_distinct(subject_id),
                         .groups = "drop")
      null_counts[k] <- sum(tmp_pc$n_subjects >= min_subj)
    }
    null_out <- data.frame(
      stratum     = stratum_label,
      observed    = nrow(pub_clusters),
      null_mean   = mean(null_counts),
      null_sd     = stats::sd(null_counts),
      null_q025   = stats::quantile(null_counts, 0.025, names = FALSE),
      null_q975   = stats::quantile(null_counts, 0.975, names = FALSE),
      p_one_sided = (sum(null_counts >= nrow(pub_clusters)) + 1) /
                     (perm_iters + 1)
    )
    log_message("  [", stratum_label, "] permutation null: obs=",
                nrow(pub_clusters), " null_mean=",
                signif(null_out$null_mean, 3), " p=",
                signif(null_out$p_one_sided, 3))

    # Complementary stat: mean within-cluster CDR3-aa Hamming identity
    # distance of PUBLIC clusters under observed labels. Low value = high
    # convergence (sequences within a public cluster are nearly identical).
    # Random V/J/L bucketing alone wouldn't produce low within-cluster
    # distances, so this stat shows whether the observed public clusters
    # represent true convergent CDR3s vs noise.
    mean_identity_dist <- function(cluster_assign, seqs) {
      pubs <- names(which(table(cluster_assign) >= 2))
      if (length(pubs) == 0L) return(NA_real_)
      vapply(pubs, function(cl) {
        s <- seqs[cluster_assign == cl]
        if (length(s) < 2L) return(NA_real_)
        cm <- utils::combn(s, 2L)
        mean(stringdist::stringdist(cm[1, ], cm[2, ],
                                    method = "hamming") /
             nchar(cm[1, ]))
      }, numeric(1)) |> mean(na.rm = TRUE)
    }
    pub_mask <- d$cluster_id %in% pub_clusters$cluster_id
    obs_identity <- if (any(pub_mask)) {
      mean_identity_dist(d$cluster_id[pub_mask], d$junction_aa[pub_mask])
    } else NA_real_
    null_out$obs_mean_identity_distance <- round(obs_identity, 4)
    null_out$null_test_direction <- "observed_below_or_above"
    null_out$interpretation <- paste(
      "Observed n_public_clusters reflects subject-label-aware convergence.",
      "If observed < null_mean, V/J/L bucketing alone explains apparent sharing.",
      "If observed > null_mean, true convergent CDR3s drive the signal.",
      "obs_mean_identity_distance is the within-cluster CDR3 Hamming",
      "fraction of public clusters; near 0 = identical CDR3s; near max_hamming_pct",
      "= drift-limited convergence.")
  }

  list(cells = cell_rows, clusters = pub_clusters, null = null_out)
}

# Attach Subject -> Etiology / Phenotype_2 from the metadata table.
.attach_subject_metadata <- function(igh, cfg) {
  meta_path <- "inputs/data/metadata.csv"
  if (!file.exists(meta_path)) {
    log_message("  metadata.csv not found at ", meta_path,
                "; etiology stratification disabled.")
    igh$Etiology <- NA_character_
    igh$Phenotype_2 <- NA_character_
    return(igh)
  }
  m <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
  # subject_id in AIRR is the Subject in metadata (verified against ARN* samples)
  subj_col <- intersect(c("Subject", "subject_id", "SubjectID"), colnames(m))[1]
  if (is.na(subj_col)) {
    log_message("  metadata.csv has no Subject column; etiology stratification disabled.")
    igh$Etiology <- NA_character_
    igh$Phenotype_2 <- NA_character_
    return(igh)
  }
  keep_cols <- c(subj_col,
                 intersect(c("Etiology", "Phenotype_2"), colnames(m)))
  m <- m[, keep_cols, drop = FALSE] |> dplyr::distinct()
  colnames(m)[1] <- "subject_id"
  igh <- dplyr::left_join(igh, m, by = "subject_id")
  igh
}

run_public_bcr_clones <- function(cfg) {
  if (!isTRUE(cfg$bcr_public$enable)) {
    log_message("bcr_public disabled in config; skipping.")
    return(invisible(FALSE))
  }
  log_message("=== Public BCR (IGH) clones across subjects ===")
  igh <- .load_igh_airr(cfg)
  if (is.null(igh) || nrow(igh) == 0L) {
    log_message("  No IGH sequences available; aborting.")
    return(invisible(FALSE))
  }
  igh <- .attach_subject_metadata(igh, cfg)

  # Optional subject_subset for smoke tests
  subset_subj <- cfg$bcr_public$subject_subset %||% NULL
  if (!is.null(subset_subj) && length(subset_subj) > 0L) {
    igh <- igh |> dplyr::filter(subject_id %in% subset_subj)
    log_message("  subject_subset active: ", length(unique(igh$subject_id)),
                " subjects kept.")
  }

  out_dir <- file.path("outputs/tables/repertoire")
  ensure_dir(out_dir)
  perm_iters <- as.integer(cfg$bcr_public$perm_null_iters %||% 0L)

  results <- list()
  null_rows <- list()

  if (isTRUE(cfg$bcr_public$require_same_etiology) &&
      "Phenotype_2" %in% colnames(igh) && any(!is.na(igh$Phenotype_2))) {
    for (cond in c("NIU", "Viral")) {
      sub <- igh |> dplyr::filter(Phenotype_2 == cond)
      if (nrow(sub) == 0L) next
      r <- .run_public_one_stratum(sub, cond, cfg, perm_iters = perm_iters)
      if (!is.null(r)) {
        results[[cond]] <- r
        if (!is.null(r$null)) null_rows[[cond]] <- r$null
      }
    }
  }
  if (isTRUE(cfg$bcr_public$also_run_cross_etiology) ||
      length(results) == 0L) {
    r <- .run_public_one_stratum(igh, "all", cfg, perm_iters = perm_iters)
    if (!is.null(r)) {
      results[["all"]] <- r
      if (!is.null(r$null)) null_rows[["all"]] <- r$null
    }
  }

  if (length(results) == 0L) {
    log_message("  No public clusters found in any stratum.")
    return(invisible(FALSE))
  }

  cells_all <- dplyr::bind_rows(lapply(results, function(x) x$cells))
  clust_all <- dplyr::bind_rows(lapply(results, function(x) x$clusters))

  cells_path <- file.path(out_dir, "BCR_public_clones.csv")
  clust_path <- file.path(out_dir, "BCR_public_cluster_summary.csv")
  utils::write.csv(cells_all, cells_path, row.names = FALSE)
  utils::write.csv(clust_all, clust_path, row.names = FALSE)
  log_message("  Wrote: ", cells_path, " (", nrow(cells_all), " rows)")
  log_message("  Wrote: ", clust_path, " (", nrow(clust_all), " clusters)")

  if (length(null_rows) > 0L) {
    null_df <- dplyr::bind_rows(null_rows)
    null_path <- file.path(out_dir, "BCR_public_clones_permnull.csv")
    utils::write.csv(null_df, null_path, row.names = FALSE)
    log_message("  Wrote: ", null_path)
  }
  invisible(TRUE)
}

# run_public_bcr_clones reads productive IGH rows from every
# outputs/tables/bcr_airr/*_airr.tsv, attaches Subject + Etiology + Phenotype_2
# from inputs/data/metadata.csv, then within each stratum (NIU, Viral, and
# optionally "all" cross-etiology) buckets sequences on V-gene + J-gene +
# CDR3-aa length and single-link clusters within each bucket using Hamming
# distance with threshold = floor(max_hamming_pct * L). A cluster is "public"
# when it contains CDR3s from min_subjects_per_cluster or more distinct
# Subjects. Outputs a cell-level table (one row per IGH sequence in a public
# cluster) and a cluster-level summary, plus a permutation null where Subject
# labels are shuffled to confirm that observed convergence exceeds the rate
# expected by bucketing alone.
