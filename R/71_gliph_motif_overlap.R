# R/71_gliph_motif_overlap.R
# Greedy overlap-layout-consensus (OLC) assembly of GLIPH motifs into
# longer "meta-motifs" — like the published VGKWY example (LNWQ, NWQG,
# KWYG, GKWY, VGKW -> assembled VGKWY).
#
# Pipeline:
#   1. Pull motif strings from gliph_cluster_properties.csv ($motif).
#   2. Strip GLIPH wildcards ("%", "{aa}") and uppercase. Drop motifs < min_motif_len.
#   3. Compute pairwise suffix(a, k) == prefix(b, k) with k maximal,
#      k >= min_overlap.
#   4. Weight each edge by shared expanded + 0.5 * shared public CDR3s
#      across the two motifs' constituent sequences.
#   5. Greedy assembly: start from highest-weight tail-of-source vertex,
#      repeatedly extend to highest-weight outgoing edge; stop on cycle,
#      length cap, or weight floor.
#   6. Consensus = position-frequency over the constituent motifs at their
#      assembly offsets.
#
# Outputs:
#   outputs/tables/repertoire/motif_overlap_edges.csv
#   outputs/tables/repertoire/meta_motif_assemblies.csv

`%||%` <- function(x, y) if (is.null(x)) y else x

.gmo_clean_motif <- function(s, min_len = 4L) {
  s <- as.character(s)
  s <- gsub("\\{[^}]+\\}", "", s)   # strip "{aa}" wildcards
  s <- gsub("%", "",        s)       # strip "%" wildcards
  s <- gsub("[^A-Z]", "",   toupper(s))
  s[nchar(s) >= min_len]
}

# For ordered (a,b), find largest k such that suffix(a, k) == prefix(b, k).
.gmo_max_overlap <- function(a, b, min_overlap = 3L) {
  n <- min(nchar(a), nchar(b))
  if (n < min_overlap) return(0L)
  for (k in n:min_overlap) {
    if (substr(a, nchar(a) - k + 1L, nchar(a)) == substr(b, 1L, k))
      return(k)
  }
  0L
}

.gmo_pairwise_overlap <- function(motifs, min_overlap = 3L) {
  if (length(motifs) < 2) return(data.frame())
  rows <- list()
  for (i in seq_along(motifs)) {
    for (j in seq_along(motifs)) {
      if (i == j) next
      k <- .gmo_max_overlap(motifs[i], motifs[j], min_overlap)
      if (k > 0L)
        rows[[length(rows) + 1]] <- data.frame(
          motif_a     = motifs[i],
          motif_b     = motifs[j],
          overlap_len = k,
          stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# Build CDR3 -> {expanded, public} lookups, then attribute clone counts to
# each motif by membership (gliph_clusters$motif -> CDR3b list).
.gmo_motif_to_clones <- function(gliph_clusters, gliph_props,
                                 trb_frame) {
  # Motif maps to all CDR3b across all clusters tagged with that motif.
  m2c <- gliph_clusters |>
    dplyr::left_join(
      dplyr::select(gliph_props, cluster_id, motif_raw = motif),
      by = "cluster_id") |>
    dplyr::select(motif_raw, CDR3b) |>
    dplyr::distinct()
  m2c
}

# For a single (motif_a, motif_b) edge, count shared CDR3s, shared expanded,
# shared public — relative to the per-CDR3 boolean lookup tables.
.gmo_weight_edges <- function(edges, motif_to_clones,
                              cdr3_expanded, cdr3_public,
                              motif_cleaner) {
  if (nrow(edges) == 0) return(edges)
  # motif_cleaner maps cleaned -> raw motif string used in gliph_clusters
  by_clean <- split(motif_to_clones, motif_cleaner[motif_to_clones$motif_raw])
  cdr3_by_clean <- lapply(by_clean, function(d) unique(d$CDR3b))

  shared_n  <- integer(nrow(edges))
  shared_ex <- integer(nrow(edges))
  shared_pu <- integer(nrow(edges))
  for (i in seq_len(nrow(edges))) {
    a_cdr3 <- cdr3_by_clean[[edges$motif_a[i]]] %||% character(0)
    b_cdr3 <- cdr3_by_clean[[edges$motif_b[i]]] %||% character(0)
    sh <- intersect(a_cdr3, b_cdr3)
    shared_n[i]  <- length(sh)
    shared_ex[i] <- sum(sh %in% cdr3_expanded)
    shared_pu[i] <- sum(sh %in% cdr3_public)
  }
  edges$shared_cdr3     <- shared_n
  edges$shared_expanded <- shared_ex
  edges$shared_public   <- shared_pu
  edges$edge_weight     <- shared_ex + 0.5 * shared_pu
  edges
}

# Greedy assembly: chase highest-weight outgoing edges from each starting
# motif; mark used edges; stop on cycle, cap length, or hit min_edge_weight floor.
.gmo_assemble <- function(edges, motifs,
                          max_len = 20L, min_edge_weight = 1) {
  if (nrow(edges) == 0) return(list())
  e <- edges[edges$edge_weight >= min_edge_weight, , drop = FALSE]
  if (nrow(e) == 0) return(list())
  e <- e[order(-e$edge_weight), , drop = FALSE]

  used <- rep(FALSE, nrow(e))
  visited_motifs <- character(0)
  assemblies <- list()
  next_id <- 1L

  while (any(!used) && length(visited_motifs) < length(motifs)) {
    seed_idx <- which(!used)[1]
    if (is.na(seed_idx)) break
    contig   <- e$motif_a[seed_idx]
    offsets  <- 0L
    consts   <- c(contig)
    used[seed_idx] <- TRUE

    # Extend forward
    cur <- e$motif_b[seed_idx]
    used_offset <- nchar(e$motif_a[seed_idx]) - e$overlap_len[seed_idx]
    contig <- paste0(contig, substr(e$motif_b[seed_idx],
                                    e$overlap_len[seed_idx] + 1L,
                                    nchar(e$motif_b[seed_idx])))
    consts  <- c(consts, cur)
    offsets <- c(offsets, used_offset)

    while (nchar(contig) < max_len) {
      cand <- which(!used & e$motif_a == cur)
      if (length(cand) == 0) break
      best <- cand[which.max(e$edge_weight[cand])[1]]
      if (e$edge_weight[best] < min_edge_weight) break
      used[best] <- TRUE
      ext <- substr(e$motif_b[best], e$overlap_len[best] + 1L,
                    nchar(e$motif_b[best]))
      if (nchar(ext) == 0) break
      used_offset <- nchar(contig) - e$overlap_len[best]
      contig <- paste0(contig, ext)
      cur <- e$motif_b[best]
      consts  <- c(consts, cur)
      offsets <- c(offsets, used_offset)
      if (cur %in% consts[-length(consts)]) break  # cycle guard
    }

    visited_motifs <- union(visited_motifs, consts)
    assemblies[[next_id]] <- list(
      meta_motif_id     = paste0("MM", sprintf("%03d", next_id)),
      consensus         = contig,
      constituent_motifs = consts,
      position_offsets  = offsets)
    next_id <- next_id + 1L
  }
  assemblies
}

# Position-frequency consensus across stacked aligned motifs.
.gmo_consensus <- function(constituents, offsets) {
  max_end <- max(offsets + nchar(constituents))
  mat <- matrix("", nrow = length(constituents), ncol = max_end)
  for (i in seq_along(constituents)) {
    cc <- strsplit(constituents[i], "")[[1]]
    cols <- seq_along(cc) + offsets[i]
    mat[i, cols] <- cc
  }
  cons <- vapply(seq_len(ncol(mat)), function(j) {
    col <- mat[, j]; col <- col[col != ""]
    if (length(col) == 0) return("")
    tab <- sort(table(col), decreasing = TRUE)
    names(tab)[1]
  }, character(1))
  paste0(cons, collapse = "")
}

# Summarize per-meta-motif: subjects, clones, dominant phenotype.
.gmo_summarize <- function(assemblies, motif_to_clones, motif_cleaner,
                           trb_frame) {
  if (length(assemblies) == 0) return(data.frame())
  rows <- lapply(assemblies, function(a) {
    cdr3s <- unique(motif_to_clones$CDR3b[
      motif_cleaner[motif_to_clones$motif_raw] %in% a$constituent_motifs])
    sub_rows <- trb_frame[trb_frame$CDR3b %in% cdr3s, , drop = FALSE]
    data.frame(
      meta_motif_id      = a$meta_motif_id,
      consensus          = a$consensus,
      constituent_motifs = paste(a$constituent_motifs, collapse = ";"),
      position_offsets   = paste(a$position_offsets, collapse = ";"),
      n_constituents     = length(a$constituent_motifs),
      n_clones           = length(cdr3s),
      n_subjects         = if (nrow(sub_rows) > 0)
                             dplyr::n_distinct(sub_rows$Subject) else 0L,
      dominant_phenotype = if (nrow(sub_rows) > 0)
        names(sort(table(sub_rows$Phenotype_2), decreasing = TRUE))[1]
        else NA_character_,
      stringsAsFactors   = FALSE)
  })
  do.call(rbind, rows)
}

run_gliph_motif_overlap <- function(cfg) {
  if (!isTRUE(cfg$steps$gliph_motif_overlap)) {
    log_message("GLIPH motif-overlap assembly disabled. Skipping.")
    return(invisible(TRUE))
  }
  log_message("Starting GLIPH motif-overlap assembly...")

  out_tables <- file.path(cfg$paths$results_tables, "repertoire")
  ensure_dir(out_tables)

  gliph_rds <- file.path(cfg$paths$results_objects, "ImmGLIPHResults.rds")
  if (!file.exists(gliph_rds)) {
    log_message("  ImmGLIPHResults.rds missing. Skipping.")
    return(invisible(FALSE))
  }
  gliph <- readRDS(gliph_rds)

  min_len <- cfg$tcr_advanced$motif$min_motif_len %||% 4L
  min_ov  <- cfg$tcr_advanced$motif$min_overlap   %||% 3L
  max_len <- cfg$tcr_advanced$motif$max_assembly_len %||% 20L
  min_w   <- cfg$tcr_advanced$motif$min_edge_weight  %||% 1

  raw_motifs <- unique(gliph$cluster_props$motif)
  cleaned    <- .gmo_clean_motif(raw_motifs, min_len = min_len)
  # Map raw -> cleaned. Keep only motifs that survived the cleaner.
  motif_cleaner <- setNames(.gmo_clean_motif(raw_motifs, min_len = 1L),
                            raw_motifs)
  motif_cleaner <- motif_cleaner[nchar(motif_cleaner) >= min_len]
  uniq_clean <- unique(unname(motif_cleaner))
  log_message("  Motifs: ", length(raw_motifs), " raw, ",
              length(uniq_clean), " unique after cleaning (min_len=", min_len, ")")
  if (length(uniq_clean) < 2) {
    log_message("  Too few motifs to assemble."); return(invisible(FALSE))
  }

  edges <- .gmo_pairwise_overlap(uniq_clean, min_overlap = min_ov)
  log_message("  Pairwise overlap edges (k>=", min_ov, "): ", nrow(edges))
  if (nrow(edges) == 0) return(invisible(FALSE))

  # Build expanded / public CDR3 lookups from the eye T cell object.
  paths_tcell <- get_target_paths(cfg, "tcell")
  tcell_rds   <- file.path(paths_tcell$results_objects,
                           "IntegratedSeuratObject.rds")
  cdr3_expanded <- character(0)
  cdr3_public   <- character(0)
  if (file.exists(tcell_rds)) {
    obj_t <- readRDS(tcell_rds)
    meta <- obj_t@meta.data
    trb  <- stringr::str_split(meta$CTaa, "_", simplify = TRUE)
    trb_aa <- if (ncol(trb) >= 2) trb[, 2] else character(nrow(meta))
    df_t <- data.frame(
      CDR3b = trb_aa,
      Subject = as.character(meta$Subject),
      CTstrict = as.character(meta$CTstrict),
      cloneSize = as.character(meta$cloneSize),
      stringsAsFactors = FALSE)
    df_t <- df_t[!is.na(df_t$CDR3b) & df_t$CDR3b != "", , drop = FALSE]
    cdr3_expanded <- unique(df_t$CDR3b[df_t$cloneSize %in%
                                       c("Large", "Hyperexpanded")])
    by_ct <- tapply(df_t$Subject, df_t$CTstrict,
                    function(x) length(unique(x[!is.na(x)])))
    public_cts <- names(by_ct)[by_ct >= 2L]
    cdr3_public <- unique(df_t$CDR3b[df_t$CTstrict %in% public_cts])
    rm(obj_t); invisible(gc())
  }
  log_message("  Expanded CDR3b set: ", length(cdr3_expanded),
              "; public CDR3b set: ",   length(cdr3_public))

  m2c <- .gmo_motif_to_clones(gliph$clusters, gliph$cluster_props, gliph$trb)
  edges <- .gmo_weight_edges(edges, m2c, cdr3_expanded, cdr3_public,
                             motif_cleaner)

  utils::write.csv(edges,
                   file.path(out_tables, "motif_overlap_edges.csv"),
                   row.names = FALSE)
  log_message("  Saved: motif_overlap_edges.csv (", nrow(edges), " rows)")

  asm <- .gmo_assemble(edges, uniq_clean,
                       max_len = max_len, min_edge_weight = min_w)
  if (length(asm) == 0) {
    log_message("  No assemblies cleared min_edge_weight=", min_w)
    return(invisible(FALSE))
  }
  log_message("  Assembled ", length(asm), " meta-motifs")

  summary_df <- .gmo_summarize(asm, m2c, motif_cleaner, gliph$trb)
  utils::write.csv(summary_df,
                   file.path(out_tables, "meta_motif_assemblies.csv"),
                   row.names = FALSE)
  log_message("  Saved: meta_motif_assemblies.csv (", nrow(summary_df), " rows)")

  saveRDS(list(edges = edges,
               assemblies = asm,
               summary = summary_df,
               motif_cleaner = motif_cleaner),
          file.path(cfg$paths$results_objects, "GliphMotifAssemblies.rds"))
  log_message("  Saved: GliphMotifAssemblies.rds")

  log_message("GLIPH motif-overlap assembly complete.")
  invisible(TRUE)
}
