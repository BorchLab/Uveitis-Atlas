# R/68_gliph_network_metrics.R
# Topology of each GLIPH convergence group, computed twice:
#   (a) membership-binary  — complete clique among CDR3s in the cluster
#   (b) tcrdist-weighted   — edges among first-order neighbors within
#                            pw_beta < cutoff; weight = 1/(1 + pw_beta)
#
# Per-cluster graph metrics: density, transitivity, mean shortest path,
# modularity + n_communities (Louvain), assortativity by Phenotype_2.
# Per-node metrics: degree, betweenness, closeness, eigen_centrality, local
# clustering coefficient. Hubs (top-k per cluster) carry cloneSize, n_subjects,
# and log10_Pgen annotations for downstream interpretation.
#
# Inputs:
#   outputs/objects/ImmGLIPHResults.rds        — gliph clusters + trb frame
#   outputs/objects/ImmLynxTcrdistResults.rds  — sce
#   outputs/objects/tcrdist_pw_beta.rds        — barcode-indexed distance matrix
#   outputs/tables/repertoire/olga_pgen_per_clone.csv
#   outputs/objects/IntegratedSeuratObject.rds — for cloneSize + n_subjects
# Outputs:
#   outputs/tables/repertoire/gliph_network_metrics.csv
#   outputs/tables/repertoire/gliph_network_node_attrs.csv
#   outputs/tables/repertoire/gliph_network_hubs.csv

`%||%` <- function(x, y) if (is.null(x)) y else x

# Build the within-cluster graph for one cluster.
#   edge_type "binary"   -> clique among unique CDR3b in the cluster
#   edge_type "tcrdist"  -> CDR3-level graph weighted by mean pw_beta across
#                           cells carrying that CDR3 (one node per CDR3),
#                           edges retained when distance < cutoff
.gnm_build_graph <- function(cluster_id, members, pw_beta, bc_by_cdr3,
                             edge_type = c("binary", "tcrdist"),
                             tcrdist_cutoff = 12.5) {
  edge_type <- match.arg(edge_type)
  cdr3s <- unique(members$CDR3b)
  if (length(cdr3s) < 2) return(NULL)

  if (edge_type == "binary") {
    el <- t(utils::combn(cdr3s, 2))
    g <- igraph::graph_from_edgelist(el, directed = FALSE)
    igraph::E(g)$weight <- 1
    return(g)
  }

  # tcrdist: CDR3-pair distance = min(pw_beta) across cell pairs (or median).
  # Min is a faithful "closest representative" — protects against subjects
  # who happen to share a CDR3 from inflating distance via biological noise.
  cdr3s <- cdr3s[cdr3s %in% names(bc_by_cdr3)]
  if (length(cdr3s) < 2) return(NULL)
  bcs <- lapply(cdr3s, function(c) intersect(bc_by_cdr3[[c]], rownames(pw_beta)))
  names(bcs) <- cdr3s
  cdr3s <- cdr3s[lengths(bcs) > 0]
  if (length(cdr3s) < 2) return(NULL)

  pairs <- utils::combn(cdr3s, 2)
  edges <- vector("list", ncol(pairs))
  for (i in seq_len(ncol(pairs))) {
    a <- pairs[1, i]; b <- pairs[2, i]
    sub <- pw_beta[bcs[[a]], bcs[[b]], drop = FALSE]
    md  <- min(sub, na.rm = TRUE)
    if (is.finite(md) && md < tcrdist_cutoff) {
      edges[[i]] <- data.frame(from = a, to = b,
                               distance = md,
                               weight   = 1 / (1 + md),
                               stringsAsFactors = FALSE)
    }
  }
  edges <- do.call(rbind, edges)
  if (is.null(edges) || nrow(edges) == 0) return(NULL)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                     vertices = data.frame(name = cdr3s))
  g
}

.gnm_node_metrics <- function(g) {
  if (is.null(g) || igraph::vcount(g) == 0) return(data.frame())
  data.frame(
    CDR3b              = igraph::V(g)$name,
    degree             = igraph::degree(g),
    betweenness        = igraph::betweenness(g, normalized = TRUE),
    closeness          = igraph::closeness(g, normalized = TRUE),
    eigen_centrality   = igraph::eigen_centrality(g)$vector,
    local_clustering   = igraph::transitivity(g, type = "local", isolates = "zero"),
    stringsAsFactors   = FALSE
  )
}

.gnm_graph_metrics <- function(g, members, edge_type) {
  if (is.null(g) || igraph::vcount(g) == 0) {
    return(data.frame(edge_type = edge_type, n_nodes = 0, n_edges = 0,
                      density = NA, transitivity = NA, mean_distance = NA,
                      modularity = NA, n_communities = NA,
                      assortativity_phenotype = NA))
  }
  louvain <- tryCatch(igraph::cluster_louvain(g, weights = igraph::E(g)$weight),
                      error = function(e) NULL)
  mod     <- if (!is.null(louvain)) igraph::modularity(louvain) else NA_real_
  n_comm  <- if (!is.null(louvain)) length(louvain) else NA_integer_

  # Per-node phenotype label: modal Phenotype_2 of the cells carrying that CDR3.
  phen_lookup <- members |>
    dplyr::group_by(CDR3b) |>
    dplyr::summarise(
      phen = names(sort(table(Phenotype_2), decreasing = TRUE))[1],
      .groups = "drop")
  v_names <- igraph::V(g)$name
  phen_v  <- phen_lookup$phen[match(v_names, phen_lookup$CDR3b)]
  phen_v[is.na(phen_v)] <- "Unknown"
  assort <- tryCatch(
    igraph::assortativity_nominal(g, as.integer(factor(phen_v))),
    error = function(e) NA_real_)

  data.frame(
    edge_type              = edge_type,
    n_nodes                = igraph::vcount(g),
    n_edges                = igraph::ecount(g),
    density                = igraph::edge_density(g),
    transitivity           = igraph::transitivity(g, type = "global"),
    mean_distance          = tryCatch(igraph::mean_distance(g),
                                       error = function(e) NA_real_),
    modularity             = mod,
    n_communities          = n_comm,
    assortativity_phenotype = assort
  )
}

# Annotate each node row with clone-level metadata (sharing, expansion, Pgen).
.gnm_hub_annotate <- function(node_df, trb_frame, per_clone_pgen,
                              seurat_meta = NULL) {
  if (nrow(node_df) == 0) return(node_df)
  # n_subjects from trb_frame (one row per cell)
  shar <- trb_frame |>
    dplyr::group_by(CDR3b) |>
    dplyr::summarise(
      n_subjects = dplyr::n_distinct(Subject),
      n_cells    = dplyr::n(),
      dominant_phenotype =
        names(sort(table(Phenotype_2), decreasing = TRUE))[1],
      .groups = "drop")
  node_df <- dplyr::left_join(node_df, shar, by = "CDR3b")

  pgen <- per_clone_pgen |>
    dplyr::group_by(CDR3b) |>
    dplyr::summarise(log10_Pgen = stats::median(log10_Pgen, na.rm = TRUE),
                     .groups = "drop")
  node_df <- dplyr::left_join(node_df, pgen, by = "CDR3b")

  if (!is.null(seurat_meta) && "CTaa" %in% colnames(seurat_meta) &&
      "cloneSize" %in% colnames(seurat_meta)) {
    trb_aa <- stringr::str_split(seurat_meta$CTaa, "_", simplify = TRUE)
    if (ncol(trb_aa) >= 2) {
      cs <- data.frame(CDR3b = trb_aa[, 2],
                       cloneSize = as.character(seurat_meta$cloneSize),
                       stringsAsFactors = FALSE)
      cs <- cs[!is.na(cs$CDR3b) & cs$CDR3b != "", , drop = FALSE]
      cs <- cs |>
        dplyr::group_by(CDR3b) |>
        dplyr::summarise(
          cloneSize = names(sort(table(cloneSize), decreasing = TRUE))[1],
          .groups = "drop")
      node_df <- dplyr::left_join(node_df, cs, by = "CDR3b")
    }
  }
  node_df
}

run_gliph_network_metrics <- function(cfg) {
  if (!isTRUE(cfg$steps$gliph_network)) {
    log_message("GLIPH network metrics disabled. Skipping.")
    return(invisible(TRUE))
  }
  log_message("Starting GLIPH network metrics...")

  obj_dir    <- cfg$paths$results_objects
  out_tables <- file.path(cfg$paths$results_tables, "repertoire")
  ensure_dir(out_tables)

  gliph_rds <- file.path(obj_dir, "ImmGLIPHResults.rds")
  pwb_rds   <- file.path(obj_dir, "tcrdist_pw_beta.rds")
  pgen_csv  <- file.path(out_tables, "olga_pgen_per_clone.csv")
  seu_rds   <- file.path(obj_dir, "IntegratedSeuratObject.rds")

  if (!file.exists(gliph_rds)) {
    log_message("  ImmGLIPHResults.rds not found. Skipping."); return(invisible(FALSE))
  }
  gliph <- readRDS(gliph_rds)

  cutoff    <- cfg$tcr_advanced$network$tcrdist_cutoff %||% 12.5
  top_hubs  <- cfg$tcr_advanced$network$top_hubs_k     %||% 20
  min_size  <- cfg$tcr_advanced$network$min_cluster_size %||% 3

  pw_beta <- if (file.exists(pwb_rds)) readRDS(pwb_rds) else NULL
  if (is.null(pw_beta))
    log_message("  WARN: tcrdist_pw_beta.rds missing — tcrdist-weighted graphs will be skipped.")
  per_clone <- if (file.exists(pgen_csv))
                 utils::read.csv(pgen_csv, stringsAsFactors = FALSE)
               else data.frame(CDR3b = character(0), log10_Pgen = numeric(0))

  seurat_meta <- NULL
  if (file.exists(seu_rds)) {
    obj <- readRDS(seu_rds)
    seurat_meta <- obj@meta.data
    rm(obj); invisible(gc())
  }

  bc_by_cdr3 <- split(gliph$trb$barcode, gliph$trb$CDR3b)

  cluster_ids <- unique(gliph$clusters$cluster_id)
  graph_rows <- list()
  node_rows  <- list()

  for (cl in cluster_ids) {
    members <- gliph$clusters[gliph$clusters$cluster_id == cl, , drop = FALSE]
    members <- dplyr::left_join(
      members,
      dplyr::select(gliph$trb, CDR3b, Phenotype_2, Subject),
      by = "CDR3b", relationship = "many-to-many")
    if (length(unique(members$CDR3b)) < min_size) next

    for (et in c("binary", "tcrdist")) {
      if (et == "tcrdist" && is.null(pw_beta)) next
      g <- .gnm_build_graph(cl, members, pw_beta, bc_by_cdr3,
                            edge_type = et, tcrdist_cutoff = cutoff)
      if (is.null(g)) next

      gm <- .gnm_graph_metrics(g, members, et)
      gm$cluster_id <- cl
      graph_rows[[paste(cl, et, sep = "::")]] <- gm

      nm <- .gnm_node_metrics(g)
      if (nrow(nm) > 0) {
        nm$cluster_id <- cl
        nm$edge_type  <- et
        node_rows[[paste(cl, et, sep = "::")]] <- nm
      }
    }
  }

  if (length(graph_rows) == 0) {
    log_message("  No graphs constructed (cluster sizes too small).")
    return(invisible(FALSE))
  }

  graphs_df <- dplyr::bind_rows(graph_rows) |>
    dplyr::select(cluster_id, edge_type, dplyr::everything())
  utils::write.csv(graphs_df,
                   file.path(out_tables, "gliph_network_metrics.csv"),
                   row.names = FALSE)
  log_message("  Saved: gliph_network_metrics.csv (", nrow(graphs_df), " rows)")

  nodes_df <- dplyr::bind_rows(node_rows)
  nodes_df <- .gnm_hub_annotate(nodes_df, gliph$trb, per_clone, seurat_meta)
  utils::write.csv(nodes_df,
                   file.path(out_tables, "gliph_network_node_attrs.csv"),
                   row.names = FALSE)
  log_message("  Saved: gliph_network_node_attrs.csv (", nrow(nodes_df), " rows)")

  hubs <- nodes_df |>
    dplyr::group_by(cluster_id, edge_type) |>
    dplyr::arrange(dplyr::desc(eigen_centrality), .by_group = TRUE) |>
    dplyr::slice_head(n = top_hubs) |>
    dplyr::ungroup()
  utils::write.csv(hubs,
                   file.path(out_tables, "gliph_network_hubs.csv"),
                   row.names = FALSE)
  log_message("  Saved: gliph_network_hubs.csv (", nrow(hubs), " rows, top-",
              top_hubs, " per cluster x edge_type)")

  log_message("GLIPH network metrics complete.")
  invisible(TRUE)
}
