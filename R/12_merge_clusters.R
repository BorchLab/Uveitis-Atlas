# R/12_merge_clusters.R
# ------------------------------------------------------------------------------
# Lightweight relabeling: with 12-20 well-annotated Leiden clusters (plus any
# SNN refinement splits from 11_annotate_celltypes.R), we apply cluster-level
# celltype labels directly as merged.celltype.cluster. No overcluster merging
# needed. Includes NK rescue as a safety net.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- NK rescue ----------------------------------------------------------------
# Safety-net cluster-level promotion: if >= majority_frac of a cluster's cells
# were marked NK by 35's per-cell recovery but the mapping row still says
# something else, promote the cluster to NK. This catches clusters that 35's
# own majority-promotion step missed (e.g. edge cases around the threshold,
# or 35 being run before this rescue was added).
#
# Assumes obj$celltype is populated — true for any object produced by the
# current 11_annotate_celltypes.R. Conservative: only promotes non-NK → NK,
# never demotes.
rescue_nk_clusters <- function(obj, mapping, majority_frac = 0.5) {
  if (!"celltype" %in% colnames(obj@meta.data) ||
      !"knn.leiden.cluster" %in% colnames(obj@meta.data)) {
    return(mapping)
  }

  leiden_ids <- as.character(obj$knn.leiden.cluster)
  cell_ct    <- as.character(obj$celltype)
  n_rescued  <- 0

  for (i in seq_len(nrow(mapping))) {
    if (mapping$celltype_broad[i] == "NK") next
    cl        <- as.character(mapping$cluster[i])
    mask      <- leiden_ids == cl
    n_cluster <- sum(mask)
    if (n_cluster == 0) next

    nk_frac <- sum(!is.na(cell_ct[mask]) & cell_ct[mask] == "NK") / n_cluster
    if (nk_frac >= majority_frac) {
      log_message(sprintf(
        "  NK rescue: cluster %s (%s -> NK, %.0f%% cells marked NK in obj$celltype)",
        cl, mapping$celltype[i], nk_frac * 100))
      mapping$celltype[i]       <- "NK"
      mapping$celltype_broad[i] <- "NK"
      n_rescued <- n_rescued + 1
    }
  }

  if (n_rescued > 0) {
    log_message(sprintf("NK rescue: reclassified %d clusters", n_rescued))
  }
  mapping
}

# --- Map fine labels to a broad lineage (for ordering) -----------------------
.label_to_broad <- function(label) {
  if (is.na(label) || label == "") return("Other")
  lbl <- tolower(label)
  if (grepl("cd4|cd8|\\bt cell|treg|mait|gdt|dnt|tfh", lbl)) return("T cell")
  if (grepl("\\bnk\\b", lbl))                                  return("NK")
  if (grepl("\\bb cell|\\bb naive|\\bb mem|\\bb inter", lbl))  return("B cell")
  if (grepl("plasma", lbl))                                     return("Plasma")
  if (grepl("mono|macro|\\bcd14\\b|\\bcd16\\b", lbl))          return("Monocyte/Macrophage")
  if (grepl("dc|dendritic", lbl))                               return("DC")
  if (grepl("neutro|granulo|\\bmast\\b|basoph|eosinoph", lbl))  return("Granulocyte")
  if (grepl("platelet|megakaryocyte", lbl))                     return("Platelet")
  if (grepl("eryth|\\brbc\\b", lbl))                            return("Eryth")
  if (grepl("mixed", lbl))                                      return("Mixed")
  "Other"
}

# ==============================================================================
# MAIN: merge_clusters_by_celltype
# ==============================================================================
merge_clusters_by_celltype <- function(cfg) {
  obj_path <- file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds")
  map_path <- file.path(cfg$paths$results_tables, "cluster_celltype_mapping.csv")

  if (!file.exists(obj_path)) {
    log_message("Integrated object not found, skipping cluster merge.")
    return(invisible(TRUE))
  }
  if (!file.exists(map_path)) {
    log_message("Cluster-celltype mapping not found, skipping cluster merge.")
    return(invisible(TRUE))
  }

  log_message("Applying cell type labels to clusters...")

  obj <- readRDS(obj_path)
  mapping <- read.csv(map_path, stringsAsFactors = FALSE)

  # --- NK rescue (safety net) ---
  mapping <- rescue_nk_clusters(obj, mapping)

  # --- Label cleanup on the mapping AND per-cell obj$celltype ----------------
  # These renames must be applied to BOTH sources consistently. 35 writes
  # obj$celltype with the raw labels ("CD4 TCM", "Plasmablast", etc.); if 37
  # only renames the mapping, the per-cell override block below sees the
  # remaining obj$celltype values as mismatches and spawns synthetic display
  # clusters for every renamed label (creating ghost "CD4 TCM" / "Plasmablast"
  # display groups with all the cells, while the intended merged groups get
  # zero cells).
  #
  # Biological rationale:
  #   - CD4 Naive vs CD4 TCM: not reliably distinguishable at scRNA-seq
  #     resolution (overlapping transcriptional programs). Collapse to one
  #     "CD4 T Naive/CM" display label.
  #   - Plasmablast: in this project these are memory B cells, not
  #     antibody-secreting plasma cells (true Plasma is a separate population).
  normalize_display_ct <- function(ct) {
    ct <- as.character(ct)
    ct[ct %in% c("CD4 Naive", "CD4 TCM")]         <- "CD4 T Naive/CM"
    ct[ct == "Plasmablast"]                        <- "B memory"
    # CD14 and CD16 monocytes share the "Monocyte" header in the final display.
    # We keep the CD14/CD16 distinction in the raw annotation (obj$celltype at
    # save time via 35), but for display purposes they collapse.
    ct[ct %in% c("CD14 Mono", "CD16 Mono", "CD4 CD14 Mono")] <- "Monocyte"
    ct[grepl("^Mixed", ct)]                        <- "Mixed"
    ct[is.na(ct) | ct == ""]                       <- "Unresolved"
    ct
  }

  mapping$celltype <- normalize_display_ct(mapping$celltype)

  if ("celltype" %in% colnames(obj@meta.data)) {
    was_plasmablast <- !is.na(obj$celltype) & as.character(obj$celltype) == "Plasmablast"
    obj$celltype    <- normalize_display_ct(obj$celltype)
    # Plasmablast was broad "Plasma"; B memory is broad "B cell" — keep the
    # per-cell broad label consistent for just those cells.
    if (any(was_plasmablast) && "celltype_broad" %in% colnames(obj@meta.data)) {
      obj$celltype_broad[was_plasmablast] <- "B cell"
    }
  }

  # --- Filter out ghost + tiny clusters ---
  # 0-cell clusters: SNN refinement can create cluster IDs that end up empty.
  # Tiny clusters (< min_cluster_size cells): likely doublets or outliers that
  # shouldn't get their own display label. Their cells fall through to the
  # fallback path below (labeled as the fallback celltype).
  min_cluster_size <- cfg$merge$min_cluster_size %||% 5
  leiden_ids_all   <- as.character(obj$knn.leiden.cluster)
  cluster_sizes    <- table(leiden_ids_all)

  ghost_rows <- !(as.character(mapping$cluster) %in% names(cluster_sizes))
  if (any(ghost_rows)) {
    log_message(sprintf("Removing %d ghost clusters with 0 cells: %s",
                        sum(ghost_rows),
                        paste(mapping$cluster[ghost_rows], collapse = ", ")))
    mapping <- mapping[!ghost_rows, , drop = FALSE]
  }

  tiny_sizes <- cluster_sizes[as.character(mapping$cluster)]
  tiny_rows  <- as.integer(tiny_sizes) < min_cluster_size
  if (any(tiny_rows)) {
    tiny_info <- sprintf("%s(%s=%d)", mapping$cluster[tiny_rows],
                         mapping$celltype[tiny_rows],
                         as.integer(tiny_sizes[tiny_rows]))
    log_message(sprintf(
      "Removing %d clusters below min_cluster_size=%d (cells reassigned via fallback): %s",
      sum(tiny_rows), min_cluster_size, paste(tiny_info, collapse = ", ")))
    mapping <- mapping[!tiny_rows, , drop = FALSE]
  }

  # --- Disambiguate duplicate cell type names across clusters ---
  # If two clusters both map to "CD4 TCM", they need distinct numbered labels.
  # Sort clusters by broad lineage, then alphabetically within, then assign numbers.
  broad_order <- c("T cell", "NK", "B cell", "Plasma",
                    "Monocyte/Macrophage", "DC", "Granulocyte",
                    "Platelet", "Eryth", "Other", "Mixed")

  mapping <- mapping %>%
    mutate(
      broad_for_sort = sapply(celltype, .label_to_broad),
      broad_for_sort = case_when(
        celltype == "Plasma"     ~ "Plasma",
        celltype == "Mixed"      ~ "Mixed",
        celltype == "Unresolved" ~ "Other",
        TRUE                     ~ broad_for_sort
      ),
      broad_rank = match(broad_for_sort, broad_order, nomatch = 99L)
    ) %>%
    arrange(broad_rank, celltype)

  # Assign display numbers: clusters with the same `celltype` string share a
  # number. After the arrange(broad_rank, celltype) above, identical labels
  # are contiguous, so match() on unique(celltype) yields a clean 1..K
  # numbering where K is the number of distinct cell types.
  mapping$display_number <- match(mapping$celltype, unique(mapping$celltype))
  mapping$display_label  <- paste0(mapping$display_number, ": ", mapping$celltype)

  # --- Direct 1:1 label assignment from annotation map ---
  label_lookup <- setNames(mapping$display_label, as.character(mapping$cluster))
  leiden_ids <- as.character(obj$knn.leiden.cluster)
  merged_labels <- unname(label_lookup[leiden_ids])

  # --- Validate: ensure no cells have NA labels ---
  # Unmapped cells get assigned to the existing "Monocyte" display cluster
  # when one is already present (so they share its display number); otherwise
  # a new display number is appended.
  n_unmapped <- sum(is.na(merged_labels))
  if (n_unmapped > 0) {
    missing_ids  <- setdiff(unique(leiden_ids), as.character(mapping$cluster))
    fallback_ct  <- "Monocyte"
    existing_row <- which(mapping$celltype == fallback_ct)
    if (length(existing_row) >= 1) {
      fallback_num   <- mapping$display_number[existing_row[1]]
      fallback_label <- mapping$display_label[existing_row[1]]
    } else {
      fallback_num   <- max(mapping$display_number) + 1L
      fallback_label <- paste0(fallback_num, ": ", fallback_ct)
    }
    log_message(sprintf("WARNING: %d cells (%d Leiden clusters: %s) not in annotation map. Labeling as '%s'.",
                        n_unmapped, length(missing_ids), paste(missing_ids, collapse = ", "), fallback_label))
    merged_labels[is.na(merged_labels)] <- fallback_label
    # Append to mapping so diagnostics table is complete
    for (mid in missing_ids) {
      mapping <- rbind(mapping, data.frame(
        cluster = mid, celltype = fallback_ct, celltype_broad = "Monocyte/Macrophage",
        vdj_confidence = NA_real_, top_score = NA_real_, n_ref_agree = NA_integer_,
        n_refs_used = NA_integer_, broad_for_sort = "Monocyte/Macrophage",
        broad_rank = match("Monocyte/Macrophage", broad_order, nomatch = 99L),
        display_number = fallback_num, display_label = fallback_label,
        stringsAsFactors = FALSE
      ))
    }
  }

  # --- Per-cell celltype overrides (from 11_annotate_celltypes.R recovery) ----
  # The per-cell NK subcluster recovery in 35 sets obj$celltype to a different
  # label than the cluster's consensus for specific barcodes (e.g. NK cells
  # inside a CD8-majority cluster). The cluster-ID lookup above would snap
  # those cells back to the cluster label, losing the recovery. Detect any
  # such overrides and reassign them to the display cluster that matches
  # their per-cell celltype — creating a synthetic display cluster if none
  # exists.
  if ("celltype" %in% colnames(obj@meta.data)) {
    cell_ct           <- as.character(obj$celltype)
    cluster_ct_lookup <- setNames(mapping$celltype, as.character(mapping$cluster))
    cluster_cts       <- unname(cluster_ct_lookup[leiden_ids])
    # Skip cells whose per-cell label is a fallback ("Unresolved") rather than
    # a genuine per-cell annotation — 35 writes "Unresolved" for cells whose
    # cluster isn't in annotation_map, and 37's own fallback has already given
    # those cells a sensible cluster-level label. Also skip when labels agree.
    override <- !is.na(cell_ct) & cell_ct != "" &
                cell_ct != "Unresolved" &
                !is.na(cluster_cts) & cell_ct != cluster_cts

    if (any(override)) {
      for (target_ct in unique(cell_ct[override])) {
        mask         <- override & cell_ct == target_ct
        target_rows  <- which(mapping$celltype == target_ct)
        if (length(target_rows) >= 1) {
          # Reassign to the first existing display cluster with this celltype
          target_label <- mapping$display_label[target_rows[1]]
          merged_labels[mask] <- target_label
          log_message(sprintf(
            "  Per-cell override: %d cells reassigned from cluster-level label to '%s'",
            sum(mask), target_label))
        } else {
          # No existing display cluster for this celltype — create a synthetic
          # row so the factor levels / diagnostics table stay consistent.
          target_broad <- .label_to_broad(target_ct)
          new_num   <- max(mapping$display_number) + 1L
          new_label <- paste0(new_num, ": ", target_ct)
          merged_labels[mask] <- new_label
          mapping <- rbind(mapping, data.frame(
            cluster        = paste0("recovered_", target_ct),
            celltype       = target_ct,
            celltype_broad = target_broad,
            vdj_confidence = NA_real_,
            top_score      = NA_real_,
            n_ref_agree    = NA_integer_,
            n_refs_used    = NA_integer_,
            broad_for_sort = target_broad,
            broad_rank     = match(target_broad, broad_order, nomatch = 99L),
            display_number = new_num,
            display_label  = new_label,
            stringsAsFactors = FALSE
          ))
          log_message(sprintf(
            "  Per-cell override: created synthetic display cluster '%s' for %d recovered cells",
            new_label, sum(mask)))
        }
      }
    }
  }

  # --- Post-hoc display-size filter ----------------------------------------
  # The earlier min_cluster_size filter operated on raw Leiden cluster sizes,
  # but the per-cell override block can shrink a display cluster by pulling
  # cells out (e.g. a Leiden cluster labeled "gdT" losing almost all its
  # cells to NK/CD8 via per-cell recovery, leaving a singleton gdT display).
  # Count final per-display-label sizes and reassign cells in too-small
  # display clusters to the fallback label.
  final_display_counts <- table(merged_labels)
  tiny_display <- names(final_display_counts)[final_display_counts < min_cluster_size]
  if (length(tiny_display) > 0) {
    fallback_for_tiny <- if (exists("fallback_label", inherits = FALSE)) {
      fallback_label
    } else {
      # Prefer an existing Monocyte display cluster as fallback; otherwise
      # use the first display cluster in the ordered mapping.
      mono_rows <- which(mapping$celltype == "Monocyte")
      if (length(mono_rows) >= 1) mapping$display_label[mono_rows[1]]
      else                        mapping$display_label[1]
    }
    to_reassign <- merged_labels %in% tiny_display
    log_message(sprintf(
      "Post-hoc display filter: reassigning %d cells from %d tiny display cluster(s) (<%d cells) to '%s': %s",
      sum(to_reassign), length(tiny_display), min_cluster_size,
      fallback_for_tiny, paste(tiny_display, collapse = ", ")))
    merged_labels[to_reassign] <- fallback_for_tiny
    # Drop the now-empty display labels from mapping so factor levels are clean
    mapping <- mapping[!mapping$display_label %in% tiny_display, , drop = FALSE]
  }

  # --- Set factor levels in display order ---
  # mapping$display_label may now contain duplicates (multiple Leiden clusters
  # sharing a display number for the same celltype). factor() levels must be
  # unique — unique() preserves first-occurrence order, which is the correct
  # display order given the preceding arrange().
  ordered_levels <- unique(mapping$display_label)

  obj$merged.celltype.cluster <- factor(merged_labels, levels = ordered_levels)
  Idents(obj) <- "merged.celltype.cluster"

  # --- Also store the plain celltype label (without number) for downstream use ---
  # Prefer per-cell obj$celltype when available (honors 35's per-cell recovery);
  # fall back to the cluster-level celltype lookup for cells without a per-cell
  # label assigned.
  celltype_lookup <- setNames(mapping$celltype, as.character(mapping$cluster))
  ct_labels <- unname(celltype_lookup[leiden_ids])
  if ("celltype" %in% colnames(obj@meta.data)) {
    cell_ct   <- as.character(obj$celltype)
    has_cell  <- !is.na(cell_ct) & cell_ct != ""
    ct_labels <- ifelse(has_cell, cell_ct, ct_labels)
  }
  ct_labels[is.na(ct_labels)] <- "Monocyte"
  obj$celltype_cluster <- ct_labels

  # --- Save diagnostic table ---
  # Per-Leiden-cluster rows so users can see which raw clusters feed each
  # display label. Cell counts are per Leiden cluster (not the display-label
  # total) to preserve that resolution.
  leiden_cell_counts <- as.data.frame(table(leiden_ids)) %>%
    setNames(c("cluster", "n_cells"))
  leiden_cell_counts$cluster <- as.character(leiden_cell_counts$cluster)

  diag_df <- mapping %>%
    dplyr::select(display_number, display_label, cluster, celltype, celltype_broad, broad_for_sort) %>%
    rename(merged_number = display_number, merged_label = display_label,
           leiden_cluster = cluster, broad_lineage = broad_for_sort)
  diag_df <- merge(diag_df, leiden_cell_counts,
                   by.x = "leiden_cluster", by.y = "cluster", all.x = TRUE)
  diag_df$n_cells[is.na(diag_df$n_cells)] <- 0
  diag_df <- diag_df %>% arrange(merged_number, leiden_cluster)

  ensure_dir(cfg$paths$results_tables)
  write.csv(diag_df,
            file.path(cfg$paths$results_tables, "merged_cluster_mapping.csv"),
            row.names = FALSE)

  # Log one line per display label with the combined cell count.
  display_summary <- as.data.frame(table(obj$merged.celltype.cluster)) %>%
    setNames(c("merged_label", "n_cells")) %>%
    arrange(match(merged_label, ordered_levels))
  log_message(sprintf("Merged -> %d display clusters (from %d Leiden clusters):",
                      nrow(display_summary), nrow(mapping)))
  for (i in seq_len(nrow(display_summary))) {
    log_message(sprintf("  %-25s %7d cells",
                        display_summary$merged_label[i], display_summary$n_cells[i]))
  }

  saveRDS(obj, obj_path)
  log_message("Saved object with merged.celltype.cluster identities.")

  invisible(TRUE)
}
