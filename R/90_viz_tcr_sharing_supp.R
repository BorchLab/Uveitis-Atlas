# R/90_viz_tcr_sharing_supp.R
# Supplemental TCR blood<->eye sharing figure (mirrors R/74b for BCR).
# Panels:
#   S-A  TCR overlap per subject (Morisita + Jaccard), NIU vs Viral.
#   S-B  Receptors of interest: top expanded eye TCR clones with VDJdb
#        antigen annotations and GLIPH-cluster overlay.
#   S-C  Substate (T cell cluster) distribution of shared / eye-only /
#        blood-only TCR clones.
#   S-D  Public-across-subjects TCR convergence (reuses the bcell-side
#        helper which is actually compartment-agnostic).
#
# No new compute: all upstream tables are produced by existing modules
# (R/51 -> TCR_*, R/69 -> vdjdb_*, R/68/R/64 -> gliph_*).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

.supp_tcr_overlap <- function(viz_dir) {
  ov_path <- "outputs/tables/repertoire/TCR_eye_blood_overlap.csv"
  if (!file.exists(ov_path)) return(invisible())
  ov <- utils::read.csv(ov_path, stringsAsFactors = FALSE)
  group_col <- if ("phenotype" %in% colnames(ov)) "phenotype" else "etiology"
  d <- ov |> dplyr::filter(.data[[group_col]] %in% c("NIU","Viral"))
  if (nrow(d) == 0L) return(invisible())
  d_long <- d |>
    dplyr::select(subject, group = .data[[group_col]], jaccard, morisita) |>
    tidyr::pivot_longer(c(jaccard, morisita), names_to = "metric",
                        values_to = "value") |>
    dplyr::filter(!is.na(value))
  p <- ggplot(d_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.75, size = 1.5) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_fill_manual(values = ETIOLOGY_GROUP_COLORS) +
    labs(title = "TCR eye<->blood overlap per subject",
         x = NULL, y = "Overlap value") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  save_pdf_png(p, file.path(viz_dir, "supp_tcr_overlap_morisita_jaccard"),
               w = 8, h = 5)
}

# Build a clone-level antigen + gliph annotation table by mapping the
# top-expanded clone_id (CTstrict) -> cells (barcodes) -> vdjdb species and
# GLIPH cluster. Returns NULL when the tcell object or annotation tables
# are missing.
.tcr_clone_annotation_table <- function(cfg) {
  paths_t <- get_target_paths(cfg, "tcell")
  obj_path <- file.path(paths_t$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("  supp TCR ROI: tcell object missing at ", obj_path)
    return(NULL)
  }
  obj <- readRDS(obj_path)
  meta <- obj@meta.data
  meta$barcode <- rownames(meta)
  if (!"CTstrict" %in% colnames(meta)) return(NULL)

  top_path <- "outputs/tables/repertoire/TCR_top_expanded_eye.csv"
  if (!file.exists(top_path)) return(NULL)
  top <- utils::read.csv(top_path, stringsAsFactors = FALSE)

  # Restrict tcell metadata to (subject, clone_id) pairs in top
  keep <- top |> dplyr::distinct(subject, clone_id)
  cells <- meta |>
    dplyr::filter(!is.na(CTstrict)) |>
    dplyr::semi_join(keep |>
                     dplyr::rename(Subject = subject, CTstrict = clone_id),
                     by = c("Subject","CTstrict")) |>
    dplyr::select(barcode, Subject, Tissue_1, CTstrict, CTaa,
                  knn.leiden.cluster, Phenotype_2)

  # VDJdb antigen species at the barcode level
  vdj_path <- "outputs/tables/repertoire/vdjdb_annotations.csv"
  if (file.exists(vdj_path)) {
    vdj <- utils::read.csv(vdj_path, stringsAsFactors = FALSE)
    cells <- dplyr::left_join(cells,
      vdj |> dplyr::select(barcode, antigen_species, panel_tier),
      by = "barcode")
  } else {
    cells$antigen_species <- NA_character_
    cells$panel_tier      <- NA_character_
  }

  # GLIPH cluster via CDR3b. gliph_clusters.csv columns: cluster_id, CDR3b,
  # TRBV, patient, motif, type.
  gliph_path <- "outputs/tables/repertoire/gliph_clusters.csv"
  if (file.exists(gliph_path)) {
    gl <- utils::read.csv(gliph_path, stringsAsFactors = FALSE)
    # Extract CDR3b from CTaa ("TRA_TRB"): take portion after "_"
    cells$cdr3b <- vapply(strsplit(as.character(cells$CTaa), "_"),
                          function(p) if (length(p) >= 2) p[2] else NA_character_,
                          character(1))
    cells <- dplyr::left_join(cells,
      gl |> dplyr::select(CDR3b, gliph_cluster = cluster_id, gliph_motif = motif),
      by = c("cdr3b" = "CDR3b"))
  }

  # Collapse to one row per (subject, clone_id): vote for antigen, gliph
  per_clone <- cells |>
    dplyr::group_by(Subject, CTstrict) |>
    dplyr::summarise(
      n_cells_total     = dplyr::n(),
      eye_cells         = sum(Tissue_1 == "Eye",   na.rm = TRUE),
      blood_cells       = sum(Tissue_1 == "Blood", na.rm = TRUE),
      antigen_consensus = {
        x <- stats::na.omit(antigen_species)
        if (length(x) == 0) NA_character_ else names(sort(table(x), decreasing = TRUE))[1]
      },
      gliph_clusters    = paste(sort(unique(stats::na.omit(gliph_cluster))), collapse = ","),
      gliph_motifs      = paste(sort(unique(stats::na.omit(gliph_motif))), collapse = ","),
      phenotype_2       = dplyr::first(stats::na.omit(Phenotype_2)),
      .groups = "drop") |>
    dplyr::rename(subject = Subject, clone_id = CTstrict)

  # Join in top's freq columns; top$etiology is already present (Etiology subtype)
  top_joined <- dplyr::left_join(top, per_clone, by = c("subject","clone_id"))
  top_joined
}

.supp_tcr_receptors_of_interest <- function(viz_dir, cfg) {
  d <- .tcr_clone_annotation_table(cfg)
  if (is.null(d) || nrow(d) == 0L) return(invisible())
  # Persist the joined table for the manuscript record
  ensure_dir("outputs/tables/eye/tcell")
  utils::write.csv(d,
    file.path("outputs/tables/eye/tcell", "tcr_receptors_of_interest.csv"),
    row.names = FALSE)

  # Strip plot: n_eye vs n_blood, color by antigen consensus, shape by phenotype
  d_plot <- d |> dplyr::filter(!is.na(antigen_consensus) |
                              found_in_blood | gliph_clusters != "")
  if (nrow(d_plot) == 0L) {
    log_message("  supp TCR ROI: no annotated receptors of interest.")
    return(invisible())
  }
  d_plot$antigen_consensus[is.na(d_plot$antigen_consensus) |
                          d_plot$antigen_consensus == ""] <- "Unannotated"
  # Color from shared palette; fall through to grey for unmapped
  pal <- c(PATHOGEN_COLORS,
           Unannotated = "grey80")
  # Shape on Phenotype_2 (NIU vs Viral); falls back to etiology subtype string
  shape_col <- if ("phenotype_2" %in% colnames(d_plot)) "phenotype_2" else "etiology"
  p <- ggplot(d_plot, aes(x = pmax(n_cells_eye, 1L),
                          y = pmax(n_cells_blood, 1L),
                          color = antigen_consensus,
                          shape = .data[[shape_col]])) +
    geom_jitter(width = 0.15, height = 0.15, size = 2.0, alpha = 0.85) +
    scale_x_log10() + scale_y_log10() +
    scale_color_manual(values = pal, name = "Antigen species") +
    geom_abline(slope = 1, linetype = "dashed", color = "grey60") +
    labs(title = "Top expanded eye TCR clones: receptors of interest",
         subtitle = "Color = VDJdb antigen consensus; shape = phenotype",
         x = "Cells in eye (log scale)",
         y = "Cells in blood (log scale)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  save_pdf_png(p,
    file.path(viz_dir, "supp_tcr_receptors_of_interest_vdjdb"),
    w = 9, h = 7)

  # GLIPH cluster counts for shared receptors (only if gliph_clusters has entries)
  facet_col <- if ("phenotype_2" %in% colnames(d)) "phenotype_2" else "etiology"
  gd <- d |>
    dplyr::filter(!is.na(gliph_clusters), gliph_clusters != "") |>
    tidyr::separate_rows(gliph_clusters, sep = ",") |>
    dplyr::mutate(shared = found_in_blood) |>
    dplyr::count(gliph_clusters, shared, !!rlang::sym(facet_col), name = "n_clones")
  if (nrow(gd) > 0L) {
    p2 <- ggplot(gd,
                 aes(x = forcats::fct_reorder(gliph_clusters, n_clones, sum),
                     y = n_clones,
                     fill = factor(shared, levels = c(TRUE, FALSE),
                                   labels = c("Shared eye+blood","Eye-only")))) +
      geom_col() +
      coord_flip() +
      facet_wrap(stats::as.formula(paste("~", facet_col))) +
      scale_fill_manual(values = c(`Shared eye+blood` = "#1F5C92",
                                   `Eye-only` = "#E2A14A"),
                        name = NULL) +
      labs(title = "GLIPH clusters containing shared eye+blood clones",
           x = "GLIPH cluster id", y = "n clones") +
      theme_classic(base_size = 9)
    save_pdf_png(p2,
      file.path(viz_dir, "supp_tcr_gliph_shared_eye_blood"),
      w = 11, h = max(6, 0.25 * dplyr::n_distinct(gd$gliph_clusters) + 3))
  }
}

.supp_tcr_substate_distribution_shared <- function(viz_dir, cfg) {
  ct_path <- "outputs/tables/repertoire/TCR_top_expanded_eye_celltype.csv"
  if (!file.exists(ct_path)) return(invisible())
  ct <- utils::read.csv(ct_path, stringsAsFactors = FALSE)
  if (nrow(ct) == 0L) return(invisible())
  top_path <- "outputs/tables/repertoire/TCR_top_expanded_eye.csv"
  if (!file.exists(top_path)) return(invisible())
  top <- utils::read.csv(top_path, stringsAsFactors = FALSE)

  shared_keys <- top |>
    dplyr::transmute(subject, clone_id, found_in_blood)
  ct <- dplyr::left_join(ct, shared_keys, by = c("subject","clone_id"))
  ct$share_group <- ifelse(ct$found_in_blood %in% TRUE,
                           "Shared eye+blood", "Eye-only")

  # Map etiology subtype -> NIU/Viral via cfg$etiology_groups
  niu_set   <- as.character(cfg$etiology_groups$niu   %||% character(0))
  viral_set <- as.character(cfg$etiology_groups$viral %||% character(0))
  ct$phenotype_2 <- dplyr::case_when(
    ct$etiology %in% niu_set   ~ "NIU",
    ct$etiology %in% viral_set ~ "Viral",
    TRUE                        ~ NA_character_)

  # Aggregate across subjects within each (share_group, phenotype_2, substate)
  agg <- ct |>
    dplyr::filter(modality == "TCR", substate != "out_of_compartment",
                  phenotype_2 %in% c("NIU","Viral")) |>
    dplyr::group_by(share_group, phenotype_2, substate) |>
    dplyr::summarise(n = sum(n_cells, na.rm = TRUE), .groups = "drop")
  if (nrow(agg) == 0L) return(invisible())
  p <- ggplot(agg, aes(x = share_group, y = n, fill = substate)) +
    geom_col(position = "fill") +
    facet_wrap(~ phenotype_2) +
    scale_fill_viridis_d(option = "viridis", name = "T substate") +
    labs(title = "T-cell substate composition of shared vs eye-only TCRs",
         x = NULL, y = "Fraction of cells") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          plot.title = element_text(face = "bold"))
  save_pdf_png(p,
    file.path(viz_dir, "supp_tcr_substate_distribution_shared"),
    w = 9, h = 6)
}

run_visualizations_tcr_sharing_supp <- function(cfg) {
  viz_dir <- file.path(viz_subdir(get_target_paths(cfg, "tcell"), "repertoire"),
                       "supp_tcr_sharing")
  ensure_dir(viz_dir)
  log_message("=== Supplemental TCR sharing figure ===")
  tryCatch(.supp_tcr_overlap(viz_dir),                 error = function(e)
           log_message("  supp TCR overlap: ", conditionMessage(e)))
  tryCatch(.supp_tcr_receptors_of_interest(viz_dir, cfg), error = function(e)
           log_message("  supp TCR ROI: ", conditionMessage(e)))
  tryCatch(.supp_tcr_substate_distribution_shared(viz_dir, cfg),
           error = function(e)
           log_message("  supp TCR substate: ", conditionMessage(e)))
  log_message("=== Supplemental TCR sharing complete ===")
  invisible(TRUE)
}

# run_visualizations_tcr_sharing_supp renders three TCR sharing supplement
# panels using only tables produced upstream by R/51 (TCR_eye_blood_overlap
# + TCR_top_expanded_eye + TCR_top_expanded_eye_celltype), R/69
# (vdjdb_annotations.csv), and R/64/R/68 (gliph_clusters.csv +
# gliph_cluster_properties.csv). The receptors-of-interest panel joins these
# at the clone level by loading the tcell IntegratedSeuratObject once to
# map clone_id (CTstrict) -> barcodes -> VDJdb antigen and GLIPH cluster.
# Output directory: outputs/viz/eye/tcell/08_repertoire/supp_tcr_sharing/.
