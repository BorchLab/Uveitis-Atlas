#!/usr/bin/env Rscript

# Run from the repository root (the directory containing R/ and config/).

# Centralized package management — installs missing packages and loads all
source("R/00_setup_packages.R")
load_all_packages()

# --- Setup + ingest --------------------------------------------------------
source("R/01_setup_utils.R")
source("R/02_ingest_data.R")
source("R/03_annotate_AIRR.R")

# --- Integration + annotation (full atlas) ---------------------------------
source("R/10_integrate_full.R")
source("R/11_annotate_celltypes.R")
source("R/12_merge_clusters.R")

# --- QC + read-only diagnostics --------------------------------------------
source("R/14_repertoire_qc.R")
source("R/15_lens_qc.R")
source("R/16_compartment_lineage_gate.R")
source("R/17_lineage_validation.R")
source("R/18_lineage_purity_audit.R")

# --- Eye subset + compartments ---------------------------------------------
source("R/20_subset_eye.R")
source("R/21_integrate_eye.R")
source("R/22_subset_compartments.R")
source("R/23_substate_labels.R")

# --- Markers / DGE / escape / differential abundance -----------------------
source("R/30_markers.R")
source("R/32_escape.R")
source("R/40_milo.R")
source("R/41_composition.R")

# --- Compartment PCA + cross-compartment ligand-receptor -------------------
source("R/45_compartment_pca.R")
source("R/46_cross_compartment_bridge.R")
source("R/47_liana_myeloid_tcell.R")
source("R/48_nichenet_myeloid_tcell.R")
source("R/49_liana_bcell.R")
source("R/49b_liana_supplemental_table.R")

# --- Repertoire + BCR lineage ----------------------------------------------
source("R/51_eye_blood_repertoire.R")
source("R/52_bcr_lineage.R")
source("R/53_public_bcr_clones.R")
source("R/54_bcell_paired_eye_blood_metrics.R")
source("R/55_tcell_paired_eye_blood_metrics.R")
source("R/56_bcell_lineage_trees.R")
source("R/57_bcell_cross_compartment_bridge.R")
source("R/58_bcell_lineage_architecture.R")

# --- TCR motif / distance / signatures -------------------------------------
source("R/64_immgliph.R")
source("R/65_immlynx_tcrdist.R")
source("R/66_immlynx_olga_sonnia.R")
source("R/67_gliph_tcrdist_correlation.R")
source("R/68_gliph_network_metrics.R")
source("R/69_vdjdb_annotation.R")
source("R/70_tcr_genex_signatures.R")
source("R/71_gliph_motif_overlap.R")
source("R/79_tcr_sharpened_bridge.R")
source("R/80_novel_tcr_scoring.R")

# --- Visualization (shared helpers first, then per-target panels) ----------
source("R/81_viz_compartment_helpers.R")
source("R/82_viz_dispatch.R")
source("R/83_viz_full.R")
source("R/84_viz_eye.R")
source("R/85_viz_myeloid.R")
source("R/86_viz_bcell.R")
source("R/88_viz_tcell.R")
source("R/89_viz_tcr_advanced.R")
source("R/90_viz_tcr_sharing_supp.R")
source("R/91_viz_cross_compartment.R")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) > 0) args[1] else "config/config.yml"
# Force UTF-8 locale + use yaml.load on a UTF-8-read text body so config keys
# resolve identically regardless of the C/system locale at startup. Without
# this, readLines on a config containing em-dashes / >= / etc. silently
# truncated the parse, leaving cfg$steps$* as NULL and making every step
# silently skip (encountered 2026-05-22 during Fig 5 production).
try(Sys.setlocale("LC_ALL", "en_US.UTF-8"), silent = TRUE)
cfg_text <- readLines(cfg_path, encoding = "UTF-8", warn = FALSE)
cfg <- yaml::yaml.load(paste(cfg_text, collapse = "\n"))

set_global_seed(cfg$seed)
options(future.globals.maxSize = cfg$globals$future_maxbytes %||% 8e9)

log_message("Pipeline start")
ensure_dir(cfg$paths$results_objects)
ensure_dir(cfg$paths$results_tables)
ensure_dir(cfg$paths$viz_dir)

# Phase 1: Core scRNA-seq (full object)
if (isTRUE(cfg$steps$ingest_data))      ingest_data(cfg)
if (isTRUE(cfg$steps$airr_annotation))  annotate_AIRR(cfg)
if (isTRUE(cfg$steps$integration))      integrate_fastmnn(cfg)
if (isTRUE(cfg$steps$annotation))       annotate_celltypes(cfg, target = "all")
if (isTRUE(cfg$steps$merge_clusters))   merge_clusters_by_celltype(cfg)
if (isTRUE(cfg$steps$markers))          run_markers(cfg, target = "all")
if (isTRUE(cfg$steps$dge))              run_dge(cfg, target = "all")
if (isTRUE(cfg$steps$escape))           run_escape_ssgsea(cfg)
if (isTRUE(cfg$steps$escape))           run_escape_differential(cfg)

# Phase 2: Eye-focus subset + re-integrate + re-annotate
if (isTRUE(cfg$steps$eye_subset))       subset_eye(cfg)
if (isTRUE(cfg$steps$eye_reintegrate))  reintegrate_eye(cfg)
if (isTRUE(cfg$steps$eye_markers))      run_markers(cfg, target = "eye")
if (isTRUE(cfg$steps$eye_annotate))     annotate_celltypes(cfg, target = "eye")
if (isTRUE(cfg$steps$eye_dge))          run_dge(cfg, target = "eye")

# Cross-tissue + within-eye downstream
if (isTRUE(cfg$steps$milo))              run_milo_da(cfg, target = "all")
if (isTRUE(cfg$steps$composition))       run_composition_testing(cfg, target = "all")
if (isTRUE(cfg$steps$eye_milo))          run_milo_da(cfg, target = "eye")
if (isTRUE(cfg$steps$eye_composition))   run_composition_testing(cfg, target = "eye")

# Phase 3: Compartment-specific deep dives within the eye sub-atlas (F3-F5)
# Builds myeloid / bcell / tcell sub-atlases by re-integrating each cluster
# group, then runs the same downstream stack as the eye branch on each.
if (isTRUE(cfg$steps$lens_ucell_diagnostic)) run_lens_ucell_diagnostic(cfg)
if (isTRUE(cfg$steps$ctaa_filter_diagnostic)) run_ctaa_filter_diagnostic(cfg)
if (isTRUE(cfg$steps$clone_definition_sensitivity)) run_clone_definition_sensitivity(cfg)
if (isTRUE(cfg$steps$lineage_validation)) run_lineage_validation(cfg)
if (isTRUE(cfg$steps$pseudobulk_N_audit)) {
  for (tgt in c("myeloid", "bcell", "tcell")) run_pseudobulk_N_audit(cfg, target = tgt)
}
if (isTRUE(cfg$steps$lineage_purity_audit)) run_lineage_purity_audit(cfg)

if (isTRUE(cfg$steps$compartment_subset)) subset_compartments(cfg)
for (cmp in c("myeloid", "bcell", "tcell")) {
  if (isTRUE(cfg$steps$compartment_markers))       run_markers(cfg, target = cmp)
  if (isTRUE(cfg$steps$compartment_dge))           run_dge(cfg, target = cmp)
  if (isTRUE(cfg$steps$compartment_escape))        run_escape_differential(cfg, target = cmp)
  if (isTRUE(cfg$steps$compartment_escape))        run_escape_etiology_breakdown(cfg, target = cmp)
  if (isTRUE(cfg$steps$compartment_custom_modules))
    run_escape_custom_modules(cfg, target = cmp)
  if (isTRUE(cfg$steps$compartment_composition))   run_composition_testing(cfg, target = cmp)
  if (isTRUE(cfg$steps$compartment_milo))          run_milo_da(cfg, target = cmp)
}

# Phase 4: Cross-compartment workflows
#   joint_substate         — stamps `substate_joint` on the eye object (still
#                            used for multi-compartment viz overlays)
#   compartment_pca        — per-substate pseudobulk PCA per compartment;
#                            writes pca_subject_scores.csv etc.
#   cross_compartment_bridge — subject-level myeloid PC1 vs T cell PC1, with
#                            bootstrap CI + permutation null (panel F)
#   liana_myeloid_tcell    — rank-aggregated LR map per etiology (panel G)
#   nichenet_myeloid_tcell — ligand-activity prediction (panel H, OFF by
#                            default; flip on after LIANA review)
#   viz_cross_compartment  — panels F-H
if (isTRUE(cfg$steps$joint_substate))            build_joint_substate_labels(cfg)
if (isTRUE(cfg$steps$compartment_pca)) {
  for (t in c("myeloid", "tcell", "bcell"))     run_compartment_pca(cfg, target = t)
}
if (isTRUE(cfg$steps$cross_compartment_bridge))  run_cross_compartment_pc1_bridge(cfg)
if (isTRUE(cfg$steps$cross_compartment_bridge))  run_bcell_cross_compartment_bridge(cfg)
if (isTRUE(cfg$steps$liana_myeloid_tcell))       run_liana_myeloid_tcell(cfg)
if (isTRUE(cfg$steps$nichenet_myeloid_tcell))    run_nichenet_myeloid_to_tcell(cfg)
if (isTRUE(cfg$steps$viz_cross_compartment))     run_visualizations(cfg, target = "cross_compartment")

# Phase 5: Repertoire (cross-tissue, full object only)
if (isTRUE(cfg$steps$bcr_lineage))      run_bcr_lineage(cfg)
if (isTRUE(cfg$steps$bcell_lineage_architecture))
  run_bcell_lineage_architecture(cfg)
if (isTRUE(cfg$steps$bcell_baseline_selection)) {
  source("R/59_bcell_baseline_selection.R")
  run_bcell_baseline_selection(cfg)
}
if (isTRUE(cfg$steps$bcell_architecture_viz)) {
  source("R/92_viz_bcell_architecture.R")
  run_viz_bcell_architecture(cfg)
  source("R/93_viz_bcell_architecture_supp.R")
  run_viz_bcell_architecture_supp(cfg)
}
if (isTRUE(cfg$steps_fig6$check_architecture) && exists("check_fig6_architecture_schema"))
  check_fig6_architecture_schema(cfg)
if (isTRUE(cfg$steps$eye_blood_repertoire) && exists("run_eye_blood_repertoire"))
  run_eye_blood_repertoire(cfg)

# Phase 6: TCR motif / distance / generation analysis (Viral vs NIU intraocular)
# All three steps are independent and config-gated. TRB filter (length<60,
# single chain) is applied per-step.
if (isTRUE(cfg$steps$gliph) && exists("run_immgliph"))
  run_immgliph(cfg)
if (isTRUE(cfg$steps$tcrdist) && exists("run_immlynx_tcrdist"))
  run_immlynx_tcrdist(cfg)
if (isTRUE(cfg$steps$olga) && exists("run_immlynx_olga"))
  run_immlynx_olga(cfg)
# Bootstrap the tcrdist pairwise-beta matrix from the saved SCE if a gliph
# distance step needs it and it is not already on disk (folded in from the
# former scripts/run_tcr_advanced.R).
if ((isTRUE(cfg$steps$gliph_tcrdist) || isTRUE(cfg$steps$gliph_network)) &&
    !file.exists(file.path(cfg$paths$results_objects, "tcrdist_pw_beta.rds")))
  source("scripts/generate_tcrdist_pw_beta.R")
if (isTRUE(cfg$steps$gliph_tcrdist) && exists("run_gliph_tcrdist_correlation"))
  run_gliph_tcrdist_correlation(cfg)
if (isTRUE(cfg$steps$gliph_network) && exists("run_gliph_network_metrics"))
  run_gliph_network_metrics(cfg)
if (isTRUE(cfg$steps$vdjdb) && exists("run_vdjdb_annotation"))
  run_vdjdb_annotation(cfg)
if (isTRUE(cfg$steps$tcr_genex_signatures) && exists("run_tcr_genex_signatures"))
  run_tcr_genex_signatures(cfg)
if (isTRUE(cfg$steps$gliph_motif_overlap) && exists("run_gliph_motif_overlap"))
  run_gliph_motif_overlap(cfg)
if (isTRUE(cfg$steps$tcr_sharpened_bridge) && exists("run_tcr_sharpened_bridge"))
  run_tcr_sharpened_bridge(cfg)
if (isTRUE(cfg$steps$novel_tcr_discovery) && exists("run_novel_tcr_scoring"))
  run_novel_tcr_scoring(cfg)
if (isTRUE(cfg$steps$viz_tcr_advanced) && exists("run_visualizations_tcr_advanced"))
  run_visualizations_tcr_advanced(cfg)
if (isTRUE(cfg$steps$compartment_visualize)) {
  for (cmp in c("myeloid", "bcell", "tcell")) {
    run_visualizations(cfg, target = cmp)
  }
}

# Phase 7: Visualization (per target)
if (isTRUE(cfg$steps$visualize))        run_visualizations(cfg, target = "all")
if (isTRUE(cfg$steps$eye_visualize))    run_visualizations(cfg, target = "eye")

# Phase 8: Figure 6 (B cell) 
if (isTRUE(cfg$steps_fig6$bcr_public_clones) && exists("run_public_bcr_clones"))
  run_public_bcr_clones(cfg)
if (isTRUE(cfg$steps_fig6$bcell_paired_metrics) &&
    exists("run_bcell_paired_eye_blood_metrics"))
  run_bcell_paired_eye_blood_metrics(cfg)
if (isTRUE(cfg$steps_fig6$tcell_paired_metrics) &&
    exists("run_tcell_paired_eye_blood_metrics"))
  run_tcell_paired_eye_blood_metrics(cfg)
if (isTRUE(cfg$steps_fig6$bcell_lineage_trees)) {
  if (exists("run_bcell_lineage_trees"))
    tryCatch(run_bcell_lineage_trees(cfg),
             error = function(e)
               log_message("  bcell_lineage_trees failed: ",
                           conditionMessage(e)))
  if (exists("run_bcell_clonal_sharing_circos"))
    tryCatch(run_bcell_clonal_sharing_circos(cfg),
             error = function(e)
               log_message("  bcell_clonal_sharing_circos failed: ",
                           conditionMessage(e)))
}
if (isTRUE(cfg$steps_fig6$liana_bcell) && exists("run_liana_bcell"))
  run_liana_bcell(cfg)
if (isTRUE(cfg$steps_fig6$liana_supplemental_table) &&
    exists("run_liana_supplemental_table"))
  run_liana_supplemental_table(cfg)
if (isTRUE(cfg$steps_fig6$viz_bcell_fig6) && exists("run_visualizations_bcell_fig6"))
  run_visualizations_bcell_fig6(cfg)
if (isTRUE(cfg$steps_fig6$viz_tcr_sharing_supp) &&
    exists("run_visualizations_tcr_sharing_supp"))
  run_visualizations_tcr_sharing_supp(cfg)

log_message("Pipeline complete")
