#!/usr/bin/env Rscript
# Re-derive tcrdist pw_beta from the saved SCE in ImmLynxTcrdistResults.rds
# without re-running the ESM-2 embedding step. Used to bootstrap the new
# R/67/68 modules when ImmLynxTcrdistResults.rds was produced by an older
# version of R/65 (before the pw_beta saver was added).
#
# Strategy:
#   The full intraocular SCE can exceed tcrdist3's 10k-clone safety threshold
#   (in this cohort ~28k unique TCRs), which causes runTCRdist to refuse the
#   dense pw_beta computation. R/67 and R/68 only need distances among GLIPH
#   cluster members and their same-Subject background — typically <5k cells.
#   So we subset the SCE to (a) all cells whose CDR3-beta appears in any
#   GLIPH convergence group AND (b) any cell from the same Subjects as those
#   GLIPH members. This keeps the matrix small and preserves the
#   Subject-stratified permutation null used in module 67.
#
# Run:   Rscript scripts/generate_tcrdist_pw_beta.R
# Out:   outputs/objects/tcrdist_pw_beta.rds

# Run from the repository root (the directory containing R/ and config/).

suppressPackageStartupMessages({
  source("R/00_setup_packages.R")
  load_all_packages()
})
source("R/01_setup_utils.R")
# .tcrdist_to_cell_matrix() expands tcrdist3's deduplicated clone-level
# matrix back to a cell-level matrix indexed by barcode (see R/65 for the
# rationale + verification).
source("R/65_immlynx_tcrdist.R")

td_rds    <- "outputs/objects/ImmLynxTcrdistResults.rds"
gliph_rds <- "outputs/objects/ImmGLIPHResults.rds"
out_path  <- "outputs/objects/tcrdist_pw_beta.rds"

if (file.exists(out_path)) {
  log_message("tcrdist_pw_beta.rds already exists; nothing to do.")
  quit(status = 0)
}
if (!file.exists(td_rds))    stop("Missing: ", td_rds)
if (!file.exists(gliph_rds)) stop("Missing: ", gliph_rds)

log_message("Loading tcrdist + GLIPH results...")
td_res <- readRDS(td_rds)
gl_res <- readRDS(gliph_rds)
sce    <- td_res$sce
if (is.null(sce)) stop("td_res$sce is NULL; cannot subset.")

cd      <- SummarizedExperiment::colData(sce)
gl_cdr3 <- unique(gl_res$clusters$CDR3b)
log_message("  Total cells in SCE: ", ncol(sce),
            "; GLIPH-cluster CDR3-beta: ", length(gl_cdr3))

# Restrict to GLIPH-cluster member cells. tcrdist3 refuses dense pw_beta
# above 10k clones, and in this cohort the same-subject "background" pool
# is essentially the whole repertoire (>50k cells) so a broader subset
# doesn't help. The Subject-stratified background in module 67 still works
# because it intersects pool with rownames(pw_beta), restricting the null
# distribution to other GLIPH-clustered TCRs from the same subjects — a
# more conservative test than draws from the full repertoire.
in_gliph <- cd$CDR3b %in% gl_cdr3
log_message("  Keeping ", sum(in_gliph), " GLIPH-member cells")
if (sum(in_gliph) < 10L)
  stop("Too few GLIPH-member cells (", sum(in_gliph), ") to compute distances.")

sce_sub <- sce[, in_gliph]
log_message("Running immLynx::runTCRdist on subset (", ncol(sce_sub), " cells)...")
td <- immLynx::runTCRdist(sce_sub, chains = "beta",
                          compute_distances = TRUE,
                          add_to_object = FALSE)

dmat <- .tcrdist_to_cell_matrix(td)
if (is.null(dmat))
  stop("Could not produce pw_beta on GLIPH-only subset (",
       ncol(sce_sub), " cells). Check immLynx / tcrdist install.")

saveRDS(dmat, out_path)
log_message("Saved: ", out_path, " (", nrow(dmat), " x ", ncol(dmat), ")")
