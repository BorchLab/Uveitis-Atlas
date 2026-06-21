# R/00_setup_packages.R
# Centralized package installation and loading for the Eye_Blood_Pairs pipeline.
# Source this file first in run_pipeline.R before any other scripts.

# ---- Package manifest --------------------------------------------------------
# Each entry: list(pkg, source) where source is "cran", "bioc", or "github/repo"
.PACKAGE_MANIFEST <- list(
  # CRAN
  list(pkg = "yaml",        source = "cran"),
  list(pkg = "dplyr",       source = "cran"),
  list(pkg = "stringr",     source = "cran"),
  list(pkg = "tidyr",       source = "cran"),
  list(pkg = "tibble",      source = "cran"),
  list(pkg = "ggplot2",     source = "cran"),
  list(pkg = "patchwork",   source = "cran"),
  list(pkg = "forcats",     source = "cran"),
  list(pkg = "Matrix",      source = "cran"),
  list(pkg = "glue",        source = "cran"),
  list(pkg = "igraph",      source = "cran"),
  list(pkg = "mclust",      source = "cran"),
  list(pkg = "cluster",     source = "cran"),
  list(pkg = "uwot",        source = "cran"),
  list(pkg = "UpSetR",      source = "cran"),
  list(pkg = "viridis",     source = "cran"),

  # Bioconductor
  list(pkg = "Seurat",                  source = "bioc"),
  list(pkg = "SingleCellExperiment",    source = "bioc"),
  list(pkg = "SummarizedExperiment",    source = "bioc"),
  list(pkg = "scDblFinder",             source = "bioc"),
  list(pkg = "batchelor",              source = "bioc"),
  list(pkg = "bluster",                source = "bioc"),
  list(pkg = "leidenAlg",              source = "bioc"),
  list(pkg = "BiocNeighbors",          source = "bioc"),
  list(pkg = "Azimuth",                source = "bioc"),
  list(pkg = "celldex",                source = "bioc"),
  list(pkg = "SingleR",                source = "bioc"),
  list(pkg = "scater",                 source = "bioc"),
  list(pkg = "miloR",                  source = "bioc"),
  list(pkg = "DESeq2",                 source = "bioc"),
  list(pkg = "clusterProfiler",        source = "bioc"),
  list(pkg = "org.Hs.eg.db",          source = "bioc"),
  list(pkg = "enrichplot",             source = "bioc"),
  list(pkg = "escape",                 source = "bioc"),
  list(pkg = "speckle",                source = "bioc"),
  list(pkg = "ComplexHeatmap",         source = "bioc"),
  list(pkg = "OmnipathR",              source = "bioc"),

  # Cross-compartment cell-cell communication (Phase 1d, Figure 4)
  list(pkg = "liana",        source = "github/saezlab/liana"),
  list(pkg = "nichenetr",    source = "github/saeyslab/nichenetr"),
  list(pkg = "circlize",     source = "cran"),

  # Plotting helpers used by panel D / loadings facets
  list(pkg = "ggh4x",        source = "cran"),
  list(pkg = "ggtext",       source = "cran"),
  list(pkg = "ggrepel",      source = "cran"),
  list(pkg = "scico",        source = "cran"),
  list(pkg = "matrixStats",  source = "cran"),

  # GitHub / specialty
  list(pkg = "scRepertoire", source = "github/ncborcherding/scRepertoire"),
  list(pkg = "alakazam",     source = "cran"),
  list(pkg = "dowser",       source = "cran"),
  list(pkg = "scoper",       source = "cran"),
  list(pkg = "shazam",       source = "cran"),
  list(pkg = "ggtree",       source = "bioc"),
  list(pkg = "scplotter",    source = "github/pwwang/scplotter"),

  # TCR motif/distance/generation analysis (slots 64–66, basilisk-managed Python)
  list(pkg = "immGLIPH",     source = "github/BorchLab/immGLIPH"),
  list(pkg = "immLynx",      source = "github/BorchLab/immLynx"),
  list(pkg = "ggseqlogo",    source = "cran"),
  list(pkg = "ggridges",     source = "cran"),
  list(pkg = "lmerTest",     source = "cran"),
  list(pkg = "ggraph",       source = "cran"),
  list(pkg = "tidygraph",    source = "cran"),

  # Advanced TCR analyses (slots 67–71, R/77 viz)
  list(pkg = "stringdist",   source = "cran"),
  list(pkg = "UCell",        source = "bioc")
)

# ---- Install missing packages ------------------------------------------------
.install_if_missing <- function(pkg, source) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))

  message(sprintf("Installing missing package: %s (from %s)...", pkg, source))

  if (source == "cran") {
    install.packages(pkg, quiet = TRUE)
  } else if (source == "bioc") {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", quiet = TRUE)
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else if (startsWith(source, "github/")) {
    repo <- sub("^github/", "", source)
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", quiet = TRUE)
    remotes::install_github(repo, upgrade = "never", quiet = TRUE)
  } else {
    warning(sprintf("Unknown source '%s' for package '%s'. Skipping.", source, pkg))
  }

  if (!requireNamespace(pkg, quietly = TRUE))
    warning(sprintf("Failed to install package: %s", pkg))
}

# ---- Main entry points -------------------------------------------------------

#' Install any missing packages and load all pipeline dependencies.
load_all_packages <- function() {
  # Install missing
  for (entry in .PACKAGE_MANIFEST) {
    .install_if_missing(entry$pkg, entry$source)
  }

  # Load all
  suppressPackageStartupMessages({
    for (entry in .PACKAGE_MANIFEST) {
      if (requireNamespace(entry$pkg, quietly = TRUE)) {
        library(entry$pkg, character.only = TRUE)
      }
    }
  })

  invisible(TRUE)
}

