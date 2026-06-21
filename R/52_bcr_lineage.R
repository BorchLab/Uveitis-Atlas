# R/52_bcr_lineage.R
# BCR lineage tree analysis using dowser
suppressPackageStartupMessages({
  library(Seurat)
  library(dowser)
  library(alakazam)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggtree)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

run_bcr_lineage <- function(cfg) {
  if (!isTRUE(cfg$bcr_lineage$enable)) {
    log_message("BCR lineage analysis disabled. Skipping.")
    return(invisible(TRUE))
  }

  log_message("Starting BCR lineage tree analysis...")

  ensure_dir(cfg$paths$results_tables)

  min_size   <- cfg$bcr_lineage$min_clone_size %||% 5
  max_trees  <- cfg$bcr_lineage$max_trees %||% 20
  tree_dir   <- viz_subdir(cfg$paths, "bcr_lineage")
  ensure_dir(tree_dir)

  # Read IMGT references
  imgt_dir <- cfg$paths$imgt_dir %||% "~/share/germlines/imgt/human/vdj"
  references <- tryCatch(dowser::readIMGT(dir = imgt_dir), error = function(e) {
    log_message("Could not read IMGT references: ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(references)) {
    log_message("IMGT references unavailable. Skipping BCR lineage.")
    return(invisible(TRUE))
  }

  # Read processed AIRR tables (clone-assigned + germline-reconstructed from step 20)
  airr_dir <- file.path(cfg$paths$results_tables, "bcr_airr")
  airr_files <- list.files(airr_dir, pattern = "_airr\\.tsv$", full.names = TRUE)

  if (length(airr_files) == 0) {
    log_message("No processed AIRR files found in ", airr_dir,
                ". Run annotate_AIRR (step 20) first. Skipping lineage.")
    return(invisible(TRUE))
  }

  bcr_db <- lapply(airr_files, read.delim, stringsAsFactors = FALSE) %>% bind_rows()
  log_message("Read ", nrow(bcr_db), " contigs from ", length(airr_files), " AIRR files")

  # Filter to productive heavy chains
  bcr_heavy <- bcr_db %>%
    filter(productive == TRUE, locus == "IGH")

  # Identify clones large enough for lineage tree building
  bcr_meta <- bcr_heavy %>%
    filter(!is.na(clone_id)) %>%
    group_by(subject_id, clone_id) %>%
    summarise(clone_size = n(), .groups = "drop") %>%
    filter(clone_size >= min_size) %>%
    arrange(desc(clone_size)) %>%
    slice_head(n = max_trees)

  log_message("Found ", nrow(bcr_meta), " clones with >= ", min_size, " cells")

  if (nrow(bcr_meta) == 0) {
    log_message("No clones meet minimum size. Skipping lineage trees.")
    return(invisible(TRUE))
  }

  lineage_summary <- list()

  for (i in seq_len(nrow(bcr_meta))) {
    cur_clone <- bcr_meta$clone_id[i]
    subject   <- bcr_meta$subject_id[i]
    n_cells   <- bcr_meta$clone_size[i]

    log_message("  Building tree for clone ", cur_clone, " (n=", n_cells, ")")

    tryCatch({
      clone_seqs <- bcr_heavy %>%
        filter(subject_id == subject, !is.na(clone_id)) %>%
        filter(as.character(clone_id) == as.character(cur_clone))

      if (nrow(clone_seqs) < 3) {
        log_message("    Too few sequences (", nrow(clone_seqs), "). Skipping.")
        next
      }

      # Reconstruct germlines if not present
      if (!"germline_alignment_d_mask" %in% colnames(clone_seqs)) {
        clone_seqs <- tryCatch(
          createGermlines(clone_seqs, references = references),
          error = function(e) { log_message("    Germline reconstruction failed."); NULL }
        )
        if (is.null(clone_seqs)) next
      }

      # Build phylogenetic tree via maximum parsimony
      clone_seqs$subject_id <- subject
      trees <- tryCatch(
        getTrees(clone_seqs, build = "pml", nproc = 1),
        error = function(e) NULL
      )

      if (!is.null(trees) && length(trees) > 0) {
        tree_file <- file.path(tree_dir, paste0("tree_", gsub("[^A-Za-z0-9_]", "_", cur_clone), ".pdf"))
        tryCatch({
          p <- plotTrees(trees)[[1]]
          ggsave(tree_file, p, width = 8, height = 6)
          log_message("    Saved: ", basename(tree_file))
        }, error = function(e) {
          log_message("    Tree plot failed: ", conditionMessage(e))
        })
      }

      lineage_summary[[i]] <- data.frame(
        clone_id = cur_clone,
        subject = subject,
        n_cells = n_cells,
        n_sequences = nrow(clone_seqs),
        tree_built = !is.null(trees),
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      log_message("  Tree building failed for ", cur_clone, ": ", conditionMessage(e))
    })
  }

  if (length(lineage_summary) > 0) {
    summary_df <- bind_rows(lineage_summary)
    write.csv(summary_df,
              file.path(cfg$paths$results_tables, "bcr_lineage_summary.csv"),
              row.names = FALSE)
    log_message("Saved: bcr_lineage_summary.csv")
  }

  log_message("BCR lineage analysis complete.")
  invisible(TRUE)
}
