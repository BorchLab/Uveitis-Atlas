# R/41_composition.R
# Differential composition testing using propeller (speckle)
suppressPackageStartupMessages({
  library(Seurat)
  library(speckle)
  library(dplyr)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

run_composition_testing <- function(cfg, target = c("all", "eye", "myeloid", "bcell", "tcell")) {
  target <- match.arg(target)
  paths  <- get_target_paths(cfg, target)

  if (!isTRUE(cfg$composition$enable)) {
    log_message("Composition testing disabled. Skipping.")
    return(invisible(TRUE))
  }

  log_message("Starting composition testing (propeller, target=", target, ")...")

  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Integrated object not found. Skipping composition testing.")
    return(invisible(TRUE))
  }

  obj <- readRDS(obj_path)
  ensure_dir(paths$results_tables)

  # Compartment objects carry the broad celltype as celltype_broad_eye
  # (preserved from the eye sub-atlas). The substate column is
  # knn.leiden.cluster. Require at least one usable grouping column.
  has_broad     <- "celltype_broad"     %in% colnames(obj[[]])
  has_broad_eye <- "celltype_broad_eye" %in% colnames(obj[[]])
  has_substate  <- "knn.leiden.cluster" %in% colnames(obj[[]])
  if (!has_broad && !has_broad_eye && !has_substate) {
    log_message("No celltype_broad / celltype_broad_eye / knn.leiden.cluster column. Skipping composition testing.")
    return(invisible(TRUE))
  }

  # Need a sample-level identifier
  if ("Subject" %in% colnames(obj[[]])) {
    sample_col <- "Subject"
  } else {
    sample_col <- "orig.ident"
  }

  meta <- obj[[]]

  # Define contrasts (target-dependent: cross-tissue on full obj, phenotype on eye obj)
  comp_contrasts <- if (target == "all") {
    list(
      list(name = "Eye_vs_Blood", col = "Tissue_1", groups = c("Eye", "Blood"))
    )
  } else {
    list(
      list(name = "Autoimmune_vs_Viral", col = "Phenotype_2", groups = c("NIU",  "Viral")),
      list(name = "Gran_vs_NonGran",     col = "Phenotype_2", groups = c("Gran", "Non_Gran"))
    )
  }

  # Run propeller for both celltype_broad AND merged.celltype.cluster when
  # available, with paired Subject blocking for contrasts that vary within
  # subject (Tissue_1). Propeller accepts a transform + design matrix via
  # propeller.ttest / propeller.anova, but the simpler `propeller()` wrapper
  # doesn't expose a formula — so for paired tests we use getTransformedProps
  # + limma::lmFit manually. Fall back to the unpaired propeller() call when
  # pairing is not applicable or rank-deficient.
  has_subject <- "Subject" %in% colnames(meta)

  run_one_contrast <- function(cname, gcol, groups, group_col_values) {
    log_message("  Composition test: ", cname, " on ", group_col_values)
    cells_keep <- which(meta[[gcol]] %in% groups)
    if (length(cells_keep) < 100) {
      log_message("  Too few cells; skipping.")
      return(NULL)
    }
    sub <- meta[cells_keep, , drop = FALSE]

    # Decide whether to block by Subject
    use_paired <- FALSE
    if (has_subject) {
      per_sbj <- tapply(as.character(sub[[gcol]]),
                        as.character(sub$Subject),
                        function(x) length(unique(x)))
      within_variation <- any(per_sbj > 1, na.rm = TRUE)
      if (within_variation) {
        # Build sample-level design matrix and check rank
        samp_df <- unique(sub[, c(sample_col, gcol, "Subject"), drop = FALSE])
        samp_df$Subject <- factor(as.character(samp_df$Subject))
        samp_df[[gcol]] <- factor(as.character(samp_df[[gcol]]))
        mm_try <- try(model.matrix(as.formula(paste("~ Subject +", gcol)),
                                   data = samp_df), silent = TRUE)
        if (!inherits(mm_try, "try-error") &&
            qr(mm_try)$rank == ncol(mm_try) &&
            nrow(mm_try) > ncol(mm_try)) {
          use_paired <- TRUE
        } else {
          log_message("  Paired design rank-deficient; using unpaired propeller.")
        }
      }
    }

    clust_vec <- sub[[group_col_values]]
    samp_vec  <- sub[[sample_col]]
    grp_vec   <- sub[[gcol]]

    if (!use_paired) {
      return(propeller(clusters = clust_vec,
                       sample   = samp_vec,
                       group    = grp_vec))
    }

    # Paired path: propeller manual pipeline
    # 1. Transform proportions (arcsin-sqrt by default)
    props <- speckle::getTransformedProps(
      clusters = clust_vec,
      sample   = samp_vec,
      transform = "asin"
    )
    # 2. Build sample-level design matching props$TransformedProps columns
    sample_ids <- colnames(props$TransformedProps)
    lookup <- sub %>%
      dplyr::distinct(!!sym(sample_col), !!sym(gcol), Subject) %>%
      as.data.frame()
    rownames(lookup) <- as.character(lookup[[sample_col]])
    lookup <- lookup[sample_ids, , drop = FALSE]
    lookup$Subject    <- factor(as.character(lookup$Subject))
    lookup[[gcol]]    <- factor(as.character(lookup[[gcol]]),
                                levels = groups)  # group1 = reference
    design_mat <- model.matrix(as.formula(paste("~ Subject +", gcol)),
                               data = lookup)
    # 3. limma linear model; last coefficient is the tested group effect
    fit <- limma::lmFit(props$TransformedProps, design = design_mat)
    fit <- limma::eBayes(fit, robust = TRUE)
    coef_name <- tail(colnames(design_mat), 1)
    tab <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
    # Match propeller()'s output structure as closely as possible
    tab$BaselineProp <- rowMeans(props$Proportions, na.rm = TRUE)
    tab <- tab[, c("BaselineProp", "logFC", "t", "P.Value", "adj.P.Val"),
               drop = FALSE]
    colnames(tab) <- c("BaselineProp.Freq", "PropMean.logFC",
                       "Tstatistic", "P.Value", "FDR")
    tab
  }

  # Eye and compartment objects skip 12_merge_clusters, so use Leiden directly.
  # Compartments inherit broad lineage as celltype_broad_eye.
  groupings <- if (target == "all") {
    c("celltype_broad", "merged.celltype.cluster")
  } else if (target == "eye") {
    c("celltype_broad", "knn.leiden.cluster")
  } else {
    c("celltype_broad_eye", "knn.leiden.cluster")
  }

  for (contrast in comp_contrasts) {
    cname <- contrast$name
    gcol  <- contrast$col
    groups <- contrast$groups

    for (gc in groupings) {
      if (!gc %in% colnames(meta)) next
      tryCatch({
        res <- run_one_contrast(cname, gcol, groups, gc)
        if (is.null(res)) next
        res$contrast <- cname
        res$grouping <- gc
        out_path <- file.path(paths$results_tables,
                              paste0("composition_test_", cname, "_", gc, ".csv"))
        write.csv(res, out_path, row.names = TRUE)
        log_message("  Saved: ", basename(out_path))
      }, error = function(e) {
        log_message("  Composition test failed for ", cname, "/", gc, ": ",
                    conditionMessage(e))
      })

      # --- Secondary: per-cluster Fisher's exact + binomial GLMM -------------
      tryCatch({
        cells_keep <- which(meta[[gcol]] %in% groups)
        if (length(cells_keep) < 100) {
          log_message("  Fisher+GLMM: too few cells for ", cname, "/", gc,
                      "; skipping.")
          next
        }
        sub <- meta[cells_keep, , drop = FALSE]
        fg <- run_fisher_glmm_per_cluster(
          meta         = sub,
          contrast_col = gcol,
          groups       = groups,
          cluster_col  = gc,
          sample_col   = if ("Subject" %in% colnames(sub)) "Subject" else sample_col
        )
        fg$contrast <- cname
        fg$grouping <- gc
        fg_path <- file.path(paths$results_tables,
                             paste0("composition_fisher_glmm_", cname,
                                    "_", gc, ".csv"))
        write.csv(fg, fg_path, row.names = FALSE)
        log_message("  Saved: ", basename(fg_path))
      }, error = function(e) {
        log_message("  Fisher+GLMM failed for ", cname, "/", gc, ": ",
                    conditionMessage(e))
      })
    }
  }

  log_message("Composition testing complete.")
  invisible(TRUE)
}
