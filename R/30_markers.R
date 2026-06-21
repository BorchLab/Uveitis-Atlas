# R/30_markers.R
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(DESeq2)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- FindAllMarkers (cluster markers) ---
run_markers <- function(cfg, target = c("all", "eye", "myeloid", "bcell", "tcell")) {
  target <- match.arg(target)
  paths  <- get_target_paths(cfg, target)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) return(invisible(TRUE))
  SeuratObj <- readRDS(obj_path)

  # --- All Markers for Cluster Differentiation ---
  DefaultAssay(SeuratObj) <- "RNA"
  all <- FindAllMarkers(SeuratObj,
                        only.pos = FALSE,
                        logfc.threshold = cfg$markers$logfc,
                        min.pct = cfg$markers$min_pct)
  ensure_dir(paths$results_tables)
  out_path <- file.path(paths$results_tables, "FindAllMarkers_objclusters.csv")
  write.csv(all, out_path, row.names = FALSE)
  log_message("Saved cluster markers: ", out_path)

  # Compartment branch also writes a curated top-N-per-substate table for the
  # F3-F5 marker dotplots (panel B). Filters: padj < 0.05, pct.diff > 0.4,
  # top 5 positive markers per substate by avg_log2FC.
  if (target %in% c("myeloid", "bcell", "tcell")) {
    gene_drop_re <- "^(LINC|MT-|ENSG|RPS|RPL|TRBC|TRBV|IGHV|IGKV|IGLV)"
    all <- all[!grepl(gene_drop_re, all$gene),]
    top <- all |>
      dplyr::filter(p_val_adj < 0.05, (pct.1 - pct.2) > 0.4, avg_log2FC > 0) |>
      dplyr::group_by(cluster) |>
      dplyr::slice_max(avg_log2FC, n = 5) |>
      dplyr::ungroup()
    top_path <- file.path(paths$results_tables, paste0(target, "_top_markers.csv"))
    write.csv(top, top_path, row.names = FALSE)
    log_message("Saved compartment top markers: ", top_path)
  }

  invisible(TRUE)
}

# --- Pseudobulk N audit (T0.3) -----------------------------------------------
#
# Read-only report: for each DGE contrast applicable to `target`, count
# n_subjects and n_cells per (cluster, group), with cohort breakdown. Use
# this BEFORE compartment DGE to flag substates that lack the power for a
# pseudobulk test. Does not modify run_dge behavior; the gate is enforced
# downstream (T2.1).
#
# Outputs:
#   outputs/tables/<target_path>/pseudobulk_N_per_group.csv
#
# Columns: contrast, cluster, group, n_subjects, n_cells,
# cohort_breakdown ("US=5|Japan=3" style), underpowered (n < 4 in either
# arm), low_power_flag (4 <= n <= 6).
run_pseudobulk_N_audit <- function(cfg, target = c("all", "eye", "myeloid", "bcell", "tcell")) {
  target <- match.arg(target)
  paths  <- get_target_paths(cfg, target)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Pseudobulk N audit (", target, "): integrated object missing; skipping.")
    return(invisible(FALSE))
  }
  obj  <- readRDS(obj_path)
  meta <- obj[[]]

  contrasts <- cfg$dge$contrasts %||% list()
  contrasts <- Filter(function(c) {
    ct <- c$target %||% "all"
    ct == target || ct == "any"
  }, contrasts)
  if (length(contrasts) == 0L) {
    log_message("Pseudobulk N audit (", target, "): no contrasts apply; skipping.")
    return(invisible(FALSE))
  }

  cluster_col <- if ("knn.leiden.cluster" %in% colnames(meta))
                   "knn.leiden.cluster" else "seurat_clusters"
  sample_col  <- if ("Subject_Timepoint" %in% colnames(meta)) "Subject_Timepoint"
                 else if ("Subject" %in% colnames(meta))      "Subject"
                 else                                          "orig.ident"

  meta$pb_cluster <- as.character(meta[[cluster_col]])
  clusters <- sort(unique(meta$pb_cluster))
  rows <- list()

  for (contrast in contrasts) {
    gcol <- contrast$group_col
    g1   <- contrast$group1
    g2   <- contrast$group2
    if (!(gcol %in% colnames(meta))) next

    for (cl in c("global", clusters)) {
      sub <- if (cl == "global") meta else meta[meta$pb_cluster == cl, , drop = FALSE]
      sub <- sub[sub[[gcol]] %in% c(g1, g2), , drop = FALSE]
      if (nrow(sub) == 0L) next

      for (grp in c(g1, g2)) {
        grp_meta   <- sub[sub[[gcol]] == grp, , drop = FALSE]
        n_subjects <- length(unique(grp_meta[[sample_col]]))
        n_cells    <- nrow(grp_meta)
        cohort_brk <- if ("Cohort" %in% colnames(grp_meta) && nrow(grp_meta) > 0L) {
          tbl <- table(grp_meta$Cohort, useNA = "ifany")
          paste(names(tbl), as.integer(tbl), sep = "=", collapse = "|")
        } else NA_character_
        rows[[length(rows) + 1L]] <- data.frame(
          contrast         = contrast$name,
          cluster          = cl,
          group            = grp,
          n_subjects       = n_subjects,
          n_cells          = n_cells,
          cohort_breakdown = cohort_brk,
          underpowered     = n_subjects < 4L,
          low_power_flag   = n_subjects >= 4L & n_subjects <= 6L,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) {
    log_message("Pseudobulk N audit (", target, "): no rows produced.")
    return(invisible(FALSE))
  }
  out_df <- dplyr::bind_rows(rows)
  ensure_dir(paths$results_tables)
  out_path <- file.path(paths$results_tables, "pseudobulk_N_per_group.csv")
  write.csv(out_df, out_path, row.names = FALSE)
  log_message("Wrote ", out_path)

  underpowered_rows <- out_df[out_df$underpowered, , drop = FALSE]
  if (nrow(underpowered_rows) > 0L) {
    log_message(sprintf(
      "Pseudobulk N audit (%s): %d (contrast, cluster, group) tuples have <4 subjects.",
      target, nrow(underpowered_rows)))
  }
  invisible(out_df)
}

# --- Pseudobulk DESeq2 + Wilcoxon DGE per contrast ---
run_dge <- function(cfg, target = c("all", "eye", "myeloid", "bcell", "tcell")) {
  target <- match.arg(target)
  paths  <- get_target_paths(cfg, target)
  obj_path <- file.path(paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Integrated object not found, skipping DGE.")
    return(invisible(TRUE))
  }
  obj <- readRDS(obj_path)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)

  contrasts <- cfg$dge$contrasts
  if (is.null(contrasts) || length(contrasts) == 0) {
    log_message("No DGE contrasts defined in config. Skipping.")
    return(invisible(TRUE))
  }

  # Filter contrasts to those tagged for this target. Defaults to "all" when
  # the target field is missing (backwards-compat with old configs).
  # target == "any" is a wildcard that fires for every target.
  contrasts <- Filter(function(c) {
    ct <- c$target %||% "all"
    ct == target || ct == "any"
  }, contrasts)
  if (length(contrasts) == 0) {
    log_message("No DGE contrasts for target=", target, ". Skipping.")
    return(invisible(TRUE))
  }

  ensure_dir(paths$results_tables)
  cluster_col <- "knn.leiden.cluster"

  for (contrast in contrasts) {
    cname    <- contrast$name
    gcol     <- contrast$group_col
    g1       <- contrast$group1
    g2       <- contrast$group2

    log_message("DGE contrast: ", cname, " (", g1, " vs ", g2, " on ", gcol, ")")

    # Subset to cells in these two groups
    cells_keep <- colnames(obj)[obj[[gcol, drop = TRUE]] %in% c(g1, g2)]
    if (length(cells_keep) < 100) {
      log_message("  Too few cells (", length(cells_keep), "), skipping contrast.")
      next
    }
    obj_sub <- subset(obj, cells = cells_keep)

    # ---- Pseudobulk DESeq2 (primary) ----
    log_message("  Running pseudobulk DESeq2...")
    tryCatch({
      pb_results <- run_pseudobulk_deseq2(
        obj_sub, group_col = gcol, group1 = g1, group2 = g2,
        cluster_col = cluster_col, cfg = cfg, target = target
      )
      write.csv(pb_results,
                file.path(paths$results_tables, paste0("DGE_pseudobulk_", cname, ".csv")),
                row.names = FALSE)
      log_message("  Saved pseudobulk results: DGE_pseudobulk_", cname, ".csv")

      # T2.1: For compartment targets running a cross-subject contrast,
      # also write within-cohort sensitivity passes. Single-cohort data
      # produces one CSV; multi-cohort produces one per cohort.
      if (target %in% c("myeloid", "bcell", "tcell") &&
          "Cohort" %in% colnames(obj_sub[[]])) {
        cohorts <- unique(stats::na.omit(obj_sub[["Cohort", drop = TRUE]]))
        if (length(cohorts) > 1L) {
          for (ck in cohorts) {
            cohort_cells <- colnames(obj_sub)[obj_sub[["Cohort", drop = TRUE]] == ck]
            if (length(cohort_cells) < 50L) next
            obj_ck <- subset(obj_sub, cells = cohort_cells)
            sens <- tryCatch(
              run_pseudobulk_deseq2(obj_ck, group_col = gcol, group1 = g1, group2 = g2,
                                    cluster_col = cluster_col, cfg = cfg,
                                    target = target, force_simple = TRUE),
              error = function(e) NULL)
            if (!is.null(sens)) {
              sens$cohort_sensitivity <- ck
              ck_slug <- tolower(gsub("[^A-Za-z0-9]+", "_", ck))
              p <- file.path(paths$results_tables,
                             sprintf("DGE_pseudobulk_%s_%s_only.csv", cname, ck_slug))
              write.csv(sens, p, row.names = FALSE)
              log_message("  Saved within-cohort sensitivity (", ck, "): ", basename(p))
            }
          }
        }
      }
    }, error = function(e) {
      log_message("  Pseudobulk DESeq2 failed: ", conditionMessage(e))
    })

    # ---- Seurat FindMarkers (cell-level, DIAGNOSTIC ONLY) ----
    # Gated behind `cfg$dge$run_wilcox` (default FALSE). This test uses each
    # CELL as an independent observation, which inflates p-values by orders
    # of magnitude when you have ~n samples but ~10^4 cells
    # (pseudoreplication). The pseudobulk DESeq2 results above are the
    # valid test for sample-level contrasts. Enable this only for quick
    # exploratory checks — do not cite these results.
    if (isTRUE(cfg$dge$run_wilcox)) {
      log_message("  [WARN] Running cell-level Wilcoxon DGE ",
                  "(pseudoreplication; diagnostic only)...")
      tryCatch({
        wilcox_results <- run_wilcox_dge(
          obj_sub, group_col = gcol, group1 = g1, group2 = g2,
          cluster_col = cluster_col, cfg = cfg
        )
        write.csv(wilcox_results,
                  file.path(paths$results_tables,
                            paste0("DGE_wilcox_", cname, ".csv")),
                  row.names = FALSE)
        log_message("  Saved Wilcoxon results: DGE_wilcox_", cname, ".csv")
      }, error = function(e) {
        log_message("  Wilcoxon DGE failed: ", conditionMessage(e))
      })
    }
  }

  # --- Build DGE summary statistics across all contrasts ---
  dge_summary_list <- list()
  for (contrast in contrasts) {
    cname <- contrast$name
    pb_path <- file.path(paths$results_tables, paste0("DGE_pseudobulk_", cname, ".csv"))
    if (!file.exists(pb_path)) next

    dge_df <- read.csv(pb_path, check.names = FALSE)
    if (nrow(dge_df) == 0 || !"padj" %in% colnames(dge_df)) next

    summary <- dge_df %>%
      filter(!is.na(padj)) %>%
      group_by(cluster) %>%
      summarise(
        n_tested = n(),
        n_sig = sum(padj < 0.05),
        n_up = sum(padj < 0.05 & log2FoldChange > 0),
        n_down = sum(padj < 0.05 & log2FoldChange < 0),
        top_gene_up = {
          up_idx <- which(padj < 0.05 & log2FoldChange > 0)
          if (length(up_idx) > 0) gene[up_idx[which.max(log2FoldChange[up_idx])]] else NA_character_
        },
        top_gene_down = {
          down_idx <- which(padj < 0.05 & log2FoldChange < 0)
          if (length(down_idx) > 0) gene[down_idx[which.min(log2FoldChange[down_idx])]] else NA_character_
        },
        .groups = "drop"
      ) %>%
      mutate(contrast = cname)

    dge_summary_list[[cname]] <- summary
  }

  if (length(dge_summary_list) > 0) {
    dge_summary <- bind_rows(dge_summary_list)
    write.csv(dge_summary,
              file.path(paths$results_tables, "dge_summary_stats.csv"),
              row.names = FALSE)
    log_message("Saved: dge_summary_stats.csv")
  }

  invisible(TRUE)
}

# --- T2.1 cohort design helpers ---------------------------------------------
#
# For compartment-level cross-subject contrasts (Autoimmune vs Viral, etc.)
# the audit (R3) wants `~ Cohort * group` as a primary design with graceful
# fallback to `~ Cohort + group` and then `~ group`. The functions below
# implement that decision tree. For the current single-cohort Paley dataset
# the selector will always return NULL and the caller keeps the legacy
# `~ group` design, but the code is ready when multi-cohort data arrives.

# Does the underlying object carry a usable Cohort column with at least
# two non-NA levels?
.has_cohort_signal <- function(obj, coldata) {
  if (!("Cohort" %in% colnames(obj[[]]))) return(FALSE)
  cv <- obj[["Cohort", drop = TRUE]]
  length(unique(stats::na.omit(cv))) >= 2L
}

# Attach Cohort to the pseudobulk coldata by sample mapping. Returns the
# coldata unchanged when Cohort cannot be resolved unambiguously per sample.
.attach_cohort <- function(obj, coldata) {
  meta <- obj[[]]
  if (!"Cohort" %in% colnames(meta)) return(coldata)
  sample_col <- if ("Subject_Timepoint" %in% colnames(meta)) "Subject_Timepoint"
                else if ("orig.ident"   %in% colnames(meta)) "orig.ident"
                else "Subject"
  if (!(sample_col %in% colnames(meta))) return(coldata)
  per_sample <- tapply(as.character(meta$Cohort),
                       as.character(meta[[sample_col]]),
                       function(x) {
                         ux <- unique(stats::na.omit(x))
                         if (length(ux) == 1L) ux else NA_character_
                       })
  coldata$Cohort <- factor(unname(per_sample[as.character(coldata$sample)]))
  coldata
}

# Choose between ~ Cohort * group / ~ Cohort + group / NULL.
#   NULL  → caller should fall back to its legacy design (~ group or paired).
# Returns: list($design = formula, $model_used = c("interaction","additive"))
.select_cohort_design <- function(coldata) {
  if (!("Cohort" %in% colnames(coldata))) return(NULL)
  if (length(unique(stats::na.omit(coldata$Cohort))) < 2L) return(NULL)
  if (anyNA(coldata$Cohort)) return(NULL)

  # Cross-tab: every (Cohort, group) cell must be populated for interaction.
  ct <- table(coldata$Cohort, coldata$group)
  all_populated <- all(ct > 0L)

  if (all_populated) {
    mm_int <- try(model.matrix(~ Cohort * group, data = coldata), silent = TRUE)
    if (!inherits(mm_int, "try-error") &&
        qr(mm_int)$rank == ncol(mm_int) &&
        nrow(mm_int) > ncol(mm_int)) {
      return(list(design = ~ Cohort * group, model_used = "interaction"))
    }
  }
  mm_add <- try(model.matrix(~ Cohort + group, data = coldata), silent = TRUE)
  if (!inherits(mm_add, "try-error") &&
      qr(mm_add)$rank == ncol(mm_add) &&
      nrow(mm_add) > ncol(mm_add)) {
    return(list(design = ~ Cohort + group, model_used = "additive"))
  }
  NULL
}

# --- Pseudobulk DESeq2 helper ---
# Aggregates counts per (sample x group) via explicit indicator-matrix
# aggregation (no Seurat name parsing — avoids the substring-matching bug
# where sample IDs that are prefixes of each other corrupt coldata).
#
# Design: adds `+ Subject` blocking when the contrast varies WITHIN subject
# (e.g. Tissue_1 for paired Eye/Blood subjects). Falls back to `~ group`
# when Subject is aliased with group (phenotype contrasts where each subject
# is a single phenotype) or when the resulting design would be rank-deficient.
#
# T2.1: when `target %in% c("myeloid","bcell","tcell")` and the contrast is
# cross-subject, the .select_cohort_design() helper above can upgrade the
# design to `~ Cohort * group` or `~ Cohort + group` if the data supports
# it. `force_simple = TRUE` bypasses cohort upgrades — used for within-
# cohort sensitivity passes.
run_pseudobulk_deseq2 <- function(obj, group_col, group1, group2,
                                   cluster_col, cfg,
                                   target = "all",
                                   force_simple = FALSE) {
  suppressPackageStartupMessages({
    library(Matrix)
  })

  # Sample-level ID for aggregation (prefer Subject_Timepoint > orig.ident > Subject)
  if ("Subject_Timepoint" %in% colnames(obj[[]])) {
    sample_col <- "Subject_Timepoint"
  } else if ("orig.ident" %in% colnames(obj[[]])) {
    sample_col <- "orig.ident"
  } else {
    sample_col <- "Subject"
  }
  has_subject <- "Subject" %in% colnames(obj[[]])

  meta <- obj[[]]
  meta$pb_cluster <- meta[[cluster_col]]
  clusters <- sort(unique(meta$pb_cluster))
  cluster_list <- c("global", as.character(clusters))

  all_results <- list()

  # --- Aggregate counts per (sample, group) using a composite key --------
  aggregate_pb <- function(obj_cl, sample_col, group_col) {
    m <- obj_cl[[]]
    # Composite key — sentinel separator avoids collisions with values
    m$pb_key <- paste(as.character(m[[sample_col]]),
                      as.character(m[[group_col]]),
                      sep = "___KEY___")
    counts_mat <- GetAssayData(obj_cl, assay = "RNA", layer = "counts")
    keys <- m$pb_key
    uk   <- unique(keys)
    # cell x key indicator matrix, sum via matrix multiply
    ind <- Matrix::sparseMatrix(
      i = seq_along(keys),
      j = match(keys, uk),
      x = 1,
      dims = c(length(keys), length(uk))
    )
    agg <- counts_mat %*% ind
    colnames(agg) <- uk

    # Coldata by exact lookup — no string parsing
    cd <- m %>%
      dplyr::distinct(pb_key, !!sym(sample_col), !!sym(group_col),
                      !!!(if (has_subject) rlang::syms("Subject") else list())) %>%
      as.data.frame()
    rownames(cd) <- cd$pb_key
    cd <- cd[colnames(agg), , drop = FALSE]
    colnames(cd)[colnames(cd) == sample_col] <- "sample"
    colnames(cd)[colnames(cd) == group_col]  <- "group"
    list(agg = agg, coldata = cd)
  }

  for (cl in cluster_list) {
    cells <- if (cl == "global") colnames(obj) else colnames(obj)[meta$pb_cluster == cl]
    if (length(cells) < 20) next

    obj_cl <- subset(obj, cells = cells)
    pb <- tryCatch(aggregate_pb(obj_cl, sample_col, group_col),
                   error = function(e) { log_message("  pb agg failed (", cl, "): ", conditionMessage(e)); NULL })
    if (is.null(pb)) next
    agg <- pb$agg
    coldata <- pb$coldata

    # Keep only the two contrast levels
    coldata <- coldata[coldata$group %in% c(group1, group2) & !is.na(coldata$group), , drop = FALSE]
    agg     <- agg[, rownames(coldata), drop = FALSE]
    if (nrow(coldata) < 4) next
    if (length(unique(coldata$group)) < 2) next
    if (any(table(coldata$group) < 2)) next

    coldata$group <- factor(coldata$group, levels = c(group2, group1))

    # --- Decide on design: paired (+ Subject) when feasible -----------------
    use_paired <- FALSE
    if (has_subject && "Subject" %in% colnames(coldata)) {
      coldata$Subject <- factor(as.character(coldata$Subject))
      # Does the contrast vary WITHIN Subject?
      per_sbj_grp <- tapply(as.character(coldata$group),
                            as.character(coldata$Subject),
                            function(x) length(unique(x)))
      within_variation <- any(per_sbj_grp > 1, na.rm = TRUE)
      if (within_variation) {
        # Rank check: ~ Subject + group on current coldata
        mm_try <- try(model.matrix(~ Subject + group, data = coldata), silent = TRUE)
        if (!inherits(mm_try, "try-error") &&
            qr(mm_try)$rank == ncol(mm_try) &&
            nrow(mm_try) > ncol(mm_try)) {
          use_paired <- TRUE
        }
      }
    }

    # --- T2.1 cohort-aware design for compartment targets -------------------
    # For compartment targets (myeloid/bcell/tcell) running a cross-subject
    # contrast (use_paired = FALSE), check whether the data supports
    # ~ Cohort * group → ~ Cohort + group → ~ group as a graceful fallback.
    # For target %in% c("all", "eye"), or force_simple = TRUE (within-cohort
    # sensitivity passes), keep the legacy behavior.
    cohort_model_meta <- list(model_used = NA_character_,
                              interaction_p = NA_real_)
    use_cohort_design <- FALSE
    cohort_design_f   <- NULL

    if (!force_simple &&
        target %in% c("myeloid", "bcell", "tcell") &&
        !use_paired &&
        .has_cohort_signal(obj, coldata)) {
      cd <- .attach_cohort(obj, coldata)
      sel <- .select_cohort_design(cd)
      if (!is.null(sel)) {
        coldata           <- cd
        use_cohort_design <- TRUE
        cohort_design_f   <- sel$design
        cohort_model_meta$model_used <- sel$model_used
      }
    }

    design_f <- if (use_cohort_design) cohort_design_f
                else if (use_paired)   ~ Subject + group
                else                    ~ group

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(as.matrix(agg)),
      colData   = coldata,
      design    = design_f
    )
    keep <- rowSums(counts(dds) >= 10) >= 2
    dds  <- dds[keep, ]
    if (nrow(dds) < 10) next

    dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) {
      log_message("  DESeq failed (", cl, ", paired=", use_paired, "): ",
                  conditionMessage(e))
      NULL
    })
    if (is.null(dds)) next

    res <- DESeq2::results(dds, contrast = c("group", group1, group2))

    # When the interaction model fit, record the LRT p on the interaction
    # term so downstream filters can flag etiology effects that differ by
    # cohort. NA for additive / simple / paired designs.
    interaction_p <- NA_real_
    if (use_cohort_design && cohort_model_meta$model_used == "interaction") {
      interaction_p <- tryCatch({
        rn <- DESeq2::resultsNames(dds)
        ix <- grep(":", rn, value = TRUE)
        if (length(ix) > 0L) {
          r2 <- DESeq2::results(dds, name = ix[1])
          stats::median(r2$pvalue, na.rm = TRUE)   # summary across genes
        } else NA_real_
      }, error = function(e) NA_real_)
    }

    model_used <- if (use_cohort_design) cohort_model_meta$model_used
                  else if (use_paired)   "paired_subject"
                  else                   "simple_group"

    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      mutate(cluster      = cl,
             contrast     = paste0(group1, "_vs_", group2),
             paired       = use_paired,
             model_used   = model_used,
             interaction_p_median = interaction_p) %>%
      filter(!is.na(padj)) %>%
      arrange(padj)

    all_results[[length(all_results) + 1]] <- res_df
  }

  if (length(all_results) == 0) {
    return(data.frame(gene = character(), cluster = character()))
  }
  bind_rows(all_results)
}

# --- Wilcoxon DGE helper ---
run_wilcox_dge <- function(obj, group_col, group1, group2,
                            cluster_col, cfg) {
  meta <- obj[[]]
  clusters <- sort(unique(meta[[cluster_col]]))
  all_results <- list()

  # Global DGE
  tryCatch({
    Idents(obj) <- group_col
    global_res <- FindMarkers(obj,
                              ident.1 = group1,
                              ident.2 = group2,
                              logfc.threshold = cfg$dge$logfc %||% 0.25,
                              min.pct = cfg$dge$min_pct %||% 0.1,
                              verbose = FALSE)
    if (nrow(global_res) > 0) {
      global_res <- global_res %>%
        rownames_to_column("gene") %>%
        mutate(cluster = "global", contrast = paste0(group1, "_vs_", group2))
      all_results[[length(all_results) + 1]] <- global_res
    }
  }, error = function(e) {
    log_message("  Global Wilcoxon failed: ", conditionMessage(e))
  })

  # Per-cluster DGE
  for (cl in clusters) {
    tryCatch({
      cells <- colnames(obj)[meta[[cluster_col]] == cl]
      if (length(cells) < 20) next

      obj_cl <- subset(obj, cells = cells)
      g1_n <- sum(obj_cl[[group_col, drop = TRUE]] == group1)
      g2_n <- sum(obj_cl[[group_col, drop = TRUE]] == group2)
      if (g1_n < 3 || g2_n < 3) next

      Idents(obj_cl) <- group_col
      cl_res <- FindMarkers(obj_cl,
                            ident.1 = group1,
                            ident.2 = group2,
                            logfc.threshold = cfg$dge$logfc %||% 0.25,
                            min.pct = cfg$dge$min_pct %||% 0.1,
                            verbose = FALSE)

      if (nrow(cl_res) > 0) {
        cl_res <- cl_res %>%
          rownames_to_column("gene") %>%
          mutate(cluster = as.character(cl),
                 contrast = paste0(group1, "_vs_", group2))
        all_results[[length(all_results) + 1]] <- cl_res
      }
    }, error = function(e) {
      log_message("  Cluster ", cl, " Wilcoxon failed: ", conditionMessage(e))
    })
  }

  if (length(all_results) == 0) {
    return(data.frame(gene = character(), cluster = character()))
  }

  bind_rows(all_results)
}
