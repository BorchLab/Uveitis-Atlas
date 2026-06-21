# R/66_immlynx_olga_sonnia.R
# Generation probability (OLGA) and post-selection inference (SoNNia) on
# intraocular TRB. Mechanistic hypothesis: viral clones are antigen-driven
# and should have lower median Pgen than NIU (which trends closer to
# unselected background). SoNNia adds a selection factor Q per V/J/length
# feature, tested for differential abundance Viral vs NIU.
#
# Cell-number scaling: per-clone Pgen is intrinsically size-invariant.
# We aggregate Pgen to one value per subject (median of log10 Pgen across
# the clones present in that subject) and fit:
#   lm(med_log10_Pgen ~ Phenotype_2 + log10(n_cells))
# which absorbs the cell-count confound.

`%||%` <- function(x, y) if (is.null(x)) y else x

run_immlynx_olga <- function(cfg) {
  if (!isTRUE(cfg$olga$enable) || !isTRUE(cfg$steps$olga)) {
    log_message("OLGA / SoNNia disabled. Skipping.")
    return(invisible(TRUE))
  }
  if (!requireNamespace("immLynx", quietly = TRUE)) {
    log_message("immLynx not installed; skipping OLGA step.")
    return(invisible(FALSE))
  }

  log_message("Starting OLGA + SoNNia on intraocular TRB...")

  obj_path <- file.path(cfg$paths$results_objects, "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("Integrated object not found. Skipping OLGA.")
    return(invisible(TRUE))
  }
  obj <- readRDS(obj_path)

  out_tables <- file.path(cfg$paths$results_tables, "repertoire")
  out_obj    <- cfg$paths$results_objects
  ensure_dir(out_tables); ensure_dir(out_obj)

  # prepare_intraocular_trb() + trb_frame_to_sce() come from R/64_immgliph.R.
  trb <- prepare_intraocular_trb(obj, cfg)
  log_message("  Intraocular TRB rows: ", nrow(trb))
  if (nrow(trb) < 50) {
    log_message("  Too few TRB sequences. Aborting OLGA.")
    return(invisible(FALSE))
  }

  sce <- trb_frame_to_sce(trb)

  # ---- OLGA generation probability ----------------------------------------
  model_arg <- cfg$olga$model %||% "humanTRB"   # immLynx model naming
  log_message("  Running runOLGA (model=", model_arg, ")...")
  sce <- tryCatch(
    immLynx::runOLGA(sce, chains = "TRB", model = model_arg,
                     use_vj_genes = isTRUE(cfg$olga$use_vj_genes),
                     return_object = TRUE,
                     column_name = "olga_pgen"),
    error = function(e) {
      log_message("  runOLGA failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(sce)) return(invisible(FALSE))

  cd <- SummarizedExperiment::colData(sce)
  # immLynx names the columns "<column_name>_<chain>" and also adds a
  # precomputed log10 variant ("..._log10_<chain>"). Match either layout.
  pgen_col <- intersect(c("olga_pgen_TRB", "olga_pgen"), colnames(cd))[1]
  log_col  <- intersect(c("olga_pgen_log10_TRB", "olga_pgen_log10"),
                        colnames(cd))[1]
  if (is.na(pgen_col)) {
    pgen_col <- grep("^olga_pgen", colnames(cd), value = TRUE)[1]
    pgen_col <- pgen_col[!grepl("log10", pgen_col)][1]
  }
  if (is.na(log_col))
    log_col <- grep("^olga_pgen.*log10", colnames(cd), value = TRUE)[1]

  if (is.na(pgen_col)) {
    log_message("  No OLGA Pgen column found in colData. Available: ",
                paste(grep("olga|pgen", colnames(cd), value = TRUE,
                           ignore.case = TRUE), collapse = ", "))
    return(invisible(FALSE))
  }
  log_message("  Reading Pgen from colData column '", pgen_col, "'",
              if (!is.na(log_col)) paste0(" (log10 from '", log_col, "')")
              else "", ".")

  pgen_vec <- as.numeric(cd[[pgen_col]])
  log10_vec <- if (!is.na(log_col)) as.numeric(cd[[log_col]])
               else log10(pmax(pgen_vec, .Machine$double.xmin))

  per_clone <- data.frame(
    barcode             = colnames(sce),
    Subject             = trb$Subject,
    Etiology            = trb$Etiology,
    Phenotype_2         = trb$Phenotype_2,
    CDR3b               = trb$CDR3b,
    TRBV                = trb$TRBV,
    TRBJ                = trb$TRBJ,
    HLA_B27_pathogenic  = trb$HLA_B27_pathogenic,
    Pgen                = pgen_vec,
    log10_Pgen          = log10_vec,
    stringsAsFactors    = FALSE
  )
  write.csv(per_clone,
            file.path(out_tables, "olga_pgen_per_clone.csv"),
            row.names = FALSE)
  log_message("  Saved: olga_pgen_per_clone.csv")

  per_subject <- per_clone |>
    dplyr::filter(is.finite(log10_Pgen)) |>
    dplyr::group_by(Subject, Etiology, Phenotype_2) |>
    dplyr::summarise(n_cells         = dplyr::n(),
                     med_log10_Pgen  = stats::median(log10_Pgen),
                     mean_log10_Pgen = mean(log10_Pgen),
                     min_log10_Pgen  = min(log10_Pgen),
                     .groups = "drop")
  write.csv(per_subject,
            file.path(out_tables, "olga_per_subject.csv"),
            row.names = FALSE)
  log_message("  Saved: olga_per_subject.csv")

  glm_res <- tryCatch({
    fit <- stats::lm(med_log10_Pgen ~ Phenotype_2 + log10(n_cells),
                     data = per_subject)
    cs  <- summary(fit)$coefficients
    data.frame(metric   = "med_log10_Pgen",
               term     = rownames(cs),
               estimate = cs[, "Estimate"],
               std_err  = cs[, "Std. Error"],
               t        = cs[, "t value"],
               p        = cs[, "Pr(>|t|)"])
  }, error = function(e) {
    log_message("  GLM failed: ", conditionMessage(e)); NULL
  })
  if (!is.null(glm_res)) {
    write.csv(glm_res,
              file.path(out_tables, "olga_glm_results.csv"),
              row.names = FALSE)
    log_message("  Saved: olga_glm_results.csv")
  }

  # ---- SoNNia selection inference -----------------------------------------
  # SoNNia requires a background CSV produced by generateOLGA. We write to
  # a tempfile alongside outputs so the run is reproducible from the table.
  sonnia_df <- NULL
  if (isTRUE(cfg$olga$sonnia)) {
    n_gen <- cfg$olga$n_generate_background %||% 100000
    bg_path <- file.path(out_tables, "olga_background_TRB.csv")
    sonia_dir <- file.path(out_obj, "sonia_output")
    n_epochs <- cfg$olga$sonnia_epochs %||% 50

    log_message("  Generating OLGA background (n=", n_gen,
                ") for SoNNia...")
    bg_ok <- tryCatch({
      bg <- immLynx::generateOLGA(n = n_gen, model = model_arg)
      utils::write.csv(bg, bg_path, row.names = FALSE)
      TRUE
    }, error = function(e) {
      log_message("  generateOLGA failed: ", conditionMessage(e)); FALSE
    })

    if (bg_ok) {
      log_message("  Running runSoNNia (n_epochs=", n_epochs, ")...")
      sonnia_sce <- tryCatch(
        immLynx::runSoNNia(sce,
                           chains          = "TRB",
                           background_file = bg_path,
                           save_folder     = sonia_dir,
                           n_epochs        = n_epochs,
                           return_object   = TRUE),
        error = function(e) {
          log_message("  runSoNNia failed: ", conditionMessage(e)); NULL
        }
      )
      if (!is.null(sonnia_sce) &&
          inherits(sonnia_sce, "SingleCellExperiment")) {
        cd2 <- SummarizedExperiment::colData(sonnia_sce)
        sel_col <- intersect(c("sonia_Q", "sonia_q", "selection_factor",
                               "Q"), colnames(cd2))[1]
        if (is.na(sel_col))
          sel_col <- grep("^son(ia|nia)_", colnames(cd2), value = TRUE)[1]
        if (!is.na(sel_col)) {
          sonnia_df <- data.frame(
            barcode     = colnames(sonnia_sce),
            Subject     = trb$Subject,
            Etiology    = trb$Etiology,
            Phenotype_2 = trb$Phenotype_2,
            CDR3b       = trb$CDR3b,
            TRBV        = trb$TRBV,
            TRBJ        = trb$TRBJ,
            Q           = as.numeric(cd2[[sel_col]]),
            stringsAsFactors = FALSE
          )
          write.csv(sonnia_df,
                    file.path(out_tables, "sonnia_selection_per_clone.csv"),
                    row.names = FALSE)
          log_message("  Saved: sonnia_selection_per_clone.csv")

          # Per-feature contrast: Phenotype_2 effect on log2 Q at the
          # V family x CDR3-length-bin level.
          feat <- sonnia_df |>
            dplyr::filter(is.finite(Q), Q > 0) |>
            dplyr::mutate(
              V_family = stringr::str_extract(TRBV, "TRBV[0-9]+"),
              len_bin  = cut(nchar(CDR3b),
                             breaks = c(0, 11, 13, 15, 17, 30),
                             labels = c("<=11", "12-13", "14-15",
                                        "16-17", ">=18"))) |>
            dplyr::group_by(V_family, len_bin) |>
            dplyr::filter(dplyr::n() >= 10) |>
            dplyr::summarise(
              n_clones = dplyr::n(),
              log2_Q   = log2(stats::median(Q)),
              p = tryCatch(
                stats::wilcox.test(log(Q) ~ Phenotype_2)$p.value,
                error = function(e) NA_real_),
              .groups = "drop"
            ) |>
            dplyr::mutate(FDR = stats::p.adjust(p, method = "BH"),
                          feature = paste(V_family, len_bin, sep = ":"))
          write.csv(feat,
                    file.path(out_tables, "sonnia_selection_features.csv"),
                    row.names = FALSE)
          log_message("  Saved: sonnia_selection_features.csv (",
                      nrow(feat), " features)")
        } else {
          log_message("  No SoNNia score column found in colData.")
        }
      }
    }
  }

  saveRDS(list(per_clone   = per_clone,
               per_subject = per_subject,
               glm         = glm_res,
               sonnia      = sonnia_df),
          file.path(out_obj, "ImmLynxOlgaResults.rds"))
  log_message("  Saved: ImmLynxOlgaResults.rds")

  log_message("OLGA + SoNNia analysis complete.")
  invisible(TRUE)
}
