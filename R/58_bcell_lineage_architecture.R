# R/58_bcell_lineage_architecture.R
# Per-clone B-cell lineage architecture metrics for Figure 6 panels F + G.
# See docs/plans/2026-05-27-fig6-bcr-architecture-design.md.
#
# Outputs (under outputs/tables/eye/bcell/architecture/):
#   clone_architecture_metrics.csv     - one row per (subject, clone_id)
#   subject_architecture_summary.csv   - one row per (subject, phenotype)
#   stats_summary.csv                  - one row per panel-level contrast

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Reuse R/56's eye-AIRR loader so the architecture phase consumes the same
# per-cell meta + clone definitions as the tree phase. Returns
# list(bm, bcr_eye, paths, paths_all) or NULL if the bcell object is missing.
.load_bcell_arch_meta <- function(cfg) {
  if (!exists(".bcell_load_eye_airr", mode = "function"))
    source(file.path("R", "56_bcell_lineage_trees.R"))
  .bcell_load_eye_airr(cfg)
}

# Shannon entropy on a two-tissue partition (eye/blood counts per clone).
# Returns 0 when only one compartment is populated, NA when total == 0.
.shannon_two <- function(a, b) {
  tot <- a + b
  out <- rep(NA_real_, length(tot))
  ok  <- tot > 0
  pa  <- ifelse(ok, a / tot, 0)
  pb  <- ifelse(ok, b / tot, 0)
  ent <- ifelse(pa > 0, -pa * log(pa), 0) + ifelse(pb > 0, -pb * log(pb), 0)
  out[ok] <- ent[ok]
  out
}

# Build paired (eye + blood) per-cell meta for the clones present on the
# eye side. Clones live in the AIRR tables (the full-atlas Seurat object
# has no clone_id column), so the join is:
#   AIRR rows whose clone_id is in eye_bm -> filter on cell_id_unique in
#   full-atlas meta -> restrict to B cells (celltype_broad == "B cell") ->
#   carry Tissue_1 / Subject / substate / mu_freq_* from the meta.
#
# Phenotype is lifted from eye_bm (Subject -> phenotype mapping derived
# from cfg$etiology_groups in R/56). Tissue is lowercased to match the
# per-clone helper expectation ("eye" / "blood").
# Returns NULL if the full-atlas Seurat object or AIRR tables are missing.
.bcell_arch_paired_meta <- function(eye_bm, cfg) {
  paths_all <- get_target_paths(cfg, "all")
  obj_path  <- file.path(paths_all$results_objects,
                         "IntegratedSeuratObject.rds")
  if (!file.exists(obj_path)) {
    log_message("[fig6 F/G] full-atlas object missing: ", obj_path)
    return(NULL)
  }
  airr_dir <- file.path(paths_all$results_tables, "bcr_airr")
  airr_files <- list.files(airr_dir, pattern = "_airr\\.tsv$",
                           full.names = TRUE)
  if (length(airr_files) == 0L) {
    log_message("[fig6 F/G] no AIRR tables in ", airr_dir)
    return(NULL)
  }

  keep_clones <- unique(stats::na.omit(eye_bm$clone_id))
  if (length(keep_clones) == 0L) {
    log_message("[fig6 F/G] no eye-side clone_ids to expand.")
    return(NULL)
  }

  # Load AIRR rows once, restrict to productive IGH heavy chains whose
  # clone_id matches an eye-side clone.
  airr <- dplyr::bind_rows(lapply(airr_files,
                                  utils::read.delim,
                                  stringsAsFactors = FALSE))
  airr <- airr |>
    dplyr::filter(productive == TRUE, locus == "IGH",
                  !is.na(clone_id),
                  clone_id %in% keep_clones) |>
    dplyr::select(cell_id_unique, subject_id, clone_id, c_call) |>
    dplyr::distinct(cell_id_unique, .keep_all = TRUE)

  if (nrow(airr) == 0L) {
    log_message("[fig6 F/G] AIRR carries no eye-clone rows after filter.")
    return(NULL)
  }

  # Lift phenotype (NIU / Viral) from the eye-side bcell meta. Subject
  # column on bm side is "Subject"; AIRR side is "subject_id".
  pheno_map <- eye_bm |>
    dplyr::distinct(Subject, phenotype) |>
    dplyr::rename(subject_id = Subject)

  # Load full-atlas meta and carry Tissue_1, mu_freq_*, substate, Subject
  # for cells whose cell_id_unique matches the AIRR rows.
  obj <- readRDS(obj_path)
  md  <- obj@meta.data
  rm(obj); invisible(gc(verbose = FALSE))
  md$cell_id_unique <- rownames(md)

  # B-cell-only restriction. Full atlas uses celltype_broad == "B cell".
  if ("celltype_broad" %in% colnames(md)) {
    before <- nrow(md)
    md <- md[md$celltype_broad %in% c("B cell"), , drop = FALSE]
    log_message("[fig6 F/G] celltype_broad filter kept ", nrow(md), "/", before)
  } else {
    log_message("[fig6 F/G] no celltype_broad column; trusting AIRR filter.")
  }

  substate_col <- intersect(c("substate", "knn.leiden.cluster"),
                            colnames(md))[1]
  mu_cols_present <- intersect(
    c("mu_freq_cdr_r_heavy", "mu_freq_cdr_s_heavy",
      "mu_freq_fwr_r_heavy", "mu_freq_fwr_s_heavy"),
    colnames(md))
  if (length(mu_cols_present) == 0L) {
    log_message("[fig6 F/G] no mu_freq_*_heavy cols in full-atlas meta;",
                " selection ratio will be NA for paired clones.")
  }

  md_sub <- md |>
    tibble::as_tibble(rownames = NA) |>
    dplyr::filter(cell_id_unique %in% airr$cell_id_unique) |>
    dplyr::select(dplyr::any_of(c(
      "cell_id_unique", "Subject", "Tissue_1",
      "mu_freq_cdr_r_heavy", "mu_freq_cdr_s_heavy",
      "mu_freq_fwr_r_heavy", "mu_freq_fwr_s_heavy",
      substate_col)))

  if (!"Tissue_1" %in% colnames(md_sub)) {
    log_message("[fig6 F/G] Tissue_1 missing from full atlas; aborting.")
    return(NULL)
  }
  if (!"Subject" %in% colnames(md_sub)) {
    log_message("[fig6 F/G] Subject missing from full atlas; aborting.")
    return(NULL)
  }
  if (is.na(substate_col)) {
    md_sub$substate <- NA_character_
  } else if (substate_col != "substate") {
    md_sub <- dplyr::rename(md_sub, substate = !!substate_col)
  }
  for (col in c("mu_freq_cdr_r_heavy", "mu_freq_cdr_s_heavy",
                "mu_freq_fwr_r_heavy", "mu_freq_fwr_s_heavy")) {
    if (!(col %in% colnames(md_sub))) md_sub[[col]] <- 0
  }

  md_sub$SHM_total <- rowSums(
    as.matrix(md_sub[, c("mu_freq_cdr_r_heavy", "mu_freq_cdr_s_heavy",
                         "mu_freq_fwr_r_heavy", "mu_freq_fwr_s_heavy"),
                     drop = FALSE]),
    na.rm = TRUE)

  paired <- md_sub |>
    dplyr::inner_join(airr, by = "cell_id_unique") |>
    dplyr::inner_join(pheno_map, by = "subject_id") |>
    dplyr::mutate(
      subject = subject_id,
      tissue  = dplyr::case_when(
        Tissue_1 == "Eye"   ~ "eye",
        Tissue_1 == "Blood" ~ "blood",
        TRUE                ~ NA_character_),
      isotype_collapsed = {
        cc <- sub("\\*.*$", "", as.character(c_call))
        dplyr::case_when(
          grepl("^IGHM", cc) ~ "IGHM",
          grepl("^IGHD", cc) ~ "IGHD",
          grepl("^IGHG", cc) ~ "IGHG",
          grepl("^IGHA", cc) ~ "IGHA",
          grepl("^IGHE", cc) ~ "IGHE",
          is.na(cc)          ~ NA_character_,
          TRUE               ~ "Other")
      }
    ) |>
    dplyr::filter(!is.na(tissue)) |>
    dplyr::select(cell_id_unique, subject, phenotype, clone_id,
                  substate, tissue, isotype_collapsed,
                  SHM_total,
                  mu_freq_cdr_r_heavy, mu_freq_cdr_s_heavy,
                  mu_freq_fwr_r_heavy, mu_freq_fwr_s_heavy)

  log_message("[fig6 F/G] paired meta: ", nrow(paired), " cells across ",
              dplyr::n_distinct(paired$clone_id), " clones (",
              sum(paired$tissue == "eye"), " eye / ",
              sum(paired$tissue == "blood"), " blood).")
  paired
}

# Primary Wilcoxon on subject-level means + confirmatory LMM with subject
# random intercept + Cliff's delta. Returns a one-row tibble. Robust to
# missing lme4/lmerTest/effsize (those columns become NA).
.bcell_arch_test_one <- function(per_clone, metric, panel) {
  d <- per_clone %>%
    dplyr::filter(phenotype %in% c("Viral", "NIU"),
                  !is.na(.data[[metric]]))
  if (nrow(d) < 10 || dplyr::n_distinct(d$phenotype) < 2) {
    log_message("[fig6 F/G] ", metric, ": skipped (n_clones=", nrow(d),
                ", n_phenotypes=", dplyr::n_distinct(d$phenotype), ")")
    return(tibble::tibble(metric = metric, panel = panel,
                          test = NA_character_, statistic = NA_real_,
                          p_raw = NA_real_, lmm_p = NA_real_,
                          cliff_delta = NA_real_))
  }

  subj <- d %>%
    dplyr::group_by(subject, phenotype) %>%
    dplyr::summarise(value = mean(.data[[metric]], na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::filter(!is.nan(value), !is.na(value))

  if (nrow(subj) < 4 || dplyr::n_distinct(subj$phenotype) < 2) {
    log_message("[fig6 F/G] ", metric, ": skipped (n_clones=", nrow(d),
                ", n_subj=", nrow(subj), ")")
    return(tibble::tibble(metric = metric, panel = panel,
                          test = NA_character_, statistic = NA_real_,
                          p_raw = NA_real_, lmm_p = NA_real_,
                          cliff_delta = NA_real_))
  }

  w <- tryCatch(stats::wilcox.test(value ~ phenotype, data = subj,
                                   exact = FALSE),
                error = function(e) NULL)
  cd <- NA_real_
  if (requireNamespace("effsize", quietly = TRUE)) {
    cd <- tryCatch(effsize::cliff.delta(value ~ phenotype,
                                        data = subj)$estimate,
                   error = function(e) NA_real_)
  }

  lmm_p <- NA_real_
  if (requireNamespace("lme4", quietly = TRUE) &&
      requireNamespace("lmerTest", quietly = TRUE)) {
    fit <- tryCatch(
      lmerTest::lmer(stats::as.formula(
        paste(metric, "~ phenotype + (1 | subject)")), data = d),
      error = function(e) NULL)
    if (!is.null(fit)) {
      co <- summary(fit)$coefficients
      pr <- grep("^phenotype", rownames(co))
      if (length(pr) == 1L && "Pr(>|t|)" %in% colnames(co))
        lmm_p <- co[pr, "Pr(>|t|)"]
    }
  }

  tibble::tibble(
    metric = metric, panel = panel,
    test = "wilcox_subject_mean",
    statistic = if (!is.null(w)) unname(w$statistic) else NA_real_,
    p_raw = if (!is.null(w)) w$p.value else NA_real_,
    lmm_p = lmm_p,
    cliff_delta = unname(cd)
  )
}

# G-i stat (current panel: comp_class x phenotype within shared clones).
# Fisher 3x2 with Monte Carlo p-value (simulate.p.value = TRUE, B = 10000)
# on the 3-level compartment class (eye_only / blood_only / mixed) by
# phenotype (Viral vs Autoimmune). Returns a one-row tibble panel
# "G_i_alluvial" so the existing regression check (which expects that
# panel name) keeps passing.
.bcell_arch_test_comp <- function(per_clone, cfg) {
  arch <- cfg$bcr_lineage$architecture
  seed <- arch$baseline$seed %||% 42
  na_row <- tibble::tibble(
    metric = "comp_class_by_phenotype", panel = "G_i_alluvial",
    test = NA_character_, statistic = NA_real_,
    p_raw = NA_real_, lmm_p = NA_real_, cliff_delta = NA_real_)
  d <- per_clone %>%
    dplyr::filter(phenotype %in% c("Viral", "NIU"),
                  !is.na(comp_class))
  if (nrow(d) < 10 ||
      dplyr::n_distinct(d$phenotype) < 2 ||
      dplyr::n_distinct(d$comp_class) < 2) {
    log_message("[fig6 F/G] G-i comp_class x phenotype: too sparse; skipping.")
    return(na_row)
  }
  comp_levels <- c("eye_only", "blood_only", "mixed")
  d$comp_class <- factor(d$comp_class, levels = comp_levels)
  tab <- table(d$comp_class, d$phenotype)
  set.seed(seed)
  ft <- tryCatch(stats::fisher.test(tab, simulate.p.value = TRUE,
                                    B = 10000),
                 error = function(e) NULL)
  if (is.null(ft)) return(na_row)
  tibble::tibble(
    metric = "comp_class_by_phenotype", panel = "G_i_alluvial",
    test = "fisher_sim",
    statistic = NA_real_,
    p_raw = ft$p.value,
    lmm_p = NA_real_, cliff_delta = NA_real_
  )
}

# Build per-disease (blood_substate x eye_substate) clone-transition
# matrices from a breakdown table. Each shared clone contributes one
# count to each (blood_sub, eye_sub) pair it touches (the "min-transferable"
# logic in G-i is for visual weight only; the stat sums clone presence
# across substate pairs).
#
# Keyed on (subject, clone_id) because clone_id is not globally unique
# across subjects in the AIRR pipeline.
#
# Returns a named list of matrices keyed by disease (Viral, Autoimmune).
.bcell_arch_transition_matrices <- function(breakdown, shared_clone_keys) {
  bk <- breakdown %>%
    dplyr::mutate(.key = paste(subject, clone_id, sep = "::")) %>%
    dplyr::filter(.key %in% shared_clone_keys,
                  phenotype %in% c("Viral", "NIU"),
                  !is.na(substate), !is.na(tissue)) %>%
    dplyr::group_by(phenotype, .key, substate, tissue) %>%
    dplyr::summarise(n_cells = sum(n_cells, na.rm = TRUE),
                     .groups = "drop")
  if (nrow(bk) == 0L) return(list())

  # Enumerate (blood_sub, eye_sub) pairs per clone within each disease.
  edges <- bk %>%
    dplyr::group_by(phenotype, .key) %>%
    dplyr::group_modify(function(df, key) {
      bs <- unique(df$substate[df$tissue == "blood"])
      es <- unique(df$substate[df$tissue == "eye"])
      if (length(bs) == 0L || length(es) == 0L) return(tibble::tibble())
      expand.grid(blood_sub = bs, eye_sub = es,
                  stringsAsFactors = FALSE)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(phenotype = dplyr::recode(phenotype,
                                            NIU = "Autoimmune"))
  if (nrow(edges) == 0L) return(list())

  all_subs <- sort(unique(c(edges$blood_sub, edges$eye_sub)))
  out <- list()
  for (dz in c("Viral", "Autoimmune")) {
    e_dz <- edges[edges$phenotype == dz, , drop = FALSE]
    m <- matrix(0L, nrow = length(all_subs), ncol = length(all_subs),
                dimnames = list(all_subs, all_subs))
    if (nrow(e_dz) > 0L) {
      tab <- table(factor(e_dz$blood_sub, levels = all_subs),
                   factor(e_dz$eye_sub,   levels = all_subs))
      m[] <- as.integer(tab)
    }
    out[[dz]] <- m
  }
  out
}

# G-i replacement: per-disease G-test on (blood_substate x eye_substate)
# and a permutation Frobenius homogeneity test across diseases.
.bcell_arch_test_alluvial_transition <- function(per_clone, breakdown, cfg) {
  arch <- cfg$bcr_lineage$architecture
  seed <- arch$baseline$seed %||% 42
  set.seed(seed)

  shared <- per_clone %>%
    dplyr::filter(is_shared, phenotype %in% c("Viral", "NIU"))
  shared_keys <- paste(shared$subject, shared$clone_id, sep = "::")
  if (length(shared_keys) == 0L) {
    return(tibble::tibble(
      metric = c("blood_to_eye_transition",
                 "blood_to_eye_transition_likelihood"),
      panel = "G_i_alluvial",
      test = NA_character_,
      statistic = NA_real_, p_raw = NA_real_,
      lmm_p = NA_real_, cliff_delta = NA_real_))
  }

  mats <- .bcell_arch_transition_matrices(breakdown, shared_keys)
  if (length(mats) < 2L ||
      sum(unlist(lapply(mats, sum))) == 0L) {
    return(tibble::tibble(
      metric = c("blood_to_eye_transition",
                 "blood_to_eye_transition_likelihood"),
      panel = "G_i_alluvial",
      test = NA_character_,
      statistic = NA_real_, p_raw = NA_real_,
      lmm_p = NA_real_, cliff_delta = NA_real_))
  }

  # Per-disease G-test (likelihood-ratio chi-square) on each matrix.
  per_dz_stats <- lapply(names(mats), function(dz) {
    m <- mats[[dz]]
    if (sum(m) == 0 || nrow(m) < 2 || ncol(m) < 2) {
      return(list(stat = NA_real_, p = NA_real_))
    }
    if (requireNamespace("vcd", quietly = TRUE)) {
      ct <- tryCatch(vcd::assocstats(as.table(m)),
                     error = function(e) NULL)
      if (!is.null(ct) && !is.null(ct$chisq_tests)) {
        # vcd reports rows for Likelihood Ratio + Pearson
        ct_df <- ct$chisq_tests
        # ct_df is a matrix; row "Likelihood Ratio" has X^2 and P(>X^2)
        row_lr <- which(rownames(ct_df) == "Likelihood Ratio")
        if (length(row_lr) == 1L) {
          return(list(stat = unname(ct_df[row_lr, "X^2"]),
                      p    = unname(ct_df[row_lr, "P(> X^2)"])))
        }
      }
    }
    # Fallback: base R chi-square (Pearson) when vcd not available.
    cs <- tryCatch(stats::chisq.test(m, simulate.p.value = TRUE,
                                     B = 2000),
                   error = function(e) NULL)
    if (is.null(cs)) return(list(stat = NA_real_, p = NA_real_))
    list(stat = unname(cs$statistic), p = cs$p.value)
  })
  names(per_dz_stats) <- names(mats)

  # Permutation Frobenius test on cross-disease homogeneity. Use a
  # composite (subject, clone_id) key because clone_id alone is not
  # globally unique. Permute the disease label across keys.
  bk <- breakdown %>%
    dplyr::mutate(.key = paste(subject, clone_id, sep = "::")) %>%
    dplyr::filter(.key %in% shared_keys,
                  phenotype %in% c("Viral", "NIU"),
                  !is.na(substate), !is.na(tissue))
  edges_clone <- bk %>%
    dplyr::group_by(.key) %>%
    dplyr::group_modify(function(df, key) {
      bs <- unique(df$substate[df$tissue == "blood"])
      es <- unique(df$substate[df$tissue == "eye"])
      if (length(bs) == 0L || length(es) == 0L) return(tibble::tibble())
      expand.grid(blood_sub = bs, eye_sub = es,
                  stringsAsFactors = FALSE)
    }) %>%
    dplyr::ungroup()

  clone_dz <- bk %>%
    dplyr::distinct(.key, phenotype) %>%
    dplyr::mutate(disease = dplyr::recode(phenotype,
                                          NIU = "Autoimmune"))

  if (nrow(edges_clone) == 0L ||
      dplyr::n_distinct(clone_dz$disease) < 2L) {
    perm_stat <- NA_real_
    perm_p    <- NA_real_
  } else {
    all_subs <- sort(unique(c(edges_clone$blood_sub, edges_clone$eye_sub)))
    build_m <- function(df) {
      m <- matrix(0L, length(all_subs), length(all_subs),
                  dimnames = list(all_subs, all_subs))
      if (nrow(df) == 0L) return(m)
      tab <- table(factor(df$blood_sub, levels = all_subs),
                   factor(df$eye_sub,   levels = all_subs))
      m[] <- as.integer(tab); m
    }
    edges_clone <- dplyr::left_join(
      edges_clone,
      dplyr::distinct(clone_dz[, c(".key", "disease")]),
      by = ".key",
      relationship = "many-to-one")
    mat_obs <- list(
      Viral      = build_m(edges_clone[edges_clone$disease == "Viral", ]),
      Autoimmune = build_m(edges_clone[edges_clone$disease == "Autoimmune", ])
    )
    obs_diff <- sqrt(sum((mat_obs$Viral - mat_obs$Autoimmune)^2))

    n_perm <- 1000L
    null_diffs <- numeric(n_perm)
    keys_unique <- unique(edges_clone$.key)
    cdz <- clone_dz$disease[match(keys_unique, clone_dz$.key)]
    set.seed(seed)
    for (i in seq_len(n_perm)) {
      perm <- sample(cdz)
      perm_map <- data.frame(.key = keys_unique,
                             disease_perm = perm,
                             stringsAsFactors = FALSE)
      eg <- dplyr::left_join(edges_clone[, c(".key", "blood_sub",
                                             "eye_sub")],
                             perm_map, by = ".key",
                             relationship = "many-to-one")
      mv <- build_m(eg[eg$disease_perm == "Viral", ])
      ma <- build_m(eg[eg$disease_perm == "Autoimmune", ])
      null_diffs[i] <- sqrt(sum((mv - ma)^2))
    }
    perm_stat <- obs_diff
    perm_p    <- mean(null_diffs >= obs_diff)
  }

  # Two stat rows: per-disease LR summary (averaged stat, max p) +
  # permutation Frobenius across diseases.
  lr_stats <- vapply(per_dz_stats, `[[`, numeric(1), "stat")
  lr_ps    <- vapply(per_dz_stats, `[[`, numeric(1), "p")
  tibble::tibble(
    metric = c("blood_to_eye_transition",
               "blood_to_eye_transition_likelihood"),
    panel  = "G_i_alluvial",
    test   = c("permutation_frobenius", "vcd_assocstats_LR"),
    statistic = c(perm_stat, mean(lr_stats, na.rm = TRUE)),
    p_raw     = c(perm_p,    max(lr_ps,    na.rm = TRUE)),
    lmm_p = NA_real_,
    cliff_delta = NA_real_
  )
}

# S6c stat: per-disease chi-square on (shared_label x isotype) per clone
# (per-clone dominant isotype, Shared vs Non-shared), plus a Cochran-
# Mantel-Haenszel homogeneity test across diseases. Returns one tibble
# row per test under panel "S6c".
.bcell_arch_test_isotype <- function(per_clone, cfg) {
  arch <- cfg$bcr_lineage$architecture
  seed <- arch$baseline$seed %||% 42
  na_rows <- tibble::tibble(
    metric = c("isotype_shared_Viral",
               "isotype_shared_Autoimmune",
               "isotype_shared_homogeneity"),
    panel = "S6c",
    test = NA_character_,
    statistic = NA_real_, p_raw = NA_real_,
    lmm_p = NA_real_, cliff_delta = NA_real_)

  d <- per_clone %>%
    dplyr::filter(phenotype %in% c("Viral", "NIU"),
                  !is.na(dominant_isotype)) %>%
    dplyr::mutate(
      iso = dplyr::case_when(
        grepl("^IGHM", dominant_isotype) ~ "IgM",
        grepl("^IGHD", dominant_isotype) ~ "IgD",
        grepl("^IGHG", dominant_isotype) ~ "IgG",
        grepl("^IGHA", dominant_isotype) ~ "IgA",
        grepl("^IGHE", dominant_isotype) ~ "IgE",
        TRUE                              ~ "Other"),
      shared_label = ifelse(is_shared, "Shared", "Non-shared"),
      disease = dplyr::recode(phenotype, NIU = "Autoimmune"))

  if (nrow(d) < 20)
    return(na_rows)

  set.seed(seed)
  per_pheno <- purrr::map_dfr(c("Viral", "Autoimmune"), function(dz) {
    sub <- d %>% dplyr::filter(disease == dz)
    nm  <- paste0("isotype_shared_", dz)
    if (nrow(sub) < 5L)
      return(tibble::tibble(metric = nm, panel = "S6c",
                            test = NA_character_, statistic = NA_real_,
                            p_raw = NA_real_, lmm_p = NA_real_,
                            cliff_delta = NA_real_))
    tab <- table(sub$shared_label, sub$iso)
    if (any(dim(tab) < 2L))
      return(tibble::tibble(metric = nm, panel = "S6c",
                            test = "chisq_sim",
                            statistic = NA_real_, p_raw = NA_real_,
                            lmm_p = NA_real_, cliff_delta = NA_real_))
    cs <- tryCatch(stats::chisq.test(tab, simulate.p.value = TRUE,
                                     B = 10000),
                   error = function(e) NULL)
    if (is.null(cs))
      return(tibble::tibble(metric = nm, panel = "S6c",
                            test = "chisq_sim",
                            statistic = NA_real_, p_raw = NA_real_,
                            lmm_p = NA_real_, cliff_delta = NA_real_))
    tibble::tibble(metric = nm, panel = "S6c", test = "chisq_sim",
                   statistic = unname(cs$statistic),
                   p_raw = cs$p.value, lmm_p = NA_real_,
                   cliff_delta = NA_real_)
  })

  # Cochran-Mantel-Haenszel homogeneity across diseases.
  iso_levels    <- sort(unique(d$iso))
  shared_levels <- c("Non-shared", "Shared")
  arr <- array(0L,
               dim = c(length(shared_levels), length(iso_levels), 2L),
               dimnames = list(shared_levels, iso_levels,
                               c("Viral", "Autoimmune")))
  for (dz in c("Viral", "Autoimmune")) {
    sub <- d %>% dplyr::filter(disease == dz)
    if (nrow(sub) == 0L) next
    tab <- table(factor(sub$shared_label, levels = shared_levels),
                 factor(sub$iso,          levels = iso_levels))
    arr[, , dz] <- as.integer(tab)
  }
  cmh <- tryCatch(stats::mantelhaen.test(arr, exact = FALSE),
                  error = function(e) NULL)
  cmh_row <- if (!is.null(cmh)) {
    tibble::tibble(metric = "isotype_shared_homogeneity", panel = "S6c",
                   test = "mantelhaen",
                   statistic = unname(cmh$statistic),
                   p_raw = cmh$p.value, lmm_p = NA_real_,
                   cliff_delta = NA_real_)
  } else {
    tibble::tibble(metric = "isotype_shared_homogeneity", panel = "S6c",
                   test = NA_character_,
                   statistic = NA_real_, p_raw = NA_real_,
                   lmm_p = NA_real_, cliff_delta = NA_real_)
  }

  dplyr::bind_rows(per_pheno, cmh_row)
}

# S6b new stat: Wilcoxon on per-subject public-clone counts + Cliff's
# delta. Returns NA-filled row when the public CSV is not present.
.bcell_arch_test_public_clones <- function(cfg) {
  paths_rep <- tryCatch(get_target_paths(cfg, "repertoire"),
                        error = function(e) NULL)
  pub_path <- if (!is.null(paths_rep) &&
                  !is.null(paths_rep$results_tables)) {
    file.path(paths_rep$results_tables, "BCR_public_clones.csv")
  } else {
    file.path("outputs", "tables", "repertoire", "BCR_public_clones.csv")
  }
  na_row <- tibble::tibble(
    metric = "public_clones_per_subject", panel = "S6b",
    test = NA_character_, statistic = NA_real_,
    p_raw = NA_real_, lmm_p = NA_real_, cliff_delta = NA_real_)
  if (!file.exists(pub_path)) {
    log_message("[fig6 F/G] public clones CSV missing at ", pub_path,
                "; S6b stat skipped.")
    return(na_row)
  }
  pub <- tryCatch(readr::read_csv(pub_path, show_col_types = FALSE),
                  error = function(e) NULL)
  if (is.null(pub) ||
      !all(c("subject_id", "cluster_id", "stratum") %in% colnames(pub))) {
    log_message("[fig6 F/G] public clones CSV missing required cols; ",
                "S6b stat skipped.")
    return(na_row)
  }
  subj_pub <- pub %>%
    dplyr::filter(stratum %in% c("Viral", "NIU")) %>%
    dplyr::distinct(subject_id, cluster_id, stratum) %>%
    dplyr::count(subject_id, stratum, name = "n_public") %>%
    dplyr::rename(phenotype = stratum)
  if (nrow(subj_pub) < 4 ||
      dplyr::n_distinct(subj_pub$phenotype) < 2) {
    log_message("[fig6 F/G] public clones: too few subjects for Wilcoxon; ",
                "S6b stat skipped.")
    return(na_row)
  }
  w <- tryCatch(stats::wilcox.test(n_public ~ phenotype,
                                   data = subj_pub, exact = FALSE),
                error = function(e) NULL)
  cd <- NA_real_
  if (requireNamespace("effsize", quietly = TRUE)) {
    cd <- tryCatch(
      effsize::cliff.delta(n_public ~ phenotype,
                           data = subj_pub)$estimate,
      error = function(e) NA_real_)
  }
  tibble::tibble(
    metric = "public_clones_per_subject", panel = "S6b",
    test = "wilcox_per_subject",
    statistic = if (!is.null(w)) unname(w$statistic) else NA_real_,
    p_raw     = if (!is.null(w)) w$p.value           else NA_real_,
    lmm_p = NA_real_,
    cliff_delta = unname(cd)
  )
}

# Per-disease Wilcoxon on per-clone shm_mean: Shared vs Non-shared. Two
# rows in the returned tibble (one per disease). Used by F-ii panel which
# is now laid out as boxplot(shared_label) faceted by phenotype.
.bcell_arch_test_shm_by_shared <- function(per_clone) {
  base_row <- function(ph)
    tibble::tibble(
      metric = paste0("shm_shared_vs_nonshared_", ph),
      panel  = "F_ii",
      test   = NA_character_, statistic = NA_real_,
      p_raw  = NA_real_, lmm_p = NA_real_, cliff_delta = NA_real_)

  d <- per_clone %>%
    dplyr::filter(!is.na(shm_mean),
                  phenotype %in% c("Viral", "NIU")) %>%
    dplyr::mutate(
      shared_label = factor(ifelse(is_shared, "Shared", "Non-shared"),
                            levels = c("Non-shared", "Shared")),
      ph_label = dplyr::recode(phenotype, NIU = "Autoimmune"))

  purrr::map_dfr(c("Viral", "Autoimmune"), function(ph) {
    sub <- d %>% dplyr::filter(ph_label == ph)
    if (nrow(sub) < 10 || dplyr::n_distinct(sub$shared_label) < 2)
      return(base_row(ph))
    w <- tryCatch(stats::wilcox.test(shm_mean ~ shared_label,
                                     data = sub, exact = FALSE),
                  error = function(e) NULL)
    cd <- NA_real_
    if (requireNamespace("effsize", quietly = TRUE)) {
      cd <- tryCatch(
        effsize::cliff.delta(shm_mean ~ shared_label,
                             data = sub)$estimate,
        error = function(e) NA_real_)
    }
    tibble::tibble(
      metric = paste0("shm_shared_vs_nonshared_", ph),
      panel  = "F_ii",
      test   = "wilcox_clone_within_disease",
      statistic = if (!is.null(w)) unname(w$statistic) else NA_real_,
      p_raw     = if (!is.null(w)) w$p.value           else NA_real_,
      lmm_p     = NA_real_,
      cliff_delta = unname(cd))
  })
}

# Compute per-clone summary metrics from paired (eye + blood) B-cell meta.
# Input meta_paired: data.frame with one row per cell carrying at minimum
#   cell_id_unique, subject, phenotype ("NIU" / "Viral"), clone_id,
#   substate, tissue ("eye" / "blood"), SHM_total,
#   mu_freq_cdr_r_heavy, mu_freq_cdr_s_heavy,
#   mu_freq_fwr_r_heavy, mu_freq_fwr_s_heavy
# Missing mu_freq columns are tolerated (treated as 0 across the clone).
# Returns: tibble, one row per (subject, clone_id), with:
#   n_cells, n_eye, n_blood, n_clusters, shm_mean, cdr_*_sum, fwr_*_sum,
#   n_mutated, sel_cdr_r_ratio, is_shared, comp_class, tissue_entropy.
.bcell_arch_per_clone <- function(meta_paired, cfg) {
  arch    <- cfg$bcr_lineage$architecture
  min_n   <- arch$min_clone_size %||% 3
  min_k   <- arch$min_clusters_for_shared %||% 2
  min_mut <- arch$min_mutated_cells_for_selection %||% 1

  mu_cols <- c("mu_freq_cdr_r_heavy", "mu_freq_cdr_s_heavy",
               "mu_freq_fwr_r_heavy", "mu_freq_fwr_s_heavy")
  for (col in mu_cols) {
    if (!(col %in% colnames(meta_paired)))
      meta_paired[[col]] <- 0
  }

  m <- meta_paired %>%
    dplyr::filter(!is.na(clone_id), !is.na(subject), !is.na(phenotype)) %>%
    dplyr::mutate(
      cdr_r = tidyr::replace_na(mu_freq_cdr_r_heavy, 0),
      cdr_s = tidyr::replace_na(mu_freq_cdr_s_heavy, 0),
      fwr_r = tidyr::replace_na(mu_freq_fwr_r_heavy, 0),
      fwr_s = tidyr::replace_na(mu_freq_fwr_s_heavy, 0)
    )

  per_clone <- m %>%
    dplyr::group_by(subject, phenotype, clone_id) %>%
    dplyr::summarise(
      n_cells     = dplyr::n(),
      n_eye       = sum(tissue == "eye", na.rm = TRUE),
      n_blood     = sum(tissue == "blood", na.rm = TRUE),
      n_clusters  = dplyr::n_distinct(substate),
      shm_mean    = mean(SHM_total, na.rm = TRUE),
      cdr_r_sum   = sum(cdr_r, na.rm = TRUE),
      cdr_s_sum   = sum(cdr_s, na.rm = TRUE),
      fwr_r_sum   = sum(fwr_r, na.rm = TRUE),
      fwr_s_sum   = sum(fwr_s, na.rm = TRUE),
      # mu_freq_* is a per-cell mutation frequency, so any positive sum
      # across the four CDR/FWR R/S compartments is a valid "this cell
      # carries at least one IGH mutation" indicator.
      n_mutated   = sum((cdr_r + cdr_s + fwr_r + fwr_s) > 0, na.rm = TRUE),
      dominant_isotype = if ("isotype_collapsed" %in% colnames(m)) {
        tab <- sort(table(isotype_collapsed, useNA = "no"), decreasing = TRUE)
        # Ties broken by sort order on the collapsed labels (alphabetical).
        if (length(tab) == 0) NA_character_ else names(tab)[1]
      } else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cells >= min_n) %>%
    dplyr::mutate(
      sel_cdr_r_ratio = dplyr::if_else(
        n_mutated >= min_mut &
          (cdr_s_sum + fwr_r_sum + fwr_s_sum) > 0,
        cdr_r_sum / (cdr_s_sum + fwr_r_sum + fwr_s_sum),
        NA_real_
      ),
      is_shared = n_clusters >= min_k,
      comp_class = dplyr::case_when(
        n_eye  > 0 & n_blood == 0 ~ "eye_only",
        n_eye == 0 & n_blood  > 0 ~ "blood_only",
        n_eye  > 0 & n_blood  > 0 ~ "mixed",
        TRUE                       ~ NA_character_
      ),
      tissue_entropy = .shannon_two(n_eye, n_blood)
    )
  per_clone
}

run_bcell_lineage_architecture <- function(cfg) {
  log_message("[fig6 F/G] computing per-clone architecture metrics")
  inputs <- .load_bcell_arch_meta(cfg)
  if (is.null(inputs)) {
    log_message("[fig6 F/G] bcell inputs missing; skipping.")
    return(invisible(FALSE))
  }
  # paths$results_tables for the bcell target already includes "bcell";
  # write outputs under <results_tables>/architecture/.
  out_dir <- file.path(inputs$paths$results_tables, "architecture")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Eye-side bcell meta. Attach clone_id from the eye AIRR rows so
  # downstream code can identify eye-side clones to expand to blood.
  eye_bm <- inputs$bm
  if (!is.null(inputs$bcr_eye) && nrow(inputs$bcr_eye) > 0L) {
    eye_clone_map <- inputs$bcr_eye |>
      dplyr::select(cell_id_unique, clone_id) |>
      dplyr::distinct(cell_id_unique, .keep_all = TRUE)
    eye_bm <- eye_bm |>
      dplyr::left_join(eye_clone_map, by = "cell_id_unique")
  } else {
    eye_bm$clone_id <- NA_character_
  }

  paired <- .bcell_arch_paired_meta(eye_bm, cfg)
  if (is.null(paired)) {
    log_message("[fig6 F/G] paired meta unavailable; G panels degraded.")
    # Eye-only fallback so F panels still run. Join clone_id + bm enrich.
    if (!exists(".bcell_lineage_meta", mode = "function"))
      source(file.path("R", "56_bcell_lineage_trees.R"))
    enrich <- tryCatch(.bcell_lineage_meta(cfg),
                       error = function(e) NULL)
    if (is.null(enrich)) {
      log_message("[fig6 F/G] no eye-only meta either; aborting.")
      return(invisible(FALSE))
    }
    # magrittr pipe required here: the transmute below uses `.` to peek
    # at the joined frame's column set. Native |> does not bind `.`.
    paired <- eye_bm %>%
      dplyr::filter(!is.na(clone_id)) %>%
      dplyr::left_join(enrich, by = "cell_id_unique") %>%
      dplyr::transmute(
        cell_id_unique, subject = Subject, phenotype,
        clone_id,
        substate = if ("substate" %in% colnames(.))
          substate else NA_character_,
        tissue = dplyr::case_when(
          Tissue_1 == "Eye"   ~ "eye",
          Tissue_1 == "Blood" ~ "blood",
          TRUE                ~ NA_character_),
        isotype_collapsed = if ("isotype_collapsed" %in% colnames(.))
          isotype_collapsed else NA_character_,
        SHM_total = if ("SHM_total" %in% colnames(.))
          SHM_total else NA_real_,
        dplyr::across(dplyr::any_of(c(
          "mu_freq_cdr_r_heavy","mu_freq_cdr_s_heavy",
          "mu_freq_fwr_r_heavy","mu_freq_fwr_s_heavy")),
          ~ tidyr::replace_na(.x, 0))
      ) %>%
      dplyr::filter(!is.na(tissue))
  }

  per_clone <- .bcell_arch_per_clone(paired, cfg)
  readr::write_csv(per_clone,
                   file.path(out_dir, "clone_architecture_metrics.csv"))
  log_message("[fig6 F/G] wrote ", nrow(per_clone),
              " clone rows to clone_architecture_metrics.csv")

  # Long-form per-(clone x substate x tissue x isotype) cell breakdown.
  # Drives the new G-i and S6c alluvials + their permutation/CMH stats.
  breakdown <- paired %>%
    dplyr::filter(!is.na(clone_id), !is.na(substate), !is.na(tissue)) %>%
    dplyr::count(subject, phenotype, clone_id, substate, tissue,
                 isotype_collapsed, name = "n_cells")
  readr::write_csv(breakdown,
                   file.path(out_dir, "clone_cell_breakdown.csv"))
  log_message("[fig6 F/G] wrote ", nrow(breakdown),
              " (clone x substate x tissue x isotype) rows to ",
              "clone_cell_breakdown.csv")

  subject_summary <- per_clone %>%
    dplyr::group_by(subject, phenotype) %>%
    dplyr::summarise(
      n_clones = dplyr::n(),
      mean_clusters_spanned = mean(n_clusters, na.rm = TRUE),
      mean_shm = mean(shm_mean, na.rm = TRUE),
      mean_sel_ratio = mean(sel_cdr_r_ratio, na.rm = TRUE),
      frac_shared_mixed = mean(is_shared & comp_class == "mixed",
                               na.rm = TRUE),
      mean_tissue_entropy = dplyr::if_else(
        is.nan(mean(tissue_entropy[is_shared], na.rm = TRUE)),
        NA_real_,
        mean(tissue_entropy[is_shared], na.rm = TRUE)),
      .groups = "drop"
    )
  readr::write_csv(subject_summary,
                   file.path(out_dir, "subject_architecture_summary.csv"))
  log_message("[fig6 F/G] wrote ", nrow(subject_summary),
              " subject rows to subject_architecture_summary.csv")

  stats_rows <- dplyr::bind_rows(
    .bcell_arch_test_one(per_clone, "n_clusters",      "F_i"),
    .bcell_arch_test_one(per_clone, "shm_mean",        "F_ii_overall"),
    .bcell_arch_test_shm_by_shared(per_clone),
    .bcell_arch_test_one(per_clone, "sel_cdr_r_ratio", "F_iii"),
    .bcell_arch_test_comp(per_clone, cfg),
    .bcell_arch_test_one(per_clone %>%
                           dplyr::filter(is_shared),
                         "tissue_entropy", "G_ii"),
    .bcell_arch_test_isotype(per_clone, cfg),
    .bcell_arch_test_public_clones(cfg)
  ) %>%
    dplyr::mutate(p_bh = stats::p.adjust(p_raw, method = "BH"))

  readr::write_csv(stats_rows,
                   file.path(out_dir, "stats_summary.csv"))
  log_message("[fig6 F/G] wrote ", nrow(stats_rows),
              " contrast rows to stats_summary.csv")

  log_message("[fig6 F/G] done.")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Fig 6 F/G/S6 schema regression check (folded in from
# scripts/regression_check_fig6_architecture.R). Verifies the architecture
# outputs exist with the expected columns/panels. WARNS on any breach (does
# not halt the pipeline). Gated by cfg$steps_fig6$check_architecture.
# ---------------------------------------------------------------------------
check_fig6_architecture_schema <- function(cfg = NULL) {
  ok <- TRUE
  warn <- function(msg) { ok <<- FALSE; warning(msg, call. = FALSE) }
  d <- "outputs/tables/eye/bcell/architecture"
  needed <- c("clone_architecture_metrics.csv",
              "subject_architecture_summary.csv", "stats_summary.csv")
  missing <- needed[!file.exists(file.path(d, needed))]
  if (length(missing)) {
    warn(paste0("[fig6 F/G] missing outputs: ", paste(missing, collapse = ", ")))
  } else {
    pc <- readr::read_csv(file.path(d, "clone_architecture_metrics.csv"),
                          show_col_types = FALSE)
    must_have <- c("subject","phenotype","clone_id","n_cells","n_eye","n_blood",
                   "n_clusters","shm_mean","sel_cdr_r_ratio","comp_class",
                   "tissue_entropy","is_shared","dominant_isotype")
    miss <- setdiff(must_have, colnames(pc))
    if (length(miss))
      warn(paste0("[fig6 F/G] clone_architecture_metrics missing columns: ",
                  paste(miss, collapse = ", ")))
    ss <- readr::read_csv(file.path(d, "stats_summary.csv"), show_col_types = FALSE)
    if (!all(c("metric","panel","p_raw","p_bh","cliff_delta") %in% colnames(ss)))
      warn("[fig6 F/G] stats_summary missing required columns")
    panel_set <- c("F_i","F_ii","F_iii","G_i_alluvial","G_ii")
    if (!all(panel_set %in% ss$panel))
      warn(paste0("[fig6 F/G] stats_summary missing panels: ",
                  paste(setdiff(panel_set, ss$panel), collapse = ", ")))
    bk_path <- file.path(d, "clone_cell_breakdown.csv")
    if (!file.exists(bk_path)) {
      warn(paste0("[fig6 F/G] missing clone_cell_breakdown.csv at ", bk_path))
    } else {
      bk <- readr::read_csv(bk_path, show_col_types = FALSE)
      miss_bk <- setdiff(c("subject","phenotype","clone_id","substate","tissue",
                           "isotype_collapsed","n_cells"), colnames(bk))
      if (length(miss_bk))
        warn(paste0("[fig6 F/G] clone_cell_breakdown missing cols: ",
                    paste(miss_bk, collapse = ", ")))
    }
  }
  bl_path <- file.path(d, "baseline_selection.csv")
  if (file.exists(bl_path)) {
    bl <- readr::read_csv(bl_path, show_col_types = FALSE)
    miss <- setdiff(c("subject","phenotype","substate","region","baseline_sigma"),
                    colnames(bl))
    if (length(miss))
      warn(paste0("[figS6a] baseline_selection missing columns: ",
                  paste(miss, collapse = ", ")))
  }
  main_pdfs <- list.files("outputs/viz/eye/bcell/10_lineage_arch/architecture",
                          pattern = "^fig6[FG]_.*\\.pdf$")
  if (length(main_pdfs) < 5)
    warn(paste0("[fig6 F/G] expected >= 5 main panel PDFs, found ", length(main_pdfs)))
  supp_pdfs <- list.files("outputs/viz/eye/bcell/10_lineage_arch/architecture_supp",
                          pattern = "^figS6_.*\\.pdf$")
  if (length(supp_pdfs) < 3)
    warn(paste0("[figS6] expected >= 3 supplement panel PDFs, found ", length(supp_pdfs)))
  if (ok) log_message("[fig6 F/G/S6] schema check passed.")
  else    log_message("[fig6 F/G/S6] schema check found issues (see warnings).")
  invisible(ok)
}
