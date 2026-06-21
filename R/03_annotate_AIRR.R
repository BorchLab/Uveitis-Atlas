# R/03_annotate_AIRR.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  # AIRR/immcantation stack
  library(alakazam)
  library(dowser)
  library(scoper)
  library(shazam)
  library(scRepertoire)
  library(ggtree)
  library(ggplot2)
})

# ------------------------- Configuration helpers -------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

safe_readRDS <- function(path) {
  if (!file.exists(path)) stop("Missing RDS: ", path)
  readRDS(path)
}

# Normalize scRepertoire barcodes to Seurat cell names.
# scRepertoire's combineTCR/combineBCR prefix the sample with an underscore:
#   "SAMPLE_AAACCGTG-1"
# Seurat RenameCells (in 10_ingest_data.R) produces "SAMPLE.AAACCGTG-1".
# Replace the LAST underscore with a dot so both sides match.
normalize_barcode <- function(barcode) {
  sub("_(?=[^_]*$)", ".", barcode, perl = TRUE)
}

# Fill missing columns so rbind/bind_rows works
fill_missing_cols <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) df[missing] <- NA
  df[, cols, drop = FALSE]
}

# Clone size binning (standardized, readable)
compute_clone_bins <- function(clone_df) {
  # clone_df: columns CTstrict, clonalFrequency
  clone_df |>
    group_by(class) |> 
    mutate(clonalProportion = clonalFrequency / sum(clonalFrequency)) |>
    mutate(
      cloneSize = cut(
        clonalProportion,
        breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 1e-1, Inf),
        labels = c(
          "Rare (0 < X <= 1e-4)",
          "Small (1e-4 < X <= 0.001)",
          "Medium (0.001 < X <= 0.01)",
          "Large (0.01 < X <= 0.1)",
          "Hyperexpanded (0.1 < X <= 1)"
        ),
        right = TRUE
      )
    )
}

# --- Robust threshold finder for BCR heavy-chain distances -------------------
robust_find_threshold <- function(
    x,
    spc_target = 0.995,
    gmm_timeout = 120,        # seconds
    max_n = 200000,
    top_trim = 0.999,
    verbose = TRUE,
    seed = 1L,
    engine = c("auto", "callr_bg", "fork", "inline"),
    density_bw = "nrd0",
    density_n  = 2048
) {
  stopifnot(is.numeric(x))
  engine <- match.arg(engine)
  
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # 1) sanitize
  x <- x[is.finite(x)]
  x <- x[x > 0]
  if (!length(x)) stop("No finite, positive distances to threshold.")
  
  # 2) trim extreme tail
  qtop <- stats::quantile(x, probs = top_trim, names = FALSE, type = 7)
  x <- pmin(x, qtop)
  
  # 3) downsample if huge
  if (length(x) > max_n) {
    set.seed(seed)
    x <- x[sample.int(length(x), max_n)]
  }
  
  # Worker (runs in child)
  gmm_worker <- function(x, spc_target) {
    ft <- shazam::findThreshold(
      x,
      progress = FALSE,
      method   = "gmm",
      model    = "gamma-norm",
      cutoff   = "optimal",
      spc      = spc_target
    )
    as.numeric(ft@threshold)
  }
  
  # Engine impls
  run_gmm_callr_bg <- function(x, spc_target, timeout) {
    if (!requireNamespace("callr", quietly = TRUE)) {
      stop("Package 'callr' is required for engine='callr_bg'. Please install.packages('callr').")
    }
    
    # Normalize timeout: processx wants milliseconds; -1 means infinite
    ms <- if (is.finite(timeout)) as.integer(max(0, timeout) * 1000L) else -1L
    
    p <- callr::r_bg(
      func = function(xx, spc) {
        ft <- shazam::findThreshold(
          xx,
          progress = FALSE,
          method   = "gmm",
          model    = "gamma-norm",
          cutoff   = "optimal",
          spc      = spc
        )
        as.numeric(ft@threshold)
      },
      args = list(xx = x, spc = spc_target),
      supervise = TRUE
    )
    
    # Be defensive about the return type of wait()
    done <- tryCatch(p$wait(ms), error = function(e) e)
    
    # If wait errored or did not return TRUE, treat as timeout
    if (!isTRUE(done)) {
      # kill the whole process tree; ignore errors
      try(p$kill(tree = TRUE, grace = 0), silent = TRUE)
      return(list(ok = FALSE, thr = NA_real_, timed_out = TRUE, err = "timeout"))
    }
    
    # If the process finished, check for an error raised in the child
    pe <- p$get_error()
    if (!is.null(pe)) {
      return(list(ok = FALSE, thr = NA_real_, timed_out = FALSE,
                  err = conditionMessage(pe)))
    }
    
    thr <- p$get_result()
    list(ok = is.finite(thr), thr = as.numeric(thr), timed_out = FALSE, err = NULL)
  }
  
  run_gmm_fork <- function(x, spc_target, timeout) {
    if (.Platform$OS.type != "unix") {
      stop("engine='fork' is only available on Unix/macOS.")
    }
    j <- parallel::mcparallel(gmm_worker(x, spc_target), mc.set.seed = FALSE)
    res <- parallel::mccollect(j, timeout = timeout)
    if (is.null(res) || is.null(res[[1]])) {
      # timed out: kill child
      try(suppressWarnings(tools::pskill(j$pid)), silent = TRUE)
      return(list(ok = FALSE, thr = NA_real_, timed_out = TRUE, err = "timeout"))
    }
    thr <- res[[1]]
    list(ok = is.finite(thr), thr = as.numeric(thr), timed_out = FALSE, err = NULL)
  }
  
  run_gmm_inline <- function(x, spc_target, timeout) {
    setTimeLimit(elapsed = timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    out <- tryCatch(gmm_worker(x, spc_target), error = identity)
    if (inherits(out, "error")) {
      to <- grepl("elapsed time limit", out$message, fixed = TRUE)
      return(list(ok = FALSE, thr = NA_real_, timed_out = to, err = out$message))
    }
    list(ok = is.finite(out), thr = as.numeric(out), timed_out = FALSE, err = NULL)
  }
  
  # Choose engine automatically (prefer callr background)
  if (engine == "auto") {
    engine <- if (requireNamespace("callr", quietly = TRUE)) "callr_bg"
    else if (.Platform$OS.type == "unix") "fork"
    else "inline"
  }
  
  msg("GMM attempt on n=%s via engine='%s' (timeout=%ss)", length(x), engine, gmm_timeout)
  
  gmm_res <- switch(
    engine,
    callr_bg = run_gmm_callr_bg(x, spc_target, gmm_timeout),
    fork     = run_gmm_fork(x, spc_target, gmm_timeout),
    inline   = run_gmm_inline(x, spc_target, gmm_timeout)
  )
  
  if (isTRUE(gmm_res$ok) && is.finite(gmm_res$thr)) {
    msg("GMM threshold: %.6f", gmm_res$thr)
    return(gmm_res$thr)
  } else {
    if (isTRUE(gmm_res$timed_out)) {
      msg("⚠️ GMM threshold estimation timed out and was killed.")
    } else if (!is.null(gmm_res$err)) {
      msg("⚠️ GMM threshold estimation failed: %s", gmm_res$err)
    } else {
      msg("⚠️ GMM threshold estimation failed for unknown reasons.")
    }
    msg("Falling back to density valley…")
  }
  
  # Density valley fallback
  d <- stats::density(x, bw = density_bw, n = density_n)
  y <- d$y; z <- d$x
  dy <- diff(y); s <- sign(dy); t <- diff(s)
  peaks_idx <- which(t == -2) + 1L
  pits_idx  <- which(t ==  2) + 1L
  
  if (length(peaks_idx) >= 2L && length(pits_idx) > 0L) {
    top2 <- peaks_idx[order(y[peaks_idx], decreasing = TRUE)][seq_len(2L)]
    rng  <- sort(range(top2))
    between <- pits_idx[pits_idx > rng[1] & pits_idx < rng[2]]
    if (length(between)) {
      thr_den <- z[between[which.min(y[between])]]
      if (is.finite(thr_den)) {
        msg("Density valley threshold: %.6f", thr_den)
        return(thr_den)
      }
    }
  }
  
  thr_q <- stats::quantile(x, 0.95, names = FALSE, type = 7)
  msg("Density valley unavailable; using conservative 95th quantile: %.6f", thr_q)
  thr_q
}




# ------------------------- IMGT references (configurable) -------------------------
read_imgt_references <- function(cfg) {
  # allow override via cfg; fall back to your original path
  imgt_dir <- cfg$paths$imgt_dir %||% "~/share/germlines/imgt/human/vdj"
  log_message("Reading IMGT references from: ", imgt_dir)
  readIMGT(dir = imgt_dir)
}

# ------------------------- TCR ingestion -------------------------
read_tcr <- function(vdj_t_dirs, sample_ids) {
  if (!length(vdj_t_dirs)) return(NULL)
  log_message("Loading TCR contigs for N runs: ", length(vdj_t_dirs))
  
  # Ensure sample_ids match directory basenames
  sample_map <- stringr::str_split(vdj_t_dirs, "/", simplify = TRUE)[,4]
  
  tcr_files <- file.path(vdj_t_dirs, "filtered_contig_annotations.csv")
  tcr_files <- tcr_files[file.exists(tcr_files)]
  if (!length(tcr_files)) return(NULL)
  
  tcr_contigs <- lapply(tcr_files, read.csv)
  # combineTCR expects a list of data.frames and a matching 'samples' vector
  combined <- scRepertoire::combineTCR(tcr_contigs, 
                                       samples = sample_map, 
                                       filterMulti = TRUE) |> bind_rows()
  
  # Normalize barcode to match Seurat cell names
  combined$barcode <- normalize_barcode(combined$barcode)
  
  # Minimal harmonized columns for downstream merge
  combined |>
    transmute(
      barcode = barcode,
      sample  = as.character(sample),
      CTaa    = CTaa,
      CTgene  = CTgene,
      CTstrict= CTstrict
    )
}

# ------------------------- BCR ingestion + cloning across runs -------------------------
bcr_pipeline <- function(vdj_b_dirs, subject, references, threads = 2L, engine = "auto") {
  if (!length(vdj_b_dirs)) return(NULL)
  log_message("Loading BCR contigs for N runs (subject ", subject, "): ", length(vdj_b_dirs))
  
  bcr_files <- file.path(vdj_b_dirs, "filtered_contig_igblast_db-pass.tsv")
  bcr_files <- bcr_files[file.exists(bcr_files)]
  if (!length(bcr_files)) return(NULL)
  
  bcr_list <- lapply(seq_along(bcr_files), function(i) {
    df <- read.delim(bcr_files[i], stringsAsFactors = FALSE)
    sample_id <- stringr::str_split(vdj_b_dirs[i], "/", simplify = TRUE)[,4]
    df |>
      mutate(
        sample_id = sample_id,
        subject_id = subject,
        cell_id_unique = paste0(sample_id, ".", .data$cell_id)
      ) |>
      filter(productive == TRUE)
  })
  
  # Remove cells with >1 heavy chain in each run
  bcr_list <- lapply(bcr_list, function(df) {
    multi_heavy_cells <- df |>
      filter(locus == "IGH") |>
      dplyr::count(cell_id, name = "n") |>
      filter(n > 1) |>
      pull(cell_id)
    df |> filter(!cell_id %in% multi_heavy_cells)
  })
  
  # Remove rows without isotype
  bcr_list <- lapply(bcr_list, \(df) df |> filter(c_call != ""))
  
  # Ensure each remaining cell has a heavy chain
  bcr_list <- lapply(bcr_list, function(df) {
    heavy_cells <- df |> filter(locus == "IGH") |> pull(cell_id)
    light_cells <- df |> filter(locus %in% c("IGK", "IGL")) |> pull(cell_id)
    no_heavy <- setdiff(light_cells, heavy_cells)
    df |> filter(!cell_id %in% no_heavy)
  })
  
  bcr_db <- bind_rows(bcr_list)
  
  if (!nrow(bcr_db)) return(NULL)
  
  # Distances / threshold / hierarchical clones
  dist_heavy <- distToNearest(
    bcr_db,
    cellIdColumn = "cell_id_unique",
    first = FALSE,
    onlyHeavy = TRUE,
    fields = "subject_id",
    nproc = threads,
    progress = TRUE
  )
  
  thr <- robust_find_threshold(
    as.numeric(dist_heavy$dist_nearest),
    spc_target = 0.995,
    gmm_timeout = 120,     
    max_n = 200000,
    top_trim = 0.999,
    verbose = TRUE, 
    engine = engine
  )
  
  results <- hierarchicalClones(
    dist_heavy,
    cell_id          = "cell_id_unique",
    threshold        = thr,
    only_heavy       = TRUE,
    split_light      = FALSE,
    summarize_clones = TRUE,
    nproc            = threads,
    verbose          = TRUE
  )
  results <- as.data.frame(results@db)
  
  
  results$sequence_id <- paste0(results$sample_id, ".", results$sequence_id)
  
  # Reconstruct germlines
  results <- createGermlines(results, reference = references)
  results$c_call[is.na(results$c_call)] <- "NA"

  # Save processed AIRR table for downstream lineage analysis
  airr_dir <- file.path("outputs", "tables", "bcr_airr")
  ensure_dir(airr_dir)
  write.table(results,
              file.path(airr_dir, paste0(subject, "_airr.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Resolve light chains into subclones
  subcl <- resolveLightChains(
    data = results,
    cell = "cell_id_unique",
    id   = "sequence_id",
    nproc = threads
  ) |>
    mutate(
      v_gene = str_split_i(v_call, "[*]", 1),
      j_gene = str_split_i(j_call, "[*]", 1),
      c_gene = str_split_i(c_call, "[*]", 1),
      chain  = ifelse(locus == "IGH", "heavy", "light")
    )
  
  # Somatic hypermutation summaries (counts + freq; combined + per-region)
  shm <- observedMutations(
    subcl,
    sequenceColumn = "sequence_alignment",
    germlineColumn = "germline_alignment_d_mask",
    regionDefinition = IMGT_VDJ,
    frequency = FALSE, combine = FALSE, nproc = threads
  ) |> 
    observedMutations(
      sequenceColumn = "sequence_alignment",
      germlineColumn = "germline_alignment_d_mask",
      regionDefinition = IMGT_VDJ,
      frequency = FALSE, combine = TRUE, nproc = threads
    ) |>
    observedMutations(
      sequenceColumn = "sequence_alignment",
      germlineColumn = "germline_alignment_d_mask",
      regionDefinition = IMGT_VDJ,
      frequency = TRUE, combine = FALSE, nproc = threads
    ) |>
    observedMutations(
      sequenceColumn = "sequence_alignment",
      germlineColumn = "germline_alignment_d_mask",
      regionDefinition = IMGT_VDJ,
      frequency = TRUE, combine = TRUE, nproc = threads
    )
  
  # Keep a single row per heavy+light pair via wide pivot
  cols_trim <- c("v_call", "d_call", "j_call", "c_call")
  shm[cols_trim] <- lapply(shm[cols_trim], \(x) gsub(",.*", "", x))
  
  shm <- shm |> arrange(chain)
  
  bcr_out <- shm |>
    mutate(d_call = replace_na(d_call, "NA")) |>
    pivot_wider(
      id_cols = c(cell_id_unique, subject_id, clone_id, sample_id),
      names_from = chain,
      values_from = c(junction_aa, v_call, d_call, j_call, c_call,
                      mu_count, mu_freq, mu_count_cdr_r, mu_count_cdr_s,
                      mu_count_fwr_r, mu_count_fwr_s, mu_freq_cdr_r, mu_freq_cdr_s,
                      mu_freq_fwr_r, mu_freq_fwr_s),
      values_fn = list
    ) |>
    unnest(cols = everything()) |>
    mutate(
      across(c(junction_aa_light, v_call_light, d_call_light, j_call_light, c_call_light),
             ~ replace_na(., "NA")),
      barcode = cell_id_unique,
      sample  = sample_id,
      CTaa    = paste(junction_aa_heavy, junction_aa_light, sep = "_"),
      heavy_gene_part = paste(v_call_heavy, d_call_heavy, j_call_heavy, c_call_heavy, sep = "."),
      light_gene_part = paste(v_call_light, j_call_light, c_call_light, sep = "."),
      CTgene  = paste(heavy_gene_part, light_gene_part, sep = "_"),
      CTstrict= paste(subject_id, clone_id, sep = "_")
    ) |>
    dplyr::select(
      barcode, sample, CTaa, CTgene, CTstrict,
      starts_with("mu_")
    )
  
  bcr_out
}

# ------------------------- Main driver -------------------------
annotate_AIRR <- function(cfg) {
  threads <- as.integer(cfg$threads %||% 2L)
  
  # 1) Load metadata
  meta <- read.csv("inputs/data/metadata.csv", row.names = 1, stringsAsFactors = FALSE)
  # Pull run ID from RNA_CRoutput robustly (take last path component)
  meta$run <- stringr::str_split(meta$RNA_CRoutput, "/", simplify = T)[,10]
  
  # 2) List candidate sample folders under parent_dir (skip dot files)
  samples <- list.files(cfg$paths$parent_dir, full.names = FALSE)
  samples <- samples[!grepl("^\\.", samples)]
  
  # 3) Load IMGT references
  references <- read_imgt_references(cfg)
  
  # 4) Loop by subject (BCR clones called across runs per subject)
  subjects <- unique(meta$Subject)
  
  for (subject in subjects) {
    log_message("---- Subject: ", subject, " ----")
    
    #problematic BCR samples
    if(subject %in% c("UV150","UV170","UV122","UV186","UV160","UV196","UV215")) {
      eng <- "callr_bg"
    } else {
      eng <- "fork"
    }
    
    meta_sub     <- meta |> filter(Subject == subject)
    sample_runs  <- intersect(samples, meta_sub$run) # runs actually present
    if (!length(sample_runs)) {
      log_message("No runs found on disk for subject: ", subject)
      next
    }
    
    # 4a) Load per-run Seurat objects (only those present)
    seurat_paths <- file.path(cfg$paths$processed_dir, paste0(sample_runs, ".rds"))
    exists_mask  <- file.exists(seurat_paths)
    if (!any(exists_mask)) {
      log_message("No Seurat RDS found for subject: ", subject)
      next
    }
    sample_runs   <- sample_runs[exists_mask]
    seurat_paths  <- seurat_paths[exists_mask]
    
    seurat_list <- lapply(seurat_paths, safe_readRDS)
    names(seurat_list) <- sample_runs
    
    # 4b) Attach study-level covariates from metadata into each Seurat object
    #     (vectorized assignment via AddMetaData; keeps object tidy)
    for (rn in sample_runs) {
      row <- meta_sub |> filter(run == rn) #|> slice(1)
      if (nrow(row) == 0) next
      
      md <- list(
        Subject          = row$Subject,
        Age              = row$Age,
        Disease_Activity = row$Disease_Activity,
        Phenotype        = row$Phenotype,
        Phenotype_2      = row$Phenotype_2,
        Phenotype_3      = row$Phenotype_3,
        Tissue_1         = row$Tissue_1,
        Tissue_2         = row$Tissue_2,
        Etiology         = row$Etiology
      )
      seurat_list[[rn]] <- AddMetaData(seurat_list[[rn]], md)
    }
    
    # 4c) Detect VDJ folders present per run
    vdj_t_dirs <- file.path(cfg$paths$parent_dir, sample_runs, "vdj_t")
    vdj_b_dirs <- file.path(cfg$paths$parent_dir, sample_runs, "vdj_b")
    vdj_t_dirs <- vdj_t_dirs[file.exists(vdj_t_dirs)]
    vdj_b_dirs <- vdj_b_dirs[file.exists(vdj_b_dirs)]
    
    # 4d) Run TCR and/or BCR pipelines conditionally
    combined_TCR <- read_tcr(vdj_t_dirs, sample_runs)
    combined_BCR <- bcr_pipeline(vdj_b_dirs, subject, references, threads = threads, engine = eng)
    
    if (is.null(combined_TCR) && is.null(combined_BCR)) {
      log_message("No TCR or BCR data available for subject: ", subject)
      # still save meta-updated Seurat objects
      for (rn in names(seurat_list)) {
        saveRDS(seurat_list[[rn]], file = file.path(cfg$paths$processed_dir, paste0(rn, ".rds")))
      }
      next
    }
    
    # 4e) Harmonize columns and combine (prefer BCR if duplicates on barcode)
    union_cols <- union(
      names(combined_TCR %||% tibble()),
      names(combined_BCR %||% tibble())
    )
    
    if (!is.null(combined_TCR)) combined_TCR <- fill_missing_cols(combined_TCR, union_cols)
    if (!is.null(combined_BCR)) combined_BCR <- fill_missing_cols(combined_BCR, union_cols)
    
    total_clones <- bind_rows(
      # put BCR first so it "wins" when deduplicating by barcode
      if (!is.null(combined_BCR)) combined_BCR,
      if (!is.null(combined_TCR)) combined_TCR
    ) |>
      filter(!is.na(barcode)) |>
      distinct(barcode, .keep_all = TRUE) |>
      as.data.frame()
    
    rownames(total_clones) <- total_clones$barcode
    
    # 4f) Add clone-level summaries back to each Seurat run and save
    for (rn in names(seurat_list)) {
      so <- seurat_list[[rn]]
      
      if(any(rownames(total_clones) %in% colnames(so))) {
        # Add per-cell clone annotations
        so <- AddMetaData(so, total_clones)
      } else {
        so[[setdiff(colnames(total_clones), c("barcode", "sample"))]] <- NA
      }
      
      
      # Compute clone summaries present in that run only (avoid cross-run bias)
      sc_meta <- so[[]] 
      
      if(any(colnames(sc_meta) %in% c("clonalFrequency", "clonalProportion", "cloneSize"))) {
        sc_meta <- sc_meta %>%
          dplyr::select(-c(clonalFrequency, clonalProportion, cloneSize))
      }
      
      # Class by CTgene prefix (IG* = BCR, TR* = TCR) — previously used
      # nchar(CTstrict) <= 10, which is fragile (locus-dependent CTstrict
      # lengths can cross the threshold).
      clones_in_run <- sc_meta |>
        filter(!is.na(CTstrict)) |>
        group_by(CTstrict) |>
        summarise(clonalFrequency = dplyr::n(),
                  CTgene_any      = dplyr::first(CTgene),
                  .groups = "drop") |>
        mutate(class = dplyr::case_when(
          grepl("^IG", CTgene_any) ~ "BCR",
          grepl("^TR", CTgene_any) ~ "TCR",
          TRUE                     ~ NA_character_
        )) |>
        dplyr::select(-CTgene_any)
      
      if (nrow(clones_in_run)) {
        binned <- compute_clone_bins(clones_in_run)
        
        clonal_join <- sc_meta |>
          left_join(binned, by = "CTstrict") |>
          dplyr::select(barcode, clonalProportion, clonalFrequency, cloneSize) |>
          filter(!is.na(barcode))
        
        rownames(clonal_join) <- clonal_join$barcode
        so <- AddMetaData(so, clonal_join)
      } else {
        # ensure columns exist (NA) if no clones present in this run
        so <- AddMetaData(so, tibble(
          barcode = Cells(so),
          clonalProportion = NA_real_,
          clonalFrequency  = NA_integer_,
          cloneSize        = NA_character_
        ) |> column_to_rownames("barcode"))
      }
      
      # Save updated Seurat object
      out_path <- file.path(cfg$paths$processed_dir, paste0(rn, ".rds"))
      saveRDS(so, file = out_path)
      log_message("Saved annotated Seurat: ", out_path)
    }
  }
  
  invisible(TRUE)
}
