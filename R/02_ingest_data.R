# R/02_ingest_data.R
suppressPackageStartupMessages({
  library(Azimuth)
  library(Seurat)
  library(Matrix)
  library(patchwork)
  library(ggplot2)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(celldex)
  library(SingleR)
})


HPCA   <- celldex::HumanPrimaryCellAtlasData()
Monaco <- celldex::MonacoImmuneData()

qc_add_mito_ribo <- function(obj,
                             mito_pattern = "^MT-",
                             ribo_pattern = "^RP[SL][-]?" # RPS/RPL, allow optional dash
) {
  obj[["mito.perc"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)
  obj[["ribo.perc"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
  obj
}

safe_run_azimuth <- function(obj) {
  res <- tryCatch({
    suppressWarnings(
      RunAzimuth(
        obj,
        reference = "pbmcref",
        verbose = FALSE,
        mapping.score.k = 50,
        k.weight = 30
      )
    )
  }, error = function(e) {
    message("⚠️ Azimuth mapping failed: ", e$message)
    return(NULL)
  })
  
  # If Azimuth succeeded
  if (!is.null(res)) return(res)
  
  # If it failed — add NA columns
  missing_cols <- c("predicted.celltype.l1.score", 
                    "predicted.celltype.l1", 
                    "predicted.celltype.l2.score",
                    "predicted.celltype.l2", 
                    "predicted.celltype.l3.score",
                    "predicted.celltype.l3")
  
  for (col in missing_cols) {
    if (!col %in% colnames(obj[[]])) {
      obj[[col]] <- NA
    }
  }
  
  return(obj)
}

ingest_data <- function(cfg) {
  options(future.globals.maxSize= 89128960000)
  ensure_dir(cfg$paths$processed_dir)
  ensure_dir(cfg$paths$qc_dir)
  
  # list sample folders
  samples <- list.files(cfg$paths$parent_dir, full.names = FALSE)
  samples <- samples[!grepl("^\\.", samples)]
  
  qc_summary_list <- list()

  for (sample in samples) {
    log_message("Ingesting:", sample)
    tenx_dir <- file.path(cfg$paths$parent_dir, sample, "sample_filtered_feature_bc_matrix")
    mat <- Read10X(tenx_dir)
    if (inherits(mat, "list")) mat <- mat[[1]]
    obj <- CreateSeuratObject(counts = mat, 
                              assay = "RNA", 
                              project = sample)
    
    # fast QC metrics using sparse ops
    log_message("Performing QC for", sample)
    obj$nCount_RNA   <- Matrix::colSums(obj@assays$RNA$counts)
    obj$nFeature_RNA <- Matrix::colSums(obj@assays$RNA$counts != 0)
    obj <- subset(obj, subset = nFeature_RNA > cfg$qc$min_features)
    cells_pre_filter <- ncol(obj)

    obj <- RenameCells(obj, new.names = paste0(sample, ".", colnames(obj)))
    obj$Cohort <- cfg$cohort
    
    obj <- qc_add_mito_ribo(obj)
    
    # plots
    p_vln <- VlnPlot(obj,
                     features = c("nCount_RNA", "nFeature_RNA", "mito.perc", "ribo.perc"),
                     pt.size = 0, cols = "steelblue", ncol = 4) &
      theme(axis.text.x = element_blank(), 
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())
    ggsave(file.path(cfg$paths$qc_dir, paste0(sample, "_vln.pdf")), p_vln, width=9, height=3)
    
    # dynamic upper bound for nCount_RNA
    lcounts <- log(obj$nCount_RNA + 1)
    cut <- round(exp(mean(lcounts) + cfg$qc$ncount_sd * sd(lcounts)))
    
    p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "mito.perc", cols = "steelblue") +
      geom_hline(yintercept = cfg$qc$max_mito, lty = 2) +
      geom_vline(xintercept = cut, lty = 2) + 
      guides(color = "none")
    p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "steelblue") +
      geom_vline(xintercept = cut, lty = 2) + 
      guides(color = "none")
    ggsave(file.path(cfg$paths$qc_dir, paste0(sample, "_scatter.pdf")), p1 + p2, width=6, height=3)
    
    obj <- subset(obj, subset = mito.perc < cfg$qc$max_mito & nCount_RNA < cut)
    cells_post_qc <- ncol(obj)

    log_message("Running Doublet Annotation for ", sample)
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce, verbose = FALSE)
    obj$doublet.class <- sce$scDblFinder.class
    obj$doublet.score <- sce$scDblFinder.score
    n_doublets <- sum(obj$doublet.class == "doublet")
    obj <- subset(obj, subset = doublet.class == "singlet")
    cells_post_doublet <- ncol(obj)
    
    log_message("Running SingleR Annotation for ", sample)
    sr1 <- SingleR(test = sce, ref = HPCA,   
                   labels = HPCA$label.fine,   
                   assay.type.test = 1)
    sr2 <- SingleR(test = sce, ref = Monaco, 
                   labels = Monaco$label.fine, 
                   assay.type.test = 1)
    obj$HPCA.label   <- sr1$labels
    obj$HPCA.pruned  <- sr1$pruned.labels
    obj$Monaco.label <- sr2$labels
    obj$Monaco.pruned<- sr2$pruned.labels
    rm(sce, sr1, sr2)
    
    message("Running Azimuth Annotation for ", sample)
    obj <- SCTransform(obj, verbose = FALSE)
    az <- safe_run_azimuth(obj)
    az_cols <- grep("predicted", colnames(az[[]]), value = TRUE)
    obj <- AddMetaData(obj, az[[az_cols]])
    
    qc_summary_list[[sample]] <- data.frame(
      sample            = sample,
      cells_pre_filter  = cells_pre_filter,
      cells_post_qc    = cells_post_qc,
      cells_post_doublet = cells_post_doublet,
      pct_doublets      = round(100 * n_doublets / cells_post_qc, 2),
      median_nFeature   = median(obj$nFeature_RNA),
      median_nCount     = median(obj$nCount_RNA),
      median_mito_perc  = round(median(obj$mito.perc), 2),
      stringsAsFactors  = FALSE
    )

    saveRDS(obj, file = file.path(cfg$paths$processed_dir, paste0(sample, ".rds")))
    rm(obj, mat, az); gc()
  }

  # ---- Write QC summary table with metadata ----
  qc_df <- do.call(rbind, qc_summary_list)
  rownames(qc_df) <- NULL

  meta <- read.csv(file.path(cfg$paths$parent_dir, "..", "metadata.csv"),
                   stringsAsFactors = FALSE)
  meta$run_folder <- basename(dirname(dirname(meta$RNA_CRoutput)))
  meta_sub <- meta[, c("run_folder", "Subject", "Tissue_1",
                        "Phenotype_2", "Etiology")]
  meta_sub <- meta_sub[!duplicated(meta_sub$run_folder), ]
  qc_df <- merge(qc_df, meta_sub,
                 by.x = "sample", by.y = "run_folder", all.x = TRUE)

  out_path <- file.path(cfg$paths$qc_dir, "qc_summary.csv")
  write.csv(qc_df, out_path, row.names = FALSE)
  log_message("QC summary written to ", out_path)

  invisible(TRUE)
}
