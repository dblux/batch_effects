#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)

# Command-line arguments
args <- commandArgs(trailingOnly=TRUE)
metric <- args[1]
n <- as.numeric(args[2])
outdir <- args[3]
# design <- args[4]


file <- "data/panc8/panc8_sel.rds"
panc8 <- readRDS(file)

# # restrict to two batches and two classes
# if (design == "two") {
#   platforms <- c("indrop", "celseq2")
#   celltypes <- c("alpha", "beta")
#   panc8 <- panc8[, panc8$celltype %in% celltypes & panc8$tech %in% platforms]
# }

set.seed(1)
repeat {
  message(sprintf("Sampling n = %d samples from panc8.", n))
  sid <- sample(seq_len(ncol(panc8)), n)
  panc8_sub <- panc8[, sid]
  nbatches <- rowSums(table(panc8_sub$celltype, panc8_sub$tech) != 0)
  # re-sample if any class contains only one batch (zero batches is ok)
  if (all(nbatches != 1)) {
    break
  }
  set.seed(NULL)
}
rm(panc8)

# Load data
if (metric == "CMS") {
  sce <- as.SingleCellExperiment(panc8_sub)
} else if (metric == "PVCA") {
  library(Biobase)
  var_metadata <- data.frame(
    labelDescription = colnames(panc8_sub@meta.data),
    row.names = colnames(panc8_sub@meta.data)
  )
  pheno_data <- new(
    "AnnotatedDataFrame",
    data = panc8_sub@meta.data,
    varMetadata = var_metadata
  )
  eset <- Biobase::ExpressionSet(
    assayData = as.matrix(GetAssayData(panc8_sub)),
    phenoData = pheno_data 
  )
} else if (metric %in% c("HVP", "HVPS")) {
  X_mat <- GetAssayData(panc8_sub)
  metadata <- panc8_sub@meta.data
} else {
  X_mat <- t(GetAssayData(panc8_sub))
  metadata <- panc8_sub@meta.data
}
rm(panc8_sub)

# Benchmarking
k <- n / 10 
message(sprintf("Benchmarking: %s (n = %d)", metric, n))
if (metric == "HVP") {
  source("R/HVP.R")
  start <- proc.time()
  obj <- .HVP(X_mat, metadata$tech, metadata$celltype)
  duration <- proc.time() - start 
} else if (metric == "HVPS") {
  source("R/HVP.R")
  start <- proc.time()
  obj <- .HVP_sparseMatrix(X_mat, metadata$tech, metadata$celltype)
  duration <- proc.time() - start 
} else if (metric == "gPCA") {
  source("R/gpca.R")
  start <- proc.time()
  # gPCA has been modified to fix error when nperm = 0
  obj <- gPCA(X_mat, metadata$tech, nperm = 0) 
  duration <- proc.time() - start 
} else if (metric == "PVCA") {
  library(pvca)
  start <- proc.time()
  obj <- pvcaBatchAssess(eset, c("tech", "celltype"), 0.6)
  duration <- proc.time() - start 
} else if (metric == "CMS") {
  library(CellMixS)
  library(scater)
  library(BiocSingular)
  set.seed(1)
  start <- proc.time()
  sce <- cms(
    # scater: runPCA - defaults to ntop = 500
    scater::runPCA(
      sce, 
      ncomponents = 50, 
      ntop = nrow(sce), 
      BSPARAM = IrlbaParam()
    ),
    k = k, 
    group = "tech"
  )
  duration <- proc.time() - start
} else if (metric == "kBET") {
  library(kBET)
  set.seed(1)
  start <- proc.time()
  obj <- kBET(
    X_mat,
    metadata$tech, 
    k0 = k,
    testSize = n, 
    do.pca = TRUE,
    heuristic = FALSE,
    n_repeat = 1
  )
  duration <- proc.time() - start
} else if (metric == "LISI") {
  library(lisi)
  library(BiocSingular)
  set.seed(1)
  start <- proc.time()
  obj <- compute_lisi(
    BiocSingular::runPCA(
      X_mat,
      rank = 50, 
      BSPARAM = IrlbaParam()
    )$x,
    metadata,
    c("tech")
  )
  duration <- proc.time() - start 
}
cat(sprintf("%s %d ", metric, n))
cat(duration, fill = TRUE)

# Save results
file <- sprintf("%s%s-%d.rds", outdir, metric, n)
message(sprintf("Saving results to: %s", file))
message("==========\n")
if (metric == "CMS") {
  saveRDS(sce@colData, file)
} else {
  saveRDS(obj, file)
}

