#!/usr/bin/env Rscript
library(pryr)
library(Seurat)
library(Matrix)

# Command-line arguments
args <- commandArgs(trailingOnly=TRUE)
metric <- args[1]
n <- as.numeric(args[2])
outdir <- args[3]


# Data
file <- "data/panc8/panc8_sel.rds"
panc8 <- readRDS(file)
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

if (metric == "cms") {
  sce <- as.SingleCellExperiment(panc8_sub)
} else if (metric == "pvca") {
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
} else if (metric == "rvp") {
  X_mat <- GetAssayData(panc8_sub)
  metadata <- panc8_sub@meta.data
} else if (metric == "rvps") {
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
if (metric == "rvp") {
  source("R/rvp.R")
  start <- proc.time()
  obj <- rvp.default(X_mat, metadata$tech, metadata$celltype)
  duration <- proc.time() - start 
} else if (metric == "rvps") {
  source("R/rvp.R")
  start <- proc.time()
  obj <- rvp.sparseMatrix(X_mat, metadata$tech, metadata$celltype)
  duration <- proc.time() - start 
} else if (metric == "gpca") {
  source("R/gpca.R")
  start <- proc.time()
  # Modified to fix error when nperm = 0
  obj <- gPCA(X_mat, metadata$tech, nperm = 0) 
  duration <- proc.time() - start 
} else if (metric == "pvca") {
  library(pvca)
  start <- proc.time()
  obj <- pvcaBatchAssess(eset, c("tech", "celltype"), 0.6)
  duration <- proc.time() - start 
} else if (metric == "cms") {
  library(CellMixS)
  start <- proc.time()
  sce <- cms(sce, k = k, group = "tech")
  duration <- proc.time() - start 
} else if (metric == "kbet") {
  library(kBET)
  start <- proc.time()
  obj <- kBET(X_mat, metadata$tech, testSize = n, n_repeat = 1)
  duration <- proc.time() - start 
} else if (metric == "lisi") {
  library(lisi)
  start <- proc.time()
  obj <- compute_lisi(X_mat, metadata, c("tech"))
  duration <- proc.time() - start 
}
cat(sprintf("%s %d ", metric, n))
cat(duration, fill = TRUE)

# Save results
file <- sprintf("%s%s-%d.rds", outdir, metric, n)
message(sprintf("Saving results to: %s", file))
message("==========\n")
if (metric == "cms") {
  saveRDS(sce@colData, file)
} else {
  saveRDS(obj, file)
}

# check size of variable
# print(pryr::object_size(X_mat, eset))
# str(as.list(.GlobalEnv))
