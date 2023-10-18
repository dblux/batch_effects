#!/usr/bin/env Rscript

source('../relapse_prediction/R/rvp.R')
source('../relapse_prediction/R/gpca.R')

# CMDLINE ARGUMENTS
# outdir: directory to save metric objects 
args = commandArgs(trailingOnly=TRUE)
metric <- args[1]
batch_delta <- as.numeric(args[2])
outdir <- args[3]

file <- sprintf('data/batchqc/large/imbalanced/imbal_large-%d.rds', batch_delta)

# METADATA 
n <- 6000
ncond <- n / 4
batch <- as.factor(rep(1:2, each = ncond * 2))
class <- rep(rep(LETTERS[1:2], each = ncond), 2)
metadata <- data.frame(batch, class)  # assign rownames below
# idx <- seq(1001, 5000) # imbalanced
idx <- c(501:2500, 3501:5500) # balanced
metadata <- metadata[idx, ]


if (metric == 'RVP') {
  X <- t(readRDS(file))[idx, ]
  rvp <- RVP(X, metadata$batch, metadata$class)
  file <- sprintf('%s%s-%d.rds', outdir, metric, batch_delta)
  saveRDS(rvp, file)
} else if (metric == 'gPCA') {
  X <- t(readRDS(file))[idx, ]
  gpca <- gPCA.batchdetect(X, metadata$batch, nperm = 0) # modified to fix error when nperm = 0
  file <- sprintf('%s%s-%d.rds', outdir, metric, batch_delta)
  saveRDS(gpca, file)
} else if (metric == 'PVCA') {
  library(pvca)

  X_mat <- as.matrix(readRDS(file))[, idx]
  rownames(metadata) <- colnames(X_mat)
  meta_metadata <- data.frame(labelDescription = colnames(metadata))
  pheno_data <- new("AnnotatedDataFrame", data = metadata, varMetadata = meta_metadata)
  eset <- Biobase::ExpressionSet(assayData = X_mat, phenoData = pheno_data)
  
  pvca_obj <- pvcaBatchAssess(eset, c('batch', 'class'), 0.6)
  file <- sprintf('%s%s-%d.rds', outdir, metric, batch_delta)
  saveRDS(pvca_obj, file)

} else if (metric == 'kBET') {
  library(kBET)
  
  X <- t(readRDS(file))[idx, ]
  kbet <- kBET(X, metadata$batch, testSize = n, n_repeat = 1, k0 = 1000)
  file <- sprintf('%s%s_k1000-%d.rds', outdir, metric, batch_delta)
  saveRDS(kbet, file)
} else if (metric == 'LISI') {
  library(lisi)

  X <- t(readRDS(file))[idx, ]
  lisi <- compute_lisi(X, metadata, c('batch'))
  file <- sprintf('%s%s-%d.rds', outdir, metric, batch_delta)
  saveRDS(lisi, file)
}

message(file)
