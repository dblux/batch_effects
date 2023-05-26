#!/usr/bin/env Rscript

source('../relapse_prediction/R/rvp.R')
source('../relapse_prediction/R/gpca.R')

# CMDLINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
metric <- args[1]
n <- as.numeric(args[2])
outdir <- args[3]

file <- sprintf('data/batchqc/sizes/batchqc-%d.rds', n)

# METADATA 
ncond <- n / 4
batch <- as.factor(rep(1:2, each = ncond * 2))
class <- rep(rep(LETTERS[1:2], each = ncond), 2)
metadata <- data.frame(batch, class)  # assign rownames below

if (metric == 'RVP') {
  X <- t(readRDS(file))
  
  start <- Sys.time()
  rvp <- RVP(X, batch, class)
  end <- Sys.time()
  delta <- difftime(end, start, units = 'secs')
  
  file <- sprintf('%s%s-%d.rds', outdir, metric, n)
  saveRDS(rvp, file)
} else if (metric == 'gPCA') {
  X <- t(readRDS(file))

  start <- Sys.time()
  gpca <- gPCA.batchdetect(X, batch, nperm = 0) # modified to fix error when nperm = 0
  end <- Sys.time()
  delta <- difftime(end, start, units = 'secs')
} else if (metric == 'PVCA') {
  library(pvca)

  X_mat <- as.matrix(readRDS(file))
  rownames(metadata) <- colnames(X_mat)

  meta_metadata <- data.frame(labelDescription = colnames(metadata))
  pheno_data <- new("AnnotatedDataFrame", data = metadata, varMetadata = meta_metadata)
  eset <- Biobase::ExpressionSet(assayData = X_mat, phenoData = pheno_data)
  rm(X_mat)
 
  start <- Sys.time()
  pvca_obj <- pvcaBatchAssess(eset, c('batch', 'class'), 0.6)
  end <- Sys.time()
  delta <- difftime(end, start, units = 'secs')
} else if (metric == 'kBET') {
  library(kBET)
  
  X <- t(readRDS(file))

  start <- Sys.time()
  kbet <- kBET(X, batch, testSize = n, n_repeat = 1)
  end <- Sys.time()
  delta <- difftime(end, start, units = 'secs')

  file <- sprintf('%s%s-%d.rds', outdir, metric, n)
  saveRDS(kbet, file)
} else if (metric == 'LISI') {
  library(lisi)

  X <- t(readRDS(file))

  start <- Sys.time()
  lisi <- compute_lisi(X, metadata, c('batch'))
  end <- Sys.time()
  delta <- difftime(end, start, units = 'secs')

  file <- sprintf('%s%s-%d.rds', outdir, metric, n)
  saveRDS(lisi, file)
}

cat(sprintf('%s,%d,%f,', metric, n, delta))

# check size of variable
# print(pryr::object_size(X_mat, eset))
# str(as.list(.GlobalEnv))
