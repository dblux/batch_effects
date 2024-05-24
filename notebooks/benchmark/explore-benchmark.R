library(pryr)
library(Seurat)
library(Matrix)

# PCA -> kNN -> Compute metric
# kbET: svd - FNN::get.knn
# kBET: Performs heuristics (find optimal k) and repeats to calculate statistics
# LISI: No PCA - RANN::nn2
# CMS: scater::runPCA - BiocNeighbors::findKNN (defaults to exact k-NN)


metric <- "RVP"
n <- 4000
outdir <- NULL

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
} else if (metric == "RVP") {
  X_mat <- GetAssayData(panc8_sub)
  metadata <- panc8_sub@meta.data
} else if (metric == "RVPS") {
  X_mat <- GetAssayData(panc8_sub)
  metadata <- panc8_sub@meta.data
} else {
  X_mat <- t(GetAssayData(panc8_sub))
  metadata <- panc8_sub@meta.data
}

# Benchmarking
k <- n / 10 
message(sprintf("Benchmarking: %s (n = %d)", metric, n))
if (metric == "RVP") {
  source("R/rvp.R")
  start <- proc.time()
  obj <- rvp.default(X_mat, metadata$tech, metadata$celltype)
  duration <- proc.time() - start 
} else if (metric == "RVPS") {
  source("R/rvp.R")
  start <- proc.time()
  obj <- rvp.sparseMatrix(X_mat, metadata$tech, metadata$celltype)
  duration <- proc.time() - start 
} else if (metric == "CMS") {
  library(CellMixS)
  library(scater)
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
  library(BiocSingular)
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
  set.seed(2)
  start <- proc.time()
  obj <- compute_lisi(
    # BiocSingular: runPCA - PCA using approximate algorithms
    BiocSingular::runPCA(
      X_mat,
      rank = 50, 
      BSPARAM = IrlbaParam()
    )$x,
    metadata, c("tech")
  )
  duration <- proc.time() - start 
}
cat(sprintf("%s %d ", metric, n))
cat(duration, fill = TRUE)


print(mean(obj$tech)) # LISI 

Rprof(memory.profiling = TRUE, interval=.002)
obj <- RVP(X_mat, metadata$tech, metadata$celltype, use.sparse = TRUE)
Rprof(NULL)
summaryRprof(memory = 'both')
