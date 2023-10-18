library(ggplot2)
library(sva)

library(Seurat)
library(SingleCellExperiment)
library(CellMixS)
library(kBET)
library(lisi)
library(pvca)

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}

# Batch correction
dsnames <- c("with-bal", "with-imbal", "without-bal", "without-imbal")
for (dsname in dsnames) {
  file <- sprintf("data/simulated/microarray/dea/phi30/%s.rds", dsname)
  print(file)
  dataset <- readRDS(file)

  # ComBat - Not modelling class covariate
  combat_data <- ComBat(
    dataset$X, dataset$metadata$batch, mod = NULL,
    par.prior = TRUE, prior.plots = FALSE
  )

  dataset$X <- combat_data
  file <- sprintf("data/simulated/microarray/dea/phi30/combat-%s.rds", dsname)
  print(file)
  saveRDS(dataset, file)

  # ComBat - Modelling class covariate
  mod <- model.matrix(~as.factor(class), data = dataset$metadata)
  combat_cov_data <- ComBat(
    dataset$X, dataset$metadata$batch, mod,
    par.prior = TRUE, ref.batch = 1
  )

  dataset$X <- combat_cov_data
  file <- sprintf("data/simulated/microarray/dea/phi30/combat_cov-%s.rds", dsname)
  print(file)
  saveRDS(dataset, file)

  # sva correction is performed by combining the surrogate variables with the
  # covariates in the design matrix, finding the coefficients of a linear model
  # wrt to all covariates and substracting the surrogate variable (multipled by
  # its coefficient) from the expression values
  mod0 <- model.matrix(~1, data = dataset$metadata)
  svars <- sva(dataset$X, mod, mod0, n.sv = 1)
  newV <- NULL
  obj <- fsva(dataset$X, mod, svars)
  sva_data <- obj$db

  dataset$X <- sva_data
  file <- sprintf("data/simulated/microarray/dea/phi30/sva-%s.rds", dsname)
  print(file)
  saveRDS(dataset, file)
}

# Evaluate batch effects
dsnames <- list.files("data/simulated/microarray/dea/phi30", full.names = TRUE)
for (dsname in dsnames) {
  print(dsname)
  dataset <- readRDS(dsname)
  results <- eval_batch(
    dataset, "batch", "class",
    k.cms = 30, k0 = 30, perplexity = 30
  )
  file <- paste0("tmp/phi30/results-", strsplit(dsname, "/")[[1]][6])
  print(file)
  saveRDS(results, file)
}

# Batch effect metrics
files <- list.files("tmp/phi30", "results-", full.names = TRUE)
for (file in files[seq(3, 16, 4)]) {
  cat(file, fill = T)
  results <- readRDS(file)

  cat(paste("k =", results$kbet$params$k0), fill = T)
  cat(paste("rvp:", results$rvp$percent.batch), fill = T)
  cat(paste("gpca:", results$gpca$delta), fill = T)
  cat(paste("pvca:", results$pvca$dat[3]), fill = T)
  cat(paste("cms:", mean(results$cms$cms)), fill = T)
  cat(paste("kbet:", results$kbet$summary$kBET.observed[1]), fill = T)
  cat(paste("blisi:", mean(results$lisi$batch)), fill = T)
  cat(fill = T)
}

# DEA
dsnames <- list.files("data/simulated/microarray/dea/phi30", full.names = TRUE)
for (dsname in dsnames) {
  print(dsname)
  dataset <- readRDS(dsname)
  class <- dataset$metadata$class
  X1 <- data.frame(dataset$X[, class == "A"])
  X2 <- data.frame(dataset$X[, class == "B"])
  print(dim(X1))
  print(dim(X2))
  pval <- calc_univariate(t.test, X1, X2)
  file <- paste0("tmp/phi30/pvalue-", strsplit(dsname, "/")[[1]][6])
  print(file)
  saveRDS(pval, file)
}

# ROC
pfiles <- list.files("tmp/phi30", "pvalue-", full.names = TRUE)
dsnames <- list.files("data/simulated/microarray/dea/phi30", full.names = TRUE)
zipnames <- Map(c, pfiles, dsnames)

for (obj in zipnames[seq(3, 16, 4)]) {
  print(obj)
  pval <- readRDS(obj[1])
  dataset <- readRDS(obj[2])
  ALPHA <- 0.05
  y <- data.frame(
    true = as.numeric(1:10000 %in% dataset$diff.features),
    proba = pval,
    pred = as.numeric(pval < ALPHA)
  )

  print(calc_recall(y$true, y$pred))
  print(calc_specificity(y$true, y$pred))
  # TODO: Check how many features are affected by both
  # roc_obj <- ggplot_roc(y, "true", "pred", return.auc = TRUE)
  # print(roc_obj$auc[1])
  # dsname <- substring_head(strsplit(obj[2], "/")[[1]][5], 4)
  # file <- paste0("tmp/phi30/fig/roc-", dsname, ".png")
  # ggsave(file, roc_obj$plot, width = 5, height = 5)
}
# Class effects are much larger than batch effects? No

# Simulation parameters
dsname <- dsnames[13]
dataset <- readRDS(dsname)
print(dataset$params)
