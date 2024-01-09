library(magrittr)
library(splatter)
library(scater)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw(base_size = 7))
source("R/subset.R")
source("R/plot.R")


# Set parameters
# Simulate 6400 cells and subsample 4000 to induce batch-class imbalance
niter <- 5
batch_scales <- c(seq(0.01, 0.09, 0.01), seq(0.1, 0.3, 0.04))
for (batch_scale in batch_scales) {
  for (i in seq_len(niter)) {
    batch_scale <- 0.01
    params <- newSplatParams()
    params <- setParams(
      params,
      nGenes = 10000,
      batchCells = c(3200, 3200),
      group.prob = c(0.5, 0.5),
      batch.facLoc = 0, # log-normal
      batch.facScale = batch_scale # log-normal
    )
    cat(paste0("Seed: ",params@seed), fill = TRUE)
    splat <- splatSimulate(params, method = "groups")
    splat <- minimiseSCE(
      splat,
      rowData.keep = TRUE,
      colData.keep = TRUE,
      metadata.keep = TRUE,
      assays.keep = "counts",
      sparsify = "auto"
    )
    file <- sprintf(
      "data/simulated/scrna/loc0/splat-scale_%.2f-%02d.rds",
      batch_scale, i
    )
    print(file)
    saveRDS(splat, file)

    # Log transform data
    splat <- logNormCounts(splat)

    # Subset
    tabnames <- list(sort(unique(splat$Group)), unique(splat$Batch))
    ct_bal <- matrix(1000, 2, 2, dimnames = tabnames)
    ct_imbal <- matrix(
      c(500, 1500, 1500, 500),
      2, 2, dimnames = tabnames
    )
    idx_bal <- idx_imbal <- numeric()
    for (g in rownames(ct_bal)) {
      for (k in colnames(ct_bal)) {
        j_kg <- which(splat$Batch == k & splat$Group == g)
        # balanced data set
        n1_kg <- ct_bal[g, k]
        j1_kg <- sample(j_kg, n1_kg)
        idx_bal <- c(idx_bal, j1_kg)
        # imbalanced data set
        n2_kg <- ct_imbal[g, k]
        j2_kg <- sample(j_kg, n2_kg)
        idx_imbal <- c(idx_imbal, j2_kg)
      }
    }
    splat_bal <- splat[, idx_bal] %>%
      minimiseSCE(colData.keep = TRUE, assays.keep = "logcounts")
    splat_imbal <- splat[, idx_imbal] %>%
      minimiseSCE(colData.keep = TRUE, assays.keep = "logcounts")
    print(table(splat_bal$Group, splat_bal$Batch))
    print(table(splat_imbal$Group, splat_imbal$Batch))

    # Save subsets 
    splat_objs <- list(
      bal = splat_bal,
      imbal = splat_imbal
    )
    file <- sprintf(
      "data/simulated/scrna/loc0/datasets-%.2f-%02d.rds",
      batch_scale, i
    )
    print(file)
    saveRDS(splat_objs, file)
  }
}

# No batch effects
niter <- 5
for (i in seq_len(niter)) {
  params <- newSplatParams()
  params <- setParams(
    params,
    nGenes = 10000,
    batchCells = 6400,
    group.prob = c(0.5, 0.5),
    batch.facLoc = 0, # log-normal
    batch.facScale = 0.01 # log-normal
  )
  cat(paste0("Seed: ",params@seed), fill = TRUE)
  splat <- splatSimulate(params, method = "groups")
  splat <- minimiseSCE(
    splat,
    rowData.keep = TRUE,
    colData.keep = TRUE,
    metadata.keep = TRUE,
    assays.keep = "counts",
    sparsify = "auto"
  )
  file <- sprintf("data/simulated/scrna/loc0/splat-scale_0.00-%02d.rds", i)
  print(file)
  saveRDS(splat, file)

  # Log transform data
  splat <- logNormCounts(splat)

  # Subset
  tabnames <- list(levels(splat$Group), c("Batch1", "Batch2")) 
  ct_bal <- matrix(1000, 2, 2, dimnames = tabnames)
  ct_imbal <- matrix(
    c(500, 1500, 1500, 500),
    2, 2, dimnames = tabnames
  )
  # Subset
  idx_bal <- idx_imbal <- numeric()
  batch_bal <- batch_imbal <- rep(NA, ncol(splat))
  for (g in rownames(ct_bal)) {
    j1_g <- j2_g <- which(splat$Group == g)
    for (k in colnames(ct_bal)) {
      # balanced data
      n1_kg <- ct_bal[g, k]
      j1_kg <- sample(j1_g, n1_kg)
      j1_g <- setdiff(j1_g, j1_kg)
      idx_bal <- c(idx_bal, j1_kg)
      batch_bal[j1_kg] <- k
      # imbalanced data
      n2_kg <- ct_imbal[g, k]
      j2_kg <- sample(j2_g, n2_kg)
      j2_g <- setdiff(j2_g, j2_kg)
      idx_imbal <- c(idx_imbal, j2_kg)
      batch_imbal[j2_kg] <- k
    }
  }
  splat_bal <- splat[, idx_bal]
  splat_bal$Batch <- batch_bal[idx_bal]
  splat_imbal <- splat[, idx_imbal]
  splat_imbal$Batch <- batch_imbal[idx_imbal]

  splat_bal <- splat_bal %>%
    minimiseSCE(colData.keep = TRUE, assays.keep = "logcounts")
  splat_imbal <- splat_imbal %>%
    minimiseSCE(colData.keep = TRUE, assays.keep = "logcounts")

  print(table(splat_bal$Group, splat_bal$Batch))
  print(table(splat_imbal$Group, splat_imbal$Batch))

  # Save subsets 
  splat_objs <- list(
    bal = splat_bal,
    imbal = splat_imbal
  )
  file <- sprintf("data/simulated/scrna/loc0/datasets-0.00-%02d.rds", i)
  print(file)
  saveRDS(splat_objs, file)
}


# Plot: Log-normalised data
batch_cols <- brewer.pal(8, "Dark2")[2:4]
class_cols <- brewer.pal(8, "Dark2")[5:8]

show.legend <- FALSE
width <- 1.5
batch_scales <- c(seq(0.00, 0.09, 0.01), seq(0.1, 0.3, 0.04))
for (batch_scale in batch_scales) {
  if (batch_scale == 0.3) {
    show.legend <- TRUE
    width <- 1.9
  }
  file <- sprintf(
    "data/simulated/scrna/loc0/datasets-%.2f-01.rds", batch_scale
  )
  dataset <- readRDS(file)
  # # Balanced
  # sce <- runPCA(dataset$bal, scale = FALSE)
  # sce_pca <- reducedDim(sce, "PCA")
  # var_pc <- attr(sce_pca, "percentVar") / 100
  # scale_title <- sprintf("Batch scale = %.02f", batch_scale)
  # colData(sce)$Batch <- factor(colData(sce)$Batch)
  # levels(colData(sce)$Batch) <- c(1, 2)
  # levels(colData(sce)$Group) <- c(1, 2)
  # ax <- ggplot_pca(
  #   sce_pca, colData(sce), col = "Batch", pch = "Group",
  #   do.pca = FALSE, var_pc = var_pc,
  #   show.legend = show.legend, plot.axis = FALSE,
  #   cex = 0.8, alpha = 0.7
  # ) +
  #   labs(title = scale_title) +
  #   theme(
  #     title = element_text(size = 6),
  #     plot.title = element_text(hjust = 0.5),
  #     axis.title.x = element_text(size = 5),
  #     axis.title.y = element_text(size = 5),
  #     legend.key.size = unit(4, "mm")
  #   ) +
  #   scale_color_manual(values = batch_cols)
  # file <- sprintf("tmp/fig/splat/pca-splat_bal_%.2f-01.jpg", batch_scale)
  # ggsave(file, ax, width = width, height = 1.5)
  # print(file)
  # Imbalanced
  sce <- runPCA(dataset$imbal, scale = FALSE)
  sce_pca <- reducedDim(sce, "PCA")
  var_pc <- attr(sce_pca, "percentVar") / 100
  scale_title <- sprintf("Batch scale = %.02f", batch_scale)
  colData(sce)$Batch <- factor(colData(sce)$Batch)
  levels(colData(sce)$Batch) <- c(1, 2)
  levels(colData(sce)$Group) <- c(1, 2)
  ax <- ggplot_pca(
    sce_pca, colData(sce), col = "Batch", pch = "Group",
    do.pca = FALSE, var_pc = var_pc,
    show.legend = show.legend, plot.axis = FALSE,
    cex = 0.8, alpha = 0.7
  ) +
    labs(title = scale_title) +
    theme(
      title = element_text(size = 6),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 5),
      axis.title.y = element_text(size = 5),
      legend.key.size = unit(4, "mm")
    ) +
    scale_color_manual(values = batch_cols)
  file <- sprintf("tmp/fig/splat/pca-splat_imbal_%.2f-01.jpg", batch_scale)
  ggsave(file, ax, width = width, height = 1.5)
  print(file)
}

# Splat - Methodology:
# Simulate batch log-normal terms for each batch
# Simulate batch cell means (BatchCellMeans)
# Simulate DE genes
# Simulate different library sizes for each sample (BaseCellMeans)
# Enforce mean-variance trend (CellMeans)
# Sample counts from Poisson distribution (TrueCounts)
# Apply dropouts
# Data is normalised and log-transformed before measuring for batch effects
