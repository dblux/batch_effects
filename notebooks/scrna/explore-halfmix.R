library(Seurat)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

src_files <- list.files("R", full.names = TRUE)
for (f in src_files) {
  source(f)
  cat(sprintf('Sourced file: %s\n', f))
}


mean_nonzero <- function(X) {
  mu <- rowSums(X) / rowSums(X != 0)
  mu[is.nan(mu)] <- 0 # nan due to division by zero
  return(mu)
}


halfmix <- readRDS("data/jurkat_293t/processed/halfmix.rds")
counts <- GetAssayData(halfmix, "counts")
table(counts[,1])

# QC
ax <- VlnPlot(
  halfmix,
  features = c("nFeature_RNA"),
  ncol = 1, group.by = "celltype",
  pt.size = .2
)
ax
# file <- "tmp/halfmix_celltype-nfeature.png"
# ggsave(file, ax, width = 6, height = 6)

# # Findings
# - Extreme variability in total counts and number of missing values even in cells in the same batch and same class
# - Extremely high number of dropouts in a sample (90%)
# - Different celltypes have different total counts 
# - Not much difference between batches (different times)
# - <1% of features have counts above 10

FeatureScatter(
  subset(halfmix, orig.ident == "zheng"),
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
  group.by = "celltype"
)

# # Dropouts
# - A third of features have no values (~12000 samples)
# - Around 2000 features have <50% dropouts
# - Dropouts do not seem to be associated with batch
# - Probability of dropouts seem to be associated with counts

# hist(halfmix$nFeature_RNA / nrow(halfmix), breaks = 20, probability = TRUE)

# Plot PCA
# Check log
# Filter cells and features
# Seurat batch correction

# Subset mixed batch and split by celltype
zheng <- subset(halfmix, orig.ident == "zheng")
list_cells <- SplitObject(zheng, "celltype")

for (i in seq_len(length(list_cells))) {
  counts <- GetAssayData(list_cells[[i]], "counts")
  metafeatures <- data.frame(
    pct_dropout = rowSums(counts == 0) / ncol(counts),
    nzmean = mean_nonzero(counts),
    row.names = rownames(counts)
  )
  list_cells[[i]]$RNA@meta.features <- metafeatures
}


ax <- ggplot() +
  geom_point(
    data = list_cells[[1]]$RNA@meta.features,
    aes(x = nzmean, y = pct_dropout),
    color = "blue", pch = 1
  ) +
  geom_point(
    data = list_cells[[2]]$RNA@meta.features,
    aes(x = nzmean, y = pct_dropout),
    color = "red", pch = 1
  ) +
  xlim(0, 10)
file <- "tmp/zheng-pct_dropout-nzmean.png"
ggsave(file, ax, width = 5, height = 5)


counts <- GetAssayData(halfmix, "counts")
metafeatures <- data.frame(
  pct_dropout = rowSums(counts == 0) / ncol(counts),
  nzmean = mean_nonzero(counts),
  row.names = rownames(counts)
)
threshold <- 10
ax <- ggplot(metafeatures) +
  geom_point(
    aes(x = nzmean, y = pct_dropout),
    pch = 1
  ) +
  xlim(0, 10)
file <- "tmp/halfmix-pct_dropout-nzmean.pdf"
ggsave(file, ax, width = 5, height = 5)


x <- counts[, 1]
hist(x[x > 0 & x < 150])
length(x[x == 8])
# head(unname(sort(x, decreasing = TRUE)), 100)

# Plot
## PCA
halfmix <- NormalizeData(halfmix)
halfmix <- FindVariableFeatures(halfmix)
halfmix <- ScaleData(halfmix, do.scale = FALSE)
halfmix <- RunPCA(halfmix)
saveRDS(halfmix, "data/293t_jurkat/processed/halfmix-pca.rds")

ax <- PCAPlot(halfmix, shuffle = TRUE)
file <- "tmp/pca-halfmix.png"
ggsave(file, ax, width = 5, height = 5)

halfmix <- RunUMAP(halfmix, dims = 1:20)
ax <- DimPlot(halfmix, reduction = "umap")
file <- "tmp/umap-halfmix.png"
ggsave(file, ax, width = 5, height = 5)
