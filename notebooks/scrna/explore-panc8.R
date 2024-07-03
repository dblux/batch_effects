library(tidyr)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())
library(Seurat)
library(SeuratData)
src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  cat(f, fill = TRUE)
  source(f)
}


panc8 <- LoadData("panc8")

# Quality control - Cells
mito_genes <- rownames(panc8)[grep("^MT-", rownames(panc8))]
ribo_genes <- rownames(panc8)[grep("^RP[SL]", rownames(panc8))]
panc8$percent.mito <-
  colSums(GetAssayData(panc8, slot = "counts")[mito_genes, ]) /
  panc8$nCount_RNA
panc8$percent.ribo <-
  colSums(GetAssayData(panc8, slot = "counts")[ribo_genes, ]) /
  panc8$nCount_RNA

VlnPlot(panc8, features = c("nFeature_RNA", "nCount_RNA"))
VlnPlot(panc8, features = c("percent.mito", "percent.ribo"))
FeatureScatter(panc8_sub, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
FeatureScatter(panc8, feature1 = "percent.mito", feature2 = "percent.ribo")

# panc8_sub <- subset(panc8, subset = nCount_RNA < 2.5e6 & percent.mito < 0.2)
panc8 <- NormalizeData(panc8)

# Filter features with many dropouts
sparse_genes <- remove_sparse(
  GetAssayData(panc8, "counts"),
  0.95, panc8$celltype,
  ret.features = TRUE
)
panc8_sel <- panc8[!(rownames(panc8) %in% sparse_genes), ]
file <- "data/panc8/panc8_sel.rds"
saveRDS(panc8_sel, file)

# Load processed data
file <- "data/panc8/panc8_sel.rds"
panc8_sel <- readRDS(file)

batches <- c("indrop3", "smartseq2", "celseq2")
panc8_3batch <- panc8[, panc8$dataset %in% batches]

y <- panc8_3batch[[]]
table(panc8_3batch$dataset, panc8_3batch$celltype)

# Subset data
# - Data sets: indrop3, smartseq2
# - Cell types: alpha, ductal
## Expt: Degrees of imbalance
# - 2:2, 2:3, 2:4, 2:5, 2:6
# - Each batch: 480 cells
metadata <- panc2fltr@meta.data[c("celltype", "dataset")]
smartseq2_alpha <- which(metadata[["celltype"]] == "alpha" & metadata[["dataset"]] == "smartseq2")
smartseq2_ductal <- which(metadata[["celltype"]] == "ductal" & metadata[["dataset"]] == "smartseq2")
indrop3_alpha <- which(metadata[["celltype"]] == "alpha" & metadata[["dataset"]] == "indrop3")
indrop3_ductal <- which(metadata[["celltype"]] == "ductal" & metadata[["dataset"]] == "indrop3")

nperbatch <- 480
rng <- seq(2, 6)
list_imbal <- list()

for (i in rng) {
  frac <- 2 / (i + 2)
  n1 <- round(frac * nperbatch)
  n2 <- round((1 - frac) * nperbatch)
  idx1 <- sample(smartseq2_alpha, n1)
  idx2 <- sample(smartseq2_ductal, n2)
  idx3 <- sample(indrop3_alpha, n2)
  idx4 <- sample(indrop3_ductal, n1)
  idx <- sort(c(idx1, idx2, idx3, idx4))
  list_imbal <- c(list_imbal, panc2fltr[, idx])
}

names(list_imbal) <- rng

obj <- list_imbal[[5]]

batch_effects <- function(obj, batch, cls) {
  X <- LayerData(obj)
  metadata <- obj@meta.data
  kbet_size <- dim(obj)[2]
  
  RVP <- rvp(obj, batch, cls)
  kbet_estimate <- kBET(
    Matrix::t(X),
    obj[[batch]],
    testSize = kbet_size, n_repeat = 1
  )
  rejection_rate <- kbet_estimate$summary$kBET.observed[1]
  lisi_results <- compute_lisi(
    Matrix::t(X),
    metadata,
    c(batch, cls)
  )
  blisi <- mean(lisi_results[[batch]])
  
  return(c(rvp = RVP, kbet = rejection_rate, lisi = blisi))
}

results <- lapply(list_imbal, batch_effects)

results1 <- results %>%
  data.frame()

results <- readRDS("../tmp/results.rds")
results1 <- data.frame(lapply(results, unlist))
colnames(results1) <- substring(colnames(results1), 2, 2)
results1[["metric"]] <- rownames(results1)
results2 <- gather(results1, "imbalance", "value", -c("metric"))

ggplot(results2) +
  geom_line(aes(x = imbalance, y = value, color = metric, group = metric))

# panc8_flt <- NormalizeData(panc8_flt)
# panc8_flt <- FindVariableFeatures(panc8_flt, selection.method = "vst", nfeatures = 2000)

# Plot: RVP

SS <- rvp$sum_squares
# sort by ss_total
SS_t <- SS[order(-SS$ss_total), ]
# sort by ss_batch
SS_b <- SS[order(-SS$ss_batch), ]
# sort by mean
counts <- GetAssayData(panc8_fltr)
SS_m <- SS[order(-rowMeans(counts)), ]
# sort by mean in top 10% cells
mean_counts_10 <- apply(counts, 1, function(x) {
  pct <- 0.1
  mean(sort(x, decreasing = TRUE)[seq_len(pct * ncol(counts))])
})
SS_m10 <- SS[order(-mean_counts_10), ]

ax_rvp <- plot.rvp(SS_m10, cex = .5)
ax_rvp

file <- "~/Dropbox/tmp/ssm10-panc8.pdf"
ggsave(file, ax_rvp, width = 6, height = 3)

panc8_fltr <- ScaleData(panc8_fltr)
panc8_fltr <- RunPCA(panc8_fltr, features = rownames(panc8_fltr), npcs = 20)
ax_pca <- DimPlot(panc8_fltr) +
  labs(subtitle = sprintf("RVP = %.4f", rvp$percentage))
ax_pca
file <- '~/Dropbox/tmp/pca-panc8_fltr.pdf'
ggsave(file, ax_pca, width = 6, height = 4)

panc8_fltr <- RunUMAP(panc8_fltr, dims = 1:20)
head(panc8_fltr@meta.data)

ax_umap <- DimPlot(panc8_fltr, reduction = "umap", group.by = "celltype")
ax_umap
file <- '~/Dropbox/tmp/umap_celltype-panc8_fltr.pdf'
ggsave(file, ax_umap, width = 6, height = 4)

# Explore SeuratData
subset(AvailableData(), species == "mouse")

# Benchmarking
file <- "data/panc8/panc8_sel.rds"
panc8 <- readRDS(file)
head(panc8@meta.data)
table(panc8$celltype, panc8$tech)
# no embeddings in panc8
str(panc8)
