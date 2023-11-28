library(biomaRt)
library(Seurat)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


# Cellbench
# - Counts data
# - QC metrics calculated using scPipe
# - Doublets identified using demuxlet

file <- "data/cellbench/cellbench-seurat.rds"
cellbench <- readRDS(file)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filter = "ensembl_gene_id",
  values = rownames(cellbench),
  mart = ensembl
)
hgnc_symbol <- mapping$hgnc_symbol[
  match(rownames(cellbench), mapping$ensembl_gene_id)]
mito_genes <- rownames(cellbench)[grep("^MT-", hgnc_symbol)]
ribo_genes <- rownames(cellbench)[grep("^RP[SL]", hgnc_symbol)]
cellbench$percent.mito <-
  colSums(GetAssayData(cellbench, slot = "counts")[mito_genes, ]) /
  cellbench$nCount_originalexp
cellbench$percent.ribo <-
  colSums(GetAssayData(cellbench, slot = "counts")[ribo_genes, ]) /
  cellbench$nCount_originalexp

# FeatureScatter(
#   cellbench_sub,
#   feature1 = "percent.mito", feature2 = "percent.ribo",
#   group.by = "batch"
# )
# FeatureScatter(
#   cellbench_sub,
#   feature1 = "nCount_originalexp", feature2 = "nFeature_originalexp",
#   group.by = "batch"
# )

cellbench_sub <- subset(
  cellbench,
  subset = demuxlet_cls == "SNG" &
    nCount_originalexp < 2e5 &
    percent.mito < 0.2 &
    percent.ribo < 0.35
)
table(cellbench_sub$celltype, cellbench_sub$batch)

# Normalise data
cellbench_sub <- NormalizeData(cellbench_sub)

# Feature selection
# - Remove mito, ribo and sparse genes
sparse_genes <- remove_sparse(
  GetAssayData(cellbench_sub, "counts"),
  0.95, cellbench_sub$celltype,
  ret.features = TRUE
)
rm_genes <- unique(c(mito_genes, ribo_genes, sparse_genes))
cellbench_sel <- cellbench_sub[!(rownames(cellbench_sub) %in% rm_genes), ]

in_all <- cellbench_sub@assays$originalexp@meta.features$in_all
is_mito_ribo <- rownames(cellbench_sub) %in% c(mito_genes, ribo_genes)
cellbench_fltr <- cellbench_sub[in_all & !is_mito_ribo, ]
cellbench_fltr <- NormalizeData(cellbench_fltr)

counts_fltr <- GetAssayData(cellbench_fltr, "counts")
# counts <- GetAssayData(cellbench_sub, "counts")
pct_zero <- rowSums(counts_fltr == 0) / ncol(counts_fltr)
nonzero_mean <- apply(counts_fltr, 1, function(x) mean(x[x > 0]))
length(pct_zero)

# Marginals
llim <- 5
rlim <- llim + 0.5
file <- sprintf("tmp/fig/cellbench-marginal_%.1f.png", llim)
pct_zero_marginal <- pct_zero[which(nonzero_mean > llim & nonzero_mean < rlim)]
png(file)
hist(
  pct_zero_marginal,
  xlab = "Percentage of zeros",
  main = sprintf(
    "Features with mean in (%.1f, %.1f); n = %d",
    llim, rlim, length(pct_zero_marginal)
  )
)
dev.off()

# Ridge plots
i <- i + 1
table(counts_fltr[i, ])
i <- 2121
ax1 <- RidgePlot(cellbench_fltr, features = rownames(cellbench_fltr)[i])
ax2 <- RidgePlot(
  cellbench_fltr,
  features = rownames(cellbench_fltr)[i],
  group.by = "celltype"
)
ax <- plot_grid(ax1, ax2)
ax
file <- sprintf("tmp/fig/cellbench-feature%d.png", i)
ggsave(file, ax, width = 12, height = 5)

ax <- VlnPlot(cellbench_fltr, features = "nCount_originalexp") 
file <- "tmp/fig/cellbench-n_count.png"
ggsave(file, ax)


file <- "tmp/fig/cellbench_fltr-dropout_mean.png"
png(file)
plot(nonzero_mean, pct_zero)
dev.off()

dropout <- na.omit(data.frame(nonzero_mean, pct_zero))
dropout_fltr <- subset(dropout, subset = nonzero_mean < 30)
ax <- ggplot(dropout_fltr, aes(x = nonzero_mean, y = pct_zero)) + geom_hex()
file <- "tmp/fig/cellbench_hex-dropout_mean.png"
ggsave(file, ax)

dropseq <- counts_fltr[, startsWith(colnames(counts_fltr), "Dropseq")]
celseq <- counts_fltr[, startsWith(colnames(counts_fltr), "CELseq")]
tenx <- counts_fltr[, startsWith(colnames(counts_fltr), "10x")]

ax <- VlnPlot(
  cellbench_fltr,
  feature = "nFeature_originalexp",
  group.by = "batch"
)
file <- "tmp/fig/cellbench_fltr-n_feature.png"
ggsave(file, ax, height = 5, width = 5)

pct_zero <- rowSums(dropseq == 0) / ncol(dropseq)
nonzero_mean <- apply(dropseq, 1, function(x) mean(x[x > 0]))
dropout_dropseq <- na.omit(data.frame(
  nonzero_mean, pct_zero, batch = "Dropseq"
))
pct_zero <- rowSums(celseq == 0) / ncol(celseq)
nonzero_mean <- apply(celseq, 1, function(x) mean(x[x > 0]))
dropout_celseq <- na.omit(data.frame(
  nonzero_mean, pct_zero, batch = "CELseq"
))
pct_zero <- rowSums(tenx == 0) / ncol(tenx)
nonzero_mean <- apply(tenx, 1, function(x) mean(x[x > 0]))
dropout_tenx <- na.omit(data.frame(
  nonzero_mean, pct_zero, batch = "10x"
))

dropout_all <- rbind(dropout_tenx, dropout_dropseq, dropout_celseq)
dropout_all_fltr <- subset(dropout_all, subset = nonzero_mean < 20)
idx <- sample(1:nrow(dropout_all_fltr), 6000)
ax <- ggplot(
  dropout_all_fltr[idx, ],
  aes(x = nonzero_mean, y = pct_zero, col = batch)
) + geom_point(pch = 1, alpha = 1)

file <- "tmp/fig/cellbench_batch-dropout_mean.png"
ggsave(file, ax)



png("tmp/fig/maqc_a-dropout.png")
hist(pct_zero)
dev.off()


# Plot
# - Difference between log and non-log PCA plots

for (idx in names(cellbench_objs)) {
  cellbench_sel1 <- cellbench_objs[[idx]] %>%
    ScaleData(do.scale = FALSE) %>%
    RunPCA(features = rownames(.))
  ax1 <- DimPlot(
    cellbench_sel1, reduction = "pca",
    group.by = "batch", shuffle = TRUE
  )
  ax2 <- DimPlot(
    cellbench_sel1, reduction = "pca",
    group.by = "celltype", shuffle = TRUE
  )
  ax <- plot_grid(ax1, ax2, rel_widths = c(1, 1))
  file <- sprintf("tmp/cellbench-pca_%s.png", idx)
  ggsave(file, ax, width = 10, height = 5)
}

# plot(1 - percent.mito, cellbench$non_mt_percent, xlim = c(0, 1), ylim = c(0, 1))
# abline(a = 0, b = 1)
# plot(1 - percent.ribo, cellbench$non_ribo_percent, xlim = c(0, 1), ylim = c(0, 1))
# abline(a = 0, b = 1)

cellbench_fltr1 <- cellbench_fltr %>%
  NormalizeData() %>%
  ScaleData(do.scale = FALSE) %>%
  RunPCA(features = rownames(.))
ax1 <- DimPlot(
  cellbench_fltr1, reduction = "pca",
  group.by = "batch", shuffle = TRUE
)
ax2 <- DimPlot(
  cellbench_fltr1, reduction = "pca",
  group.by = "celltype", shuffle = TRUE
)
ax <- plot_grid(ax1, ax2, rel_widths = c(1, 1))
ax
file <- "tmp/fig/cellbench-pca_log_fvf.png"
ggsave(file, ax, width = 10, height = 5)

