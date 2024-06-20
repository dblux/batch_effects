library(Seurat)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw(base_size = 7))
library(biomaRt)
library(kBET)
library(lisi)
src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}

# Dataset: Cellbench
# - Counts data
# - QC metrics calculated using scPipe

file <- "data/cellbench/cellbench-seurat.rds"
cellbench <- readRDS(file)

# QC
# - Doublets identified using demuxlet
# workaround for SSL cert error
httr::set_config(httr::config(ssl_verifypeer = FALSE))
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

# Feature selection: Remove mito, ribo and sparse genes
sparse_genes <- remove_sparse(
  GetAssayData(cellbench_sub, "counts"),
  0.95, cellbench_sub$celltype,
  ret.features = TRUE
)
rm_genes <- unique(c(mito_genes, ribo_genes, sparse_genes))
cellbench_sel <- cellbench_sub[!(rownames(cellbench_sub) %in% rm_genes), ]
# in_all <- cellbench_sel@assays$originalexp@meta.features$in_all
# table(in_all)

# Subset
# with batch effects
tabnames <- list(
  unique(cellbench_sel$celltype),
  unique(cellbench_sel$batch)
)
ct_bal <- matrix(40, 3, 3, dimnames = tabnames)
ct_imbal <- matrix(
  c(20, 60, 40, 40, 20, 60, 60, 40, 20),
  3, 3, dimnames = tabnames
)
# with batch effects
idx_bal <- idx_imbal <- numeric()
for (g in rownames(ct_bal)) {
  for (k in colnames(ct_bal)) {
    j_kg <- which(cellbench_sel$batch == k & cellbench_sel$celltype == g)
    if (length(j_kg) == 0)
      next
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
cellbench_bal <- cellbench_sel[, idx_bal]
cellbench_imbal <- cellbench_sel[, idx_imbal]
table(cellbench_bal$celltype, cellbench_bal$batch)
table(cellbench_imbal$celltype, cellbench_imbal$batch)


# without batch effects
cellbench_b1 <- subset(cellbench_sel, subset = batch == "10x")
idx_bal <- idx_imbal <- numeric()
batch_bal <- batch_imbal <- rep(NA, ncol(cellbench_b1))
for (g in rownames(ct_bal)) {
  j1_g <- j2_g <- which(cellbench_b1$celltype == g)
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
cellbench_b1_bal <- cellbench_b1[, idx_bal]
cellbench_b1_bal$batch <- batch_bal[idx_bal]
cellbench_b1_imbal <- cellbench_b1[, idx_imbal]
cellbench_b1_imbal$batch <- batch_imbal[idx_imbal]
table(cellbench_b1_bal$celltype, cellbench_b1_bal$batch)
table(cellbench_b1_imbal$celltype, cellbench_b1_imbal$batch)

# Evaluate batch effects

cellbench_objs <- list(
  negctrl_bal = cellbench_b1_bal,
  negctrl_imbal = cellbench_b1_imbal,
  bal = cellbench_bal,
  imbal = cellbench_imbal
)
saveRDS(cellbench_objs, "tmp/cellbench-datasets.rds")


file <- "tmp/cellbench-k60.rds"
objs <- readRDS(file)
for (idx in names(objs)) {
  cat(idx, fill = T)
  obj <- objs[[idx]]
  cat(obj$percent.batch, fill = T)
  cat(obj$p.value, fill = T)
}
for (idx in names(objs)) {
  cat(idx, fill = T)
  obj <- objs[[idx]]
  cat(obj$kbet$summary$kBET.observed[1], fill = T)
  cat(obj$kbet$params$k0, fill = T)
  cat(paste("blisi:", mean(obj$lisi$batch)), fill = T)
  cat(fill = T)
}

# Plot
- Difference between log and non-log PCA plots

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

cellbench_sub1 <- cellbench_sub %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(do.scale = FALSE) %>%
  RunPCA()
ax1 <- DimPlot(
  cellbench_sub1, reduction = "pca",
  group.by = "batch", shuffle = TRUE
)
ax2 <- DimPlot(
  cellbench_sub1, reduction = "pca",
  group.by = "celltype", shuffle = TRUE
)
ax <- plot_grid(ax1, ax2, rel_widths = c(1, 1))
file <- "tmp/fig/cellbench-pca_log_fvf.png"
ggsave(file, ax, width = 10, height = 5)

# Plot: Permutation tests
file <- "tmp/scrna/cellbench/cellbench-datasets.rds"
datasets <- readRDS(file)

rvp_with <- rvp(datasets$bal, "batch", "celltype", nperm = 1000)
rvp_without <- rvp(datasets$negctrl_bal, "batch", "celltype", nperm = 1000)

file <- "tmp/scrna/cellbench/rvp_with-permutations.rds"
rvp_with <- readRDS(file)
file <- "tmp/scrna/cellbench/rvp_without-permutations.rds"
rvp_without <- readRDS(file)

ax1 <- data.frame(rvp = rvp_with$null.distribution) %>%
  ggplot() +
    geom_histogram(
      aes(x = rvp), fill = "lightblue", col = "black", size = 0.2
    ) +
    geom_vline(xintercept = rvp_with$RVP, col = "red") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # legend.title = element_blank(),
      # legend.key.size = unit(4, "mm"),
      # legend.spacing.y = unit(1, "pt")
    ) +
    labs(
      title = "Permutation null distribution (n = 1000)",
      subtitle = "Data with batch effects",
      x = "HVP", y = "Count"
    ) +
    annotate(
      geom = "text", x = rvp_with$RVP - 0.01, y = 650,
      label = sprintf("Observed (p = %.2f)", rvp_with$p.value),
      color = "red", cex = 1.7, angle = 90
    )
file <- "tmp/fig/permutations-with.pdf"
ggsave(file, ax1, width = 3.6, height = 2)

ax2 <- data.frame(rvp = rvp_without$null.distribution) %>%
  ggplot() +
    geom_histogram(
      aes(x = rvp), fill = "lightblue", col = "black", size = 0.2
    ) +
    geom_vline(xintercept = rvp_without$RVP, col = "red") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # legend.title = element_blank(),
      # legend.key.size = unit(4, "mm"),
      # legend.spacing.y = unit(1, "pt")
    ) +
    labs(
      title = "Permutation null distribution (n = 1000)",
      subtitle = "Data without batch effects",
      x = "HVP", y = "Count"
    ) +
    annotate(
      geom = "text", x = rvp_without$RVP - 2e-4, y = 65,
      label = sprintf("Observed (p = %.2f)", rvp_without$p.value),
      color = "red", cex = 1.7, angle = 90
    )
file <- "tmp/fig/permutations-without.pdf"
ggsave(file, ax2, width = 3.6, height = 2)
