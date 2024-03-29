library(Seurat)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())
library(biomaRt)
library(kBET)
library(lisi)

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


halfmix <- readRDS("data/jurkat_293t/processed/halfmix.rds")


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filter = "ensembl_gene_id",
  values = rownames(halfmix),
  mart = ensembl
)
hgnc_symbol <- mapping$hgnc_symbol[
  match(rownames(halfmix), mapping$ensembl_gene_id)]
mito_genes <- rownames(halfmix)[grep("^MT-", hgnc_symbol)]
ribo_genes <- rownames(halfmix)[grep("^RP[SL]", hgnc_symbol)]
halfmix$percent.mito <-
  colSums(GetAssayData(halfmix, slot = "counts")[mito_genes, ]) /
  halfmix$nCount_RNA
halfmix$percent.ribo <-
  colSums(GetAssayData(halfmix, slot = "counts")[ribo_genes, ]) /
  halfmix$nCount_RNA

# QC
- Only ribo genes present, no mito genes

# FeatureScatter(
#   halfmix_sub,
#   feature1 = "percent.mito", feature2 = "percent.ribo",
#   group.by = "batch"
# )
# FeatureScatter(
#   halfmix_sub,
#   feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
#   group.by = "batch"
# )
halfmix_sub <- subset(
  halfmix,
  subset = nFeature_RNA > 1500 &
    nCount_RNA < 3e4 &
    percent.mito < 0.1 &
    percent.ribo > 0.2
)
table(halfmix_sub$celltype, halfmix_sub$batch)

# Normalise data

halfmix_sub <- NormalizeData(halfmix_sub)

# Feature selection: Remove ribo and sparse genes

sparse_genes <- remove_sparse(
  GetAssayData(halfmix_sub, "counts"),
  0.95, halfmix_sub$celltype,
  ret.features = TRUE
)
rm_genes <- unique(c(mito_genes, ribo_genes, sparse_genes))
halfmix_sel <- halfmix_sub[!(rownames(halfmix_sub) %in% rm_genes), ]
dim(halfmix_sel)

# Subset
# with batch effects
tabnames <- list(unique(halfmix_sel$celltype), unique(halfmix_sel$batch))
ct_bal <- matrix(
  c(600, 600, 0, 600, 600, 0),
  2, 3, dimnames = tabnames
)
ct_imbal1 <- matrix(
  c(300, 900, 0, 300, 900, 0),
  2, 3, dimnames = tabnames
)
ct_imbal2 <- matrix(
  c(900, 300, 0, 900, 300, 0),
  2, 3, dimnames = tabnames
)
# with batch effects
idx_bal <- idx_imbal1 <- idx_imbal2 <- numeric()
for (g in rownames(ct_bal)) {
  for (k in colnames(ct_bal)) {
    j_kg <- which(halfmix_sel$batch == k & halfmix_sel$celltype == g)
    if (length(j_kg) == 0)
      next
    # balanced data set
    n1_kg <- ct_bal[g, k]
    j1_kg <- sample(j_kg, n1_kg)
    idx_bal <- c(idx_bal, j1_kg)
    # imbalanced data set
    n2_kg <- ct_imbal1[g, k]
    j2_kg <- sample(j_kg, n2_kg)
    idx_imbal1 <- c(idx_imbal1, j2_kg)
    # imbalanced data set
    n3_kg <- ct_imbal2[g, k]
    j3_kg <- sample(j_kg, n3_kg)
    idx_imbal2 <- c(idx_imbal2, j3_kg)
  }
}
halfmix_bal <- halfmix_sel[, idx_bal]
halfmix_imbal1 <- halfmix_sel[, idx_imbal1]
halfmix_imbal2 <- halfmix_sel[, idx_imbal2]
table(halfmix_bal$celltype, halfmix_bal$batch)
table(halfmix_imbal1$celltype, halfmix_imbal1$batch)
table(halfmix_imbal2$celltype, halfmix_imbal2$batch)


# without batch effects
tabnames <- list(unique(halfmix_sel$celltype), 1:3)
ct_bal <- matrix(
  c(600, 600, 0, 600, 600, 0),
  2, 3, dimnames = tabnames
)
ct_imbal <- matrix(
  c(300, 900, 0, 300, 900, 0),
  2, 3, dimnames = tabnames
)
halfmix_b1 <- subset(halfmix_sel, subset = batch == "zheng")

idx_bal <- idx_imbal <- numeric()
batch_bal <- batch_imbal <- rep(NA, ncol(halfmix_b1))
for (g in rownames(ct_bal)) {
  j1_g <- j2_g <- which(halfmix_b1$celltype == g)
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
halfmix_b1_bal <- halfmix_b1[, idx_bal]
halfmix_b1_bal$batch <- batch_bal[idx_bal]
halfmix_b1_imbal <- halfmix_b1[, idx_imbal]
halfmix_b1_imbal$batch <- batch_imbal[idx_imbal]
table(halfmix_b1_bal$celltype, halfmix_b1_bal$batch)
table(halfmix_b1_imbal$celltype, halfmix_b1_imbal$batch)

# Evaluate batch effects

halfmix_objs <- list(
  negctrl_bal = halfmix_b1_bal,
  negctrl_imbal = halfmix_b1_imbal,
  bal = halfmix_bal,
  imbal1 = halfmix_imbal1,
  imbal2 = halfmix_imbal2
)
saveRDS(halfmix_objs, "tmp/halfmix-datasets.rds")


file <- "tmp/halfmix-results_k0.rds"
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
## PCA
halfmix_objs <- readRDS("tmp/halfmix-datasets.rds")
for (idx in names(halfmix_objs)) {
  halfmix_sel1 <- halfmix_objs[[idx]] %>%
    ScaleData(do.scale = FALSE) %>%
    RunPCA(features = rownames(.))
  ax1 <- DimPlot(
    halfmix_sel1, reduction = "pca",
    group.by = "batch", shuffle = TRUE
  )
  ax2 <- DimPlot(
    halfmix_sel1, reduction = "pca",
    group.by = "celltype", shuffle = TRUE
  )
  ax <- plot_grid(ax1, ax2, rel_widths = c(1, 1))
  file <- sprintf("tmp/halfmix-pca_%s.png", idx)
  ggsave(file, ax, width = 10, height = 5)
}


halfmix_sub1 <- halfmix_sub %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(do.scale = FALSE) %>%
  RunPCA()
ax1 <- DimPlot(
  halfmix_sub1, reduction = "pca",
  group.by = "batch", shuffle = TRUE
)
ax2 <- DimPlot(
  halfmix_sub1, reduction = "pca",
  group.by = "celltype", shuffle = TRUE
)
ax <- plot_grid(ax1, ax2, rel_widths = c(1, 1))
file <- "tmp/fig/halfmix-pca_log_fvf.png"
ggsave(file, ax, width = 10, height = 5)


halfmix <- RunUMAP(halfmix, dims = 1:20)
ax <- DimPlot(halfmix, reduction = "umap")
file <- "tmp/umap-halfmix.png"
ggsave(file, ax, width = 5, height = 5)
