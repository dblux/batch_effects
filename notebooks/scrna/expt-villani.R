library(Seurat)
library(scater)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

library(CellMixS)
library(kBET)
library(lisi)

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}

# villani
# - TPM values (Smartseq2)
# - Two batches: **Different plates**
# - Four classes: Human dendritic cell lines

file <- "data/villani/processed/villani-seurat.rds"
villani <- readRDS(file)

# QC - Filter cells
# - No count data present, hence no nCounts, etc.

villani$nFeature_tpm <- colSums(GetAssayData(villani) != 0)
mito_genes <- rownames(villani)[grep("^MT-", rownames(villani))]
print(mito_genes) # no mito genes
ribo_genes <- rownames(villani)[grep("^RP[SL]", rownames(villani))]
total_ribo <- colSums(GetAssayData(villani[ribo_genes, ]))
villani$percent.ribo <- total_ribo / colSums(GetAssayData(villani))

# FeatureScatter(
#   villani_sub,
#   feature1 = "nFeature_tpm", feature2 = "percent.ribo",
#   group.by = "batch"
# )

villani_sub <- subset(
  villani, subset = percent.ribo > 0.05 &
    nFeature_tpm < 8000 &
    nFeature_tpm > 2000
)
table(villani_sub$batch, villani_sub$celltype)

# # QC - Filter features
# - Filter out ribo genes
# - Filter out sparse genes

sparse_genes <- remove_sparse(
  GetAssayData(villani_sub, "data"),
  0.95, villani_sub$celltype,
  ret.features = TRUE
)
rm_genes <- union(ribo_genes, sparse_genes)
villani_sel <- villani_sub[!(rownames(villani_sub) %in% rm_genes), ]

logcnts <- GetAssayData(villani_sel)
pct_zero <- rowSums(logcnts == 0) / ncol(logcnts) 
pdf("tmp/fig/villani-missingness.pdf")
hist(pct_zero)
dev.off()

sum(logcnts == 0) / prod(dim(villani_sel))

# Log transform data
villani_sel@assays$log_tpm <- villani_sel@assays$tpm
villani_sel <- SetAssayData(
  villani_sel, slot = "data", assay = "log_tpm",
  new.data = log1p(GetAssayData(villani_sel))
)
DefaultAssay(villani_sel) <- "log_tpm"

# Subset
tabnames <- list(unique(villani$celltype), unique(villani$batch))
ct_bal <- matrix(30, 4, 2, dimnames = tabnames)
ct_imbal <- matrix(
  c(15, 45, 20, 40, 45, 15, 40, 20),
  4, 2, dimnames = tabnames
)
idx_bal <- idx_imbal <- numeric()
for (g in rownames(ct_bal)) {
  for (k in colnames(ct_bal)) {
    j_kg <- which(villani_sel$batch == k & villani_sel$celltype == g)
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
villani_bal <- villani_sel[, idx_bal]
villani_imbal <- villani_sel[, idx_imbal]
table(villani_bal$celltype, villani_bal$batch)
table(villani_imbal$celltype, villani_imbal$batch)

villani_b1 <- subset(villani_sel, subset = batch == 1)
table(villani_b1$celltype, villani_b1$batch)

idx_bal <- idx_imbal <- numeric()
batch_bal <- batch_imbal <- rep(NA, ncol(villani_b1))
for (g in rownames(ct_bal)) {
  j1_g <- j2_g <- which(villani_b1$celltype == g)
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
villani_b1_bal <- villani_b1[, idx_bal]
villani_b1_bal$batch <- batch_bal[idx_bal]
villani_b1_imbal <- villani_b1[, idx_imbal]
villani_b1_imbal$batch <- batch_imbal[idx_imbal]
table(villani_b1_bal$celltype, villani_b1_bal$batch)
table(villani_b1_imbal$celltype, villani_b1_imbal$batch)

# Save data sets 
villani_objs <- list(
  negctrl_bal = villani_b1_bal,
  negctrl_imbal = villani_b1_imbal,
  bal = villani_bal,
  imbal = villani_imbal
)
# saveRDS(villani_objs, "tmp/villani-datasets.rds")

# Analyse saved results
objs <- readRDS("tmp/halfmix-results_k500.rds")
for (idx in names(objs)) {
  cat(idx, fill = T)
  obj <- objs[[idx]]
  cat(paste("k =", obj$kbet$params$k0), fill = T)
  cat(paste("cms:", mean(obj$cms$cms)), fill = T)
  cat(paste("kbet:", obj$kbet$summary$kBET.observed[1]), fill = T)
  cat(paste("blisi:", mean(obj$lisi$batch)), fill = T)
  cat(fill = T)
}


# Plot
villani_sel1 <- villani_sel %>%
  ScaleData(do.scale = FALSE) %>%
  RunPCA()
ax1 <- DimPlot(
  villani_sel1, reduction = "pca", group.by = "batch", shuffle = TRUE
)
ax2 <- DimPlot(
  villani_sel1, reduction = "pca", group.by = "celltype", shuffle = TRUE
)
ax <- plot_grid(ax1, ax2, rel_widths = c(0.9, 1))
file <- "tmp/pca-villani.png"
ggsave(file, ax, width = 10, height = 5)


for (idx in names(villani_objs)) {
  villani_sel1 <- villani_objs[[idx]] %>%
    ScaleData(do.scale = FALSE) %>%
    RunPCA(features = rownames(.))
  ax1 <- DimPlot(
    villani_sel1, reduction = "pca",
    group.by = "batch", shuffle = TRUE
  )
  ax2 <- DimPlot(
    villani_sel1, reduction = "pca",
    group.by = "celltype", shuffle = TRUE
  )
  ax <- plot_grid(ax1, ax2, rel_widths = c(1, 1))
  file <- sprintf("tmp/villani-pca_%s.png", idx)
  ggsave(file, ax, width = 10, height = 5)
}


ax1 <- DimPlot(
  villani_sel, reduction = "pca",
  group.by = "batch", shuffle = TRUE
)
ax2 <- DimPlot(
  villani_sel, reduction = "pca",
  group.by = "celltype", shuffle = TRUE
)
ax <- plot_grid(ax1, ax2, rel_widths = c(0.9, 1))
ax
# file <- "tmp/pca-villani_sel.png"
# ggsave(file, ax, width = 10, height = 5)
