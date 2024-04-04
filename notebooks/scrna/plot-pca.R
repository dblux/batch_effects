library(ggplot2)
library(magrittr)
library(RColorBrewer)
library(Seurat)
theme_set(theme_bw(base_size = 7))
# source files
src_files <- list.files("R", full.names = TRUE)
for (f in src_files) {
  source(f)
  cat(sprintf("Sourced file: %s\n", f))
}


batch_cols <- brewer.pal(8, "Dark2")[2:4]
class_cols <- brewer.pal(8, "Dark2")[5:8]
custom_theme <- theme(
  title = element_text(size = 6),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  legend.key.size = unit(2.1, "mm"),
  legend.spacing.x = unit(0.5, "mm"),
  legend.position = "bottom", 
  legend.key = element_rect(fill = "transparent"),
  legend.background = element_rect(fill = "transparent")
)


# Plot: Log-normalised values
# villani
dataname <- "villani"
file <- sprintf("tmp/scrna/%s/%s-datasets.rds", dataname, dataname)
datasets <- readRDS(file)
print(datasets$bal)

# with batch effects
ax_batch <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Villani et al. (scRNA-seq)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "celltype",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Cell type"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme +
  guides(col = guide_legend(nrow = 2))
file <- sprintf("tmp/fig/%s/pca-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# without batch effects
ax_batch <- ggplot_pca(
  GetAssayData(datasets$negctrl_bal, layer = "data"),
  datasets$negctrl_bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Villani et al. (scRNA-seq)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  GetAssayData(datasets$negctrl_bal, layer = "data"),
  datasets$negctrl_bal@meta.data,
  col = "celltype",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Cell type"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme +
  guides(col = guide_legend(nrow = 2))
file <- sprintf("tmp/fig/%s/pca_without-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# halfmix
dataname <- "halfmix"
file <- sprintf("tmp/scrna/%s/%s-datasets.rds", dataname, dataname)
datasets <- readRDS(file)
# with
datasets$bal@meta.data$batch <- 
  factor(datasets$bal@meta.data$batch)
datasets$bal@meta.data$celltype <- 
  factor(datasets$bal@meta.data$celltype)
levels(datasets$bal@meta.data$batch) <- c("293T", "Jurkat", "Mix")
levels(datasets$bal@meta.data$celltype) <- c("293T", "Jurkat")
# without
datasets$negctrl_bal@meta.data$batch <- 
  factor(datasets$negctrl_bal@meta.data$batch)
datasets$negctrl_bal@meta.data$celltype <- 
  factor(datasets$negctrl_bal@meta.data$celltype)
levels(datasets$negctrl_bal@meta.data$batch) <- c("293T", "Jurkat", "Mix")
levels(datasets$negctrl_bal@meta.data$celltype) <- c("293T", "Jurkat")
print(datasets$negctrl_bal)

# with
ax_batch <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Zheng et al. (scRNA-seq)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "celltype",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Cell type"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# without
ax_batch <- ggplot_pca(
  GetAssayData(datasets$negctrl_bal, layer = "data"),
  datasets$negctrl_bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Zheng et al. (scRNA-seq)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  GetAssayData(datasets$negctrl_bal, layer = "data"),
  datasets$negctrl_bal@meta.data,
  col = "celltype",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Cell type"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca_without-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# cellbench
dataname <- "cellbench"
file <- sprintf("tmp/scrna/%s/%s-datasets.rds", dataname, dataname)
datasets <- readRDS(file)
print(datasets$bal)
print(colnames(datasets$bal@meta.data))

# with
ax_batch <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Tian et al. (scRNA-seq)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "celltype",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Cell type"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# without
ax_batch <- ggplot_pca(
  GetAssayData(datasets$negctrl_bal, layer = "data"),
  datasets$negctrl_bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Tian et al. (scRNA-seq)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  GetAssayData(datasets$negctrl_bal, layer = "data"),
  datasets$negctrl_bal@meta.data,
  col = "celltype",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Cell type"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca_without-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# MAQC
dataname <- "maqc"
file_maqc <- "data/MAQC-I/processed/mas5_original-ref.tsv"
raw_maqc <- read.table(file_maqc, sep = "\t", header = T, row.names = 1)

# MAQC metadata
# Class A - Universal Human Reference RNA (UHRR)
# Class B - Human Brain Reference RNA (HBRR)
batch <- as.factor(rep(1:6, each = 10))
class <- rep(rep(c("UHRR", "HBRR"), each = 5), 6)
metadata_maqc <- data.frame(
  batch, class,
  row.names = colnames(raw_maqc)
)

log_maqc <- raw_maqc %>%
  scale_trimmed() %>%
  remove_affymetrix() %>%
  remove_sparse(0.3, metadata_maqc$class) %>%
  log2_transform()

# subsetting maqc
maqc_bal <- log_maqc[, c(1:4, 6:9, 11:14, 16:19)]
maqc_without <- log_maqc[, c(1:4, 6:9)]
metadata_bal <- metadata_maqc[colnames(maqc_bal), ]
metadata_without <- metadata_maqc[colnames(maqc_without), ]
metadata_without[c(3, 4, 7, 8), 1] <- 2

# with
ax_batch <- ggplot_pca(
  maqc_bal,
  metadata_bal,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "MAQC (microarray)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  maqc_bal,
  metadata_bal,
  col = "class",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Class"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# without
ax_batch <- ggplot_pca(
  maqc_without,
  metadata_without,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "MAQC (microarray)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  maqc_without,
  metadata_without,
  col = "class",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Class"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca_without-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# Yeoh et al.
dataname <- "yeoh"
file1 <- "data/GSE67684/processed/metadata/sid-metadata_v2.tsv"
file2 <- "data/GSE67684/processed/metadata/pid-metadata_v7.tsv"
metadata_sid <- read.table(file1, sep = "\t")
metadata_pid <- read.table(file2, sep = "\t", row.names = 1, quote = '"')

## Data
# Removed outliers, patients with timepoints from different batches and batch 5
SUBSET_RPATH <- "data/GSE67684/processed/subset_yeoh.tsv"
raw_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

# Metadata
metadata_sid$label <- as.factor(metadata_sid$label)
levels(metadata_sid$label) <- c("Remission", "Relapse")
metadata_sid$batch <- as.factor(metadata_sid$batch) 
metadata_pid$label <- as.factor(metadata_pid$label)
levels(metadata_pid$label) <- c("Remission", "Relapse")
metadata_yeoh <- metadata_sid[colnames(raw_yeoh), ]

# SCALE->REMOVE->FILTER->LOG
log_yeoh <- raw_yeoh %>%
  scale_trimmed() %>%
  remove_affymetrix() %>%
  remove_sparse(0.3, metadata_yeoh$class) %>%
  log2_transform()

bal_sids <- get_sid(c(
  "P022", "P023", "P024", "P025",
  "P099", "P106", "P120", "P121"
))
yeoh_bal <- log_yeoh[, bal_sids]
metadata_telaml1 <- metadata_yeoh[bal_sids, ]

without_sids <- get_sid(c(
  "P022", "P023", "P024", "P025",
  "P026", "P027", "P035", "P036"
))
yeoh_without <- log_yeoh[, without_sids]
metadata_without <- metadata_sid[without_sids, ]
metadata_without[c(1:4, 9:12), 1] <- 1

# metadata_telaml1 <- subset(
#   metadata_yeoh,
#   subtype == "TEL-AML1" &
#   label == "Remission" &
#   batch %in% c(2, 9)
# )

# imbal_sids <- c(
#   "P022_D0", "P023_D0", "P024_D0", "P025_D0", "P026_D0",
#   "P022_D8", "P023_D8", "P024_D8", "P099_D0", "P106_D0",
#   "P120_D0", "P099_D8", "P106_D8", "P120_D8", "P121_D8", "P127_D8"
# )
# yeoh_imbal <- yeoh[, imbal_sids]

# with
ax_batch <- ggplot_pca(
  yeoh_bal,
  metadata_telaml1,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Yeoh et al. (microarray)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  yeoh_bal,
  metadata_telaml1,
  col = "class",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Class"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# without
ax_batch <- ggplot_pca(
  yeoh_without,
  metadata_without,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Yeoh et al. (microarray)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  yeoh_without,
  metadata_without,
  col = "class",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Class"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca_without-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# westlake
dataname <- "westlake"
# with - bal
file1 <- "data/westlake/processed/sparse/with-bal.csv"
file2 <- "data/westlake/processed/metadata/with-balanced.csv"
with_bal <- read.csv(file1, row.names = 1)
metadata_with_bal <- read.csv(file2, row.names = 1)
# with - imbal
file1 <- "data/westlake/processed/sparse/with-imbal.csv"
file2 <- "data/westlake/processed/metadata/with-severe1.csv"
with_imbal <- read.csv(file1, row.names = 1)
metadata_with_imbal <- read.csv(file2, row.names = 1)
# without - bal
file1 <- "data/westlake/processed/sparse/without-bal.csv"
file2 <- "data/westlake/processed/metadata/without-balanced1.csv"
without_bal <- read.csv(file1, row.names = 1)
metadata_without_bal <- read.csv(file2, row.names = 1)
# without - imbal
file1 <- "data/westlake/processed/sparse/without-imbal.csv"
file2 <- "data/westlake/processed/metadata/without-severe1.csv"
without_imbal <- read.csv(file1, row.names = 1)
metadata_without_imbal <- read.csv(file2, row.names = 1)
# change numeric to factor
metadata_with_bal$machine <- as.factor(metadata_with_bal$machine)
metadata_with_imbal$machine <- as.factor(metadata_with_imbal$machine)
metadata_without_bal$machine <- as.factor(metadata_without_bal$machine)
metadata_without_imbal$machine <- as.factor(metadata_without_imbal$machine)

# table(metadata_with_bal$class, metadata_with_bal$machine)
# table(metadata_with_imbal$class, metadata_with_imbal$machine)
# table(metadata_without_bal$class, metadata_without_bal$machine)
# table(metadata_without_imbal$class, metadata_without_imbal$machine)

# plot: PCA
batch_cols <- brewer.pal(8, "Dark2")[2:4]
class_cols <- brewer.pal(8, "Dark2")[5:8]
custom_theme <- theme(
  title = element_text(size = 6),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  legend.key.size = unit(2.1, "mm"),
  legend.spacing.x = unit(0.5, "mm"),
  legend.position = "bottom", 
  legend.key = element_rect(fill = "transparent"),
  legend.background = element_rect(fill = "transparent")
)

# with
ax_batch <- ggplot_pca(
  with_bal,
  metadata_with_bal,
  col = "machine",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Wang et al. (proteomics)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  with_bal,
  metadata_with_bal,
  col = "class",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Class"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)

# without
ax_batch <- ggplot_pca(
  without_bal,
  metadata_without_bal,
  col = "machine",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Wang et al. (proteomics)",
    color = "Batch"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
ax_celltype <- ggplot_pca(
  without_bal,
  metadata_without_bal,
  col = "class",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    color = "Class"
  ) +
  scale_color_manual(values = class_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca_without-%s.pdf", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 3.6, height = 2)
print(file)
