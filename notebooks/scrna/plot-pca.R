library(magrittr)
library(ggplot2)
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


# Plot: Log-normalised values
# villani
dataname <- "villani"
file <- sprintf("tmp/scrna/%s/%s-datasets.rds", dataname, dataname)
datasets <- readRDS(file)
print(datasets$bal)

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
ax_batch <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Villani et al.",
METADATA_SID <- "../../relapse_prediction/data/GSE67684/processed/metadata/sid-metadata_v2.tsv"
METADATA_PID <- "../../relapse_prediction/data/GSE67684/processed/metadata/pid-metadata_v7.tsv"
metadata_sid <- read.table(METADATA_SID, sep = "\t")
metadata_pid <- read.table(METADATA_PID, sep = "\t", row.names = 1, quote = '"')

## Data
# Removed outliers, patients with timepoints from different batches and batch 5
SUBSET_RPATH <- "../../relapse_prediction/data/GSE67684/processed/subset_yeoh.tsv"
raw_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

# Metadata
metadata_sid$label <- as.factor(metadata_sid$label)
levels(metadata_sid$label) <- c('Remission', 'Relapse')
metadata_sid$batch_info <- as.factor(metadata_sid$batch_info) 
metadata_pid$label <- as.factor(metadata_pid$label)
levels(metadata_pid$label) <- c('Remission', 'Relapse')

# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- normaliseMeanScaling(raw_yeoh)
selected_yeoh <- removeProbesets(scaled_yeoh)
yeoh <- log2_transform(filterProbesets(selected_yeoh, 0.7, metadata_sid))

# All features
yeoh_allps <- log2_transform(scaled_yeoh)
yeoh_unfltr <- log2_transform(selected_yeoh)

metadata_yeoh <- metadata_sid[colnames(yeoh), ]
metadata_telaml1 <- subset(
  metadata_sid,
  subtype == 'TEL-AML1' &
  label == 'Remission' &
  batch_info %in% c(2, 9)
)

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
file <- sprintf("tmp/fig/%s/pca-%s.jpg", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 4, height = 2.3)
print(file)

# halfmix
dataname <- "halfmix"
file <- sprintf("tmp/scrna/%s/%s-datasets.rds", dataname, dataname)
datasets <- readRDS(file)
datasets$bal@meta.data$batch <- 
  factor(datasets$bal@meta.data$batch)
datasets$bal@meta.data$celltype <- 
  factor(datasets$bal@meta.data$celltype)
levels(datasets$bal@meta.data$batch) <- c("293T", "Jurkat", "Mix")
levels(datasets$bal@meta.data$celltype) <- c("293T", "Jurkat")
print(datasets$bal)
print(colnames(datasets$bal@meta.data))

custom_theme <- theme(
  title = element_text(size = 6),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  legend.key.size = unit(2.1, "mm"),
  legend.spacing.x = unit(0.5, "mm"),
  legend.position = "bottom", 
  # legend.justification = c("left", "bottom"),
  # legend.box.just = "left",
  legend.key = element_rect(fill = "transparent"),
  legend.background = element_rect(fill = "transparent")
)
ax_batch <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Zheng et al.",
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
file <- sprintf("tmp/fig/%s/pca-%s.jpg", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 4, height = 2.3)
print(file)

# cellbench
dataname <- "cellbench"
file <- sprintf("tmp/scrna/%s/%s-datasets.rds", dataname, dataname)
datasets <- readRDS(file)
print(datasets$bal)
print(colnames(datasets$bal@meta.data))

custom_theme <- theme(
  title = element_text(size = 6),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  legend.key.size = unit(2.1, "mm"),
  legend.spacing.x = unit(0.5, "mm"),
  legend.position = "bottom", 
  # legend.justification = c("left", "bottom"),
  # legend.box.just = "left",
  legend.key = element_rect(fill = "transparent"),
  legend.background = element_rect(fill = "transparent")
)
ax_batch <- ggplot_pca(
  GetAssayData(datasets$bal, layer = "data"),
  datasets$bal@meta.data,
  col = "batch",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Tian et al.",
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
file <- sprintf("tmp/fig/%s/pca-%s.jpg", dataname, dataname)
ggsave(file, ax_batch + ax_celltype, width = 4, height = 2.3)
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
metadata_without[3, 1] <- 2
metadata_without[6, 1] <- 1

custom_theme <- theme(
  title = element_text(size = 6),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  legend.key.size = unit(2.1, "mm"),
  legend.spacing.x = unit(0.5, "mm"),
  legend.position = "bottom", 
  # legend.justification = c("left", "bottom"),
  # legend.box.just = "left",
  legend.key = element_rect(fill = "transparent"),
  legend.background = element_rect(fill = "transparent")
)
ax <- ggplot_pca(
  maqc_bal,
  metadata_bal,
  col = "batch",
  pch = "class",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "MAQC",
    color = "Batch",
    pch = "Class"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca-%s.jpg", dataname, dataname)
ggsave(file, ax, width = 2, height = 2.3)
print(file)

# Yeoh et al.
METADATA_SID <- "../../relapse_prediction/data/GSE67684/processed/metadata/sid-metadata_v2.tsv"
METADATA_PID <- "../../relapse_prediction/data/GSE67684/processed/metadata/pid-metadata_v7.tsv"
metadata_sid <- read.table(METADATA_SID, sep = "\t")
metadata_pid <- read.table(METADATA_PID, sep = "\t", row.names = 1, quote = '"')

## Data
# Removed outliers, patients with timepoints from different batches and batch 5
SUBSET_RPATH <- "../../relapse_prediction/data/GSE67684/processed/subset_yeoh.tsv"
raw_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

# Metadata
metadata_sid$label <- as.factor(metadata_sid$label)
levels(metadata_sid$label) <- c('Remission', 'Relapse')
metadata_sid$batch_info <- as.factor(metadata_sid$batch_info) 
metadata_pid$label <- as.factor(metadata_pid$label)
levels(metadata_pid$label) <- c('Remission', 'Relapse')

# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- normaliseMeanScaling(raw_yeoh)
selected_yeoh <- removeProbesets(scaled_yeoh)
yeoh <- log2_transform(filterProbesets(selected_yeoh, 0.7, metadata_sid))

# All features
yeoh_allps <- log2_transform(scaled_yeoh)
yeoh_unfltr <- log2_transform(selected_yeoh)

metadata_yeoh <- metadata_sid[colnames(yeoh), ]
metadata_telaml1 <- subset(
  metadata_sid,
  subtype == 'TEL-AML1' &
  label == 'Remission' &
  batch_info %in% c(2, 9)
)


# westlake
dataname <- "westlake"
file1 <- "data/westlake/processed/balanced.csv"
file2 <- "data/westlake/processed/metadata/balanced.csv"
data <- read.csv(file1, row.names = 1)
metadata <- read.csv(file2, row.names = 1)
metadata$machine <- as.factor(metadata$machine)

stopifnot(identical(colnames(data), rownames(metadata)))
custom_theme <- theme(
  title = element_text(size = 6),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  legend.key.size = unit(2.1, "mm"),
  legend.spacing.x = unit(0.5, "mm"),
  legend.position = "bottom", 
  # legend.justification = c("left", "bottom"),
  # legend.box.just = "left",
  legend.key = element_rect(fill = "transparent"),
  legend.background = element_rect(fill = "transparent")
)
ax <- ggplot_pca(
  data,
  metadata,
  col = "machine",
  pch = "class",
  show.legend = TRUE,
  plot.axis = FALSE,
  cex = 0.8, alpha = 0.7
) +
  labs(
    title = "Wang et al.",
    color = "Batch",
    pch = "Class"
  ) +
  scale_color_manual(values = batch_cols) +
  custom_theme
file <- sprintf("tmp/fig/%s/pca-%s.jpg", dataname, dataname)
ggsave(file, ax, width = 2, height = 2.3)
print(file)
