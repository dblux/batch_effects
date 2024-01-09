library(magrittr)
library(ggplot2)
library(cowplot)
library(Seurat)
theme_set(theme_bw(base_size = 7))

# source files
src_files <- list.files("R", full.names = TRUE)
for (f in src_files) {
  source(f)
  cat(sprintf("Sourced file: %s\n", f))
}


# Load real data
# Measure batch effects and without batch effects
dataname <- "halfmix"
file <- sprintf("tmp/scrna/%s/%s-datasets.rds", dataname, dataname)
print(file)
datasets <- readRDS(file)

# With batch effects
rvp_obj <- rvp(datasets$bal, "batch", "celltype", ret.percent = FALSE)

ax <- plot_rvp(rvp_obj) +
  labs(title = dataname)
file <- sprintf("tmp/fig/rvp_total-%s.jpg", dataname)
ggsave(file, ax, width = 6, height = 2)

# ss_ord <- ss[rev(order(ss$ss_batch / ss$ss_total)), ]
# ax <- plot_rvp(ss_ord)
# file <- sprintf("tmp/fig/rvp-%s.jpg", dataname)
# ggsave(file, ax, width = 6, height = 2)

# Without batch effects
rvp_obj <- rvp(datasets$negctrl_bal, "batch", "celltype", ret.percent = FALSE)

ax <- plot_rvp(rvp_obj)
file <- sprintf("tmp/fig/rvp_total-%s_without.jpg", dataname)
ggsave(file, ax, width = 6, height = 2)

# ss_ord <- ss[rev(order(ss$ss_batch / ss$ss_total)), ]
# ax <- plot_rvp(ss_ord)
# file <- sprintf("tmp/fig/rvp-%s_without.jpg", dataname)
# ggsave(file, ax, width = 6, height = 2)


# Microarray
file_maqc <- "data/MAQC-I/processed/mas5_original-ref.tsv"
raw_maqc <- read.table(file_maqc, sep = "\t", header = T, row.names = 1)

# MAQC metadata
# Class A - Universal Human Reference RNA (UHRR)
# Class B - Human Brain Reference RNA (HBRR)
batch <- as.factor(rep(1:6, each = 10))
class <- rep(rep(LETTERS[1:2], each = 5), 6)
metadata_maqc <- data.frame(
  batch, class,
  row.names = colnames(raw_maqc)
)

# SCALE->REMOVE->FILTER->LOG
log_maqc <- raw_maqc %>%
  scale_trimmed() %>%
  remove_affymetrix() %>%
  remove_sparse(0.3, metadata_maqc$class) %>%
  log2_transform()

# subsetting maqc
maqc_bal <- log_maqc[, c(1:4, 6:9, 11:14, 16:19)]
maqc_fakebatch <- log_maqc[, c(1:4, 6:9)]

metadata_bal <- metadata_maqc[colnames(maqc_bal), ]
fake_metadata <- metadata_maqc[colnames(maqc_fakebatch), ]
fake_metadata[3, 1] <- 2
fake_metadata[6, 1] <- 1

dataname <- "maqc"

# With batch effects
rvp_obj <- rvp(
  t(maqc_bal), metadata_bal$batch, metadata_bal$class,
  ret.percent = FALSE
)

ax <- plot_rvp(rvp_obj)
file <- sprintf("tmp/fig/rvp_total-%s.jpg", dataname)
ggsave(file, ax, width = 6, height = 2)

# ss_ord <- ss[rev(order(ss$ss_batch / ss$ss_total)), ]
# print(head(ss_ord))
# ax <- plot_rvp(ss_ord)
# file <- sprintf("tmp/fig/rvp-%s.jpg", dataname)
# print(file)
# ggsave(file, ax, width = 6, height = 2)

# Without batch effects
rvp_obj <- rvp(
  t(maqc_fakebatch), fake_metadata$batch, fake_metadata$class,
  ret.percent = FALSE
)

ax <- plot_rvp(rvp_obj)
file <- sprintf("tmp/fig/rvp_total-%s_without.jpg", dataname)
ggsave(file, ax, width = 6, height = 2)
