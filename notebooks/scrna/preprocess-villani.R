library(Seurat)

file <- "data/villani/raw/GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt.gz"
raw <- read.table(file, sep = "\t")
tpm <- raw[, 1:768]

# Metadata
info <- strsplit(colnames(tpm), split = "_")
metadata <- data.frame(do.call(rbind, info))
colnames(metadata) <- c("celltype", "patient", "cell")
batch <- ifelse(metadata$patient %in% c("P7", "P8", "P9", "P10"), 1, 2)
metadata$batch <- as.factor(batch)
rownames(metadata) <- colnames(tpm)

villani <- CreateSeuratObject(tpm, meta.data = metadata)
ofile <- "data/villani/processed/villani-seurat.rds"
saveRDS(villani, ofile)