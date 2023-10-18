library(Seurat)
library(magrittr)


file <- "data/villani/processed/villani-sce.rds"
villani <- readRDS(file)
vill <- as.Seurat(villani, counts = NULL, data = "tpm")
vill <- RenameAssays(vill, originalexp = "tpm")
colnames(vill@meta.data)[2] <- "celltype"
Idents(vill) <- vill@meta.data$batch
vill <- vill %>%
  FindVariableFeatures() %>%
  ScaleData(do.scale = FALSE) %>%
  RunPCA()
ofile <- "data/villani/processed/villani-seurat.rds"
saveRDS(vill, ofile)