library(SingleCellExperiment)
library(Seurat)
library(ggVennDiagram)

varnames <- load("data/cellbench/cellbench.RData")
list_sce <- list(sce_sc_10x_qc, sce_sc_CELseq2_qc, sce_sc_Dropseq_qc)
platforms <- c("10x", "CELseq2", "Dropseq")
list_seurat <- lapply(
  list_sce, as.Seurat,
  counts = "counts", data = "logcounts"
)
for (i in seq_along(platforms)) {
  Idents(list_seurat[[i]]) <- platforms[i]
}

featuresets <- lapply(list_seurat, rownames)
# ggVennDiagram(featuresets)
venndata <- process_data(Venn(featuresets))
list_features <- venndata@region$item
names(list_features) <- venndata@region$name
print(str(list_features))

cellbench <- merge(
  list_seurat[[1]],
  y = list_seurat[2:3],
  add.cell.ids = c("10x", "CELseq2", "Dropseq"),
  project = "CellBench"
)
cellbench@meta.data$orig.ident <- cellbench@active.ident
# different features present in each platform
cellbench@assays$originalexp@meta.features$in_all <-
  rownames(cellbench) %in% list_features[[7]]

file <- "data/cellbench/cellbench-seurat.rds"
saveRDS(cellbench, file)