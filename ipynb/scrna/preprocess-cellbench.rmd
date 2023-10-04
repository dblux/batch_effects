```{r}
library(SingleCellExperiment)
library(Seurat)
library(ggVennDiagram)

varnames <- load("data/cellbench/cellbench.RData")
list_sce <- list(sce_sc_10x_qc, sce_sc_CELseq2_qc, sce_sc_Dropseq_qc)
list_seurat <- lapply(list_sce, as.Seurat, counts = "counts")

platforms <- c("10x", "CELseq2", "Dropseq")
for (i in seq_along(platforms)) {
  Idents(list_seurat[[i]]) <- platforms[i]
}

# Different features are present in each platform
# Identifying features present in all platforms
featuresets <- lapply(list_seurat, rownames)
# ggVennDiagram(featuresets)
venndata <- process_data(Venn(featuresets))
list_features <- venndata@region$item
names(list_features) <- venndata@region$name

cellbench <- merge(
  list_seurat[[1]],
  y = list_seurat[2:3],
  add.cell.ids = c("10x", "CELseq2", "Dropseq"),
  merge.data = FALSE,
  project = "CellBench"
)
```
```{r}
cellbench@meta.data$orig.ident <- cellbench@active.ident
cellbench@meta.data$batch <- cellbench@active.ident
cellbench@assays$originalexp@meta.features$in_all <-
  rownames(cellbench) %in% list_features[[7]]
colnames(cellbench@meta.data)[16] <- "celltype"
```
```{r}
file <- "data/cellbench/cellbench-seurat.rds"
saveRDS(cellbench, file)
```