library(magrittr)
library(Seurat) # v4.3.0
library(kBET)


# Data set: CellBench
# Counts data
# QC metrics calculated using scPipe
# Doublets identified using demuxlet
file <- "data/cellbench/cellbench-seurat.rds"
raw_cellbench <- readRDS(file)
print(dim(raw_cellbench))
print(table(raw_cellbench$batch, raw_cellbench$celltype))
print(str(raw_cellbench@meta.data))

# QC: Filtering cells
FeatureScatter(
  raw_cellbench,
  feature1 = "non_mt_percent", feature2 = "non_ribo_percent"
)
raw_cellbench <- subset(
  raw_cellbench,
  subset = demuxlet_cls == "SNG" &
    non_mt_percent > 0.8 &
    non_ribo_percent > 0.75
)

# Normalize data
cellbench <- NormalizeData(raw_cellbench)

# Batch effects detection
## PCA
cellbench <- cellbench %>%
  FindVariableFeatures() %>%
  ScaleData(do.scale = FALSE) %>%
  RunPCA(verbose = FALSE)
# cellbench <- FindVariableFeatures(cellbench)
# cellbench <- ScaleData(cellbench, do.scale = FALSE)
# cellbench <-  RunPCA(cellbench)
print(str(cellbench))

pca_batch <- DimPlot(cellbench, reduction = "pca", group.by = "batch")
pca_celltype <- DimPlot(cellbench, reduction = "pca", group.by = "celltype")
pca_batch + pca_celltype

## UMAP
cellbench <- RunUMAP(cellbench, dims = 1:50)
umap_batch <- DimPlot(cellbench, reduction = "umap", group.by = "batch")
umap_celltype <- DimPlot(cellbench, reduction = "umap", group.by = "celltype")
umap_batch + umap_celltype

# Batch effects metrics
X_pca <- Embeddings(cellbench, reduction = "pca")
kbet <- kBET(
  X_pca,
  cellbench$batch,
  do.pca = FALSE,
  k0 = NULL,
  testSize = ncol(cellbench),
  n_repeat = 1,
  verbose = TRUE
)
print(str(kbet))
rejection_rate <- kbet$summary$kBET.observed[1]
print(rejection_rate)

# Correction of batch effects 
## Seurat
list_cellbench <- SplitObject(raw_cellbench, split.by = "batch")
list_cellbench <- lapply(
  list_cellbench,
  function(obj) FindVariableFeatures(NormalizeData(obj))
)
features <- SelectIntegrationFeatures(list_cellbench) 
anchor_set <- FindIntegrationAnchors(
  list_cellbench, anchor.features = features
)
corr_cellbench <- IntegrateData(anchor_set)
print(corr_cellbench)

# Downstream analysis
## Plots
corr_cellbench <- corr_cellbench %>%
  ScaleData(do.scale = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:50)

pca_celltype <- DimPlot(
  corr_cellbench, reduction = "pca", group.by = "celltype"
)
pca_celltype
umap_batch <- DimPlot(corr_cellbench, reduction = "umap", group.by = "batch")
umap_celltype <- DimPlot(
  corr_cellbench, reduction = "umap", group.by = "celltype"
)
umap_batch + umap_celltype

# Batch effects metric: kBET
X_pca <- Embeddings(corr_cellbench, reduction = "pca")
kbet <- kBET(
  X_pca,
  cellbench$batch,
  do.pca = FALSE,
  k0 = NULL,
  testSize = ncol(cellbench),
  n_repeat = 1,
  verbose = TRUE
)
print(str(kbet))
rejection_rate <- kbet$summary$kBET.observed[1]
print(rejection_rate)

## Clustering
corr_cellbench <- corr_cellbench %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.2)
umap_celltype <- DimPlot(
  corr_cellbench, reduction = "umap", group.by = "seurat_clusters"
)
umap_celltype

# Differential expression analysis
markers <- FindMarkers(
  corr_cellbench, group.by = "seurat_clusters",
  ident.1 = "0", ident.2 = "1", verbose = FALSE
)
print(head(markers))

FeaturePlot(
  corr_cellbench,
  features = rownames(markers)[1:6]
)
