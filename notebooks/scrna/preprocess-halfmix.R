library(Seurat)
library(magrittr)


jurkat_counts <- Read10X(
  "data/jurkat_293t/raw/jurkat/filtered_matrices_mex/hg19",
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)
jurkat <- CreateSeuratObject(counts = jurkat_counts)

hek_counts <- Read10X(
  "data/jurkat_293t/raw/293t/filtered_matrices_mex/hg19",
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)
hek <- CreateSeuratObject(counts = hek_counts)

both_counts <- Read10X(
  "data/jurkat_293t/raw/jurkat_293t/filtered_matrices_mex/hg19",
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)
both <- CreateSeuratObject(counts = both_counts)

clusters <- read.csv(paste0(
  "data/jurkat_293t/raw/jurkat_293t/analysis_csv/",
  "kmeans/2_clusters/clusters.csv"
))

# Edit metadata
jurkat$celltype <- "jurkat"
jurkat$orig.ident <- "jurkat"
jurkat$batch <- jurkat$orig.ident

hek$celltype <- "293t"
hek$orig.ident <- "293t"
hek$batch <- hek$orig.ident

print(identical(clusters$Barcode, colnames(both)))
both$orig.ident <- "zheng"
both$batch <- both$orig.ident
# cluster 1 consists of 293t cells
both$celltype <- ifelse(clusters$Cluster == 1, "293t", "jurkat")

# Merge datasets
# roughly 60 colnames are duplicated
halfmix <- merge(both, c(hek, jurkat), project = "Jurkat-HEK293T")
Idents(halfmix) <- "orig.ident"


saveRDS(halfmix, "data/jurkat_293t/processed/halfmix.rds")
