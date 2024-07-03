library(scater)
library(harmony)


########## Data: Splat ##########

data <- readRDS("data/simulated/scrna/loc0/datasets-0.08-01.rds")
bal <- data$bal

# Plot: PCA (before correction)
bal <- runPCA(bal)
ax <- plotPCA(bal, colour_by = "Batch", shape_by = "Group")
file <- "tmp/fig/pca-before.pdf"
ggsave(file, ax, width = 5, height = 4)

# Run: Harmony (sce object)
bal <- RunHarmony(bal, group.by.vars = "Batch")
ax <- plotReducedDim(bal, "HARMONY", colour_by = "Batch", shape_by = "Group")
file <- "tmp/fig/pca-after.pdf"
ggsave(file, ax, width = 5, height = 4)

pca_bal <- reducedDim(bal)
y_bal <- colData(bal)
y_bal$cell_id <- rownames(y_bal)
# Run: Harmony (embeddings)
splat_corr <- RunHarmony(pca_bal, y_bal, "Batch")

ggplot(bal) +
  geom_point(aes(x = PC1, y = PC2))

########## Data: Jurkat ##########

data(cell_lines)
v <- cell_lines$scaled_pcs
metadata <- cell_lines$meta_data
corr_pca <- RunHarmony(v, metadata, "dataset")
