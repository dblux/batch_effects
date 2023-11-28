library(ggplot2)
theme_set(theme_bw(base_size = 8))

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}

# Simulate DEA data set
## With.Balanced
crosstab <- matrix(80, 2, 2)
data <- simulate_microarray(
  10000, crosstab, delta = 0.7,
  phi = 0.3, zeta = 0.7
)
file <- "data/simulated/microarray/dea/phi30/with-bal.rds"
saveRDS(data, file)

## With.Imbalanced
npercond <- c(40, 120, 120, 40)
crosstab <- matrix(npercond, 2, 2)
data <- simulate_microarray(
  10000, crosstab, delta = 0.7,
  phi = 0.3, zeta = 0.7
)
file <- "data/simulated/microarray/dea/phi30/with-imbal.rds"
saveRDS(data, file)


crosstab <- matrix(160, 2, 1)
data <- simulate_microarray(
  10000, crosstab, delta = 0.7,
  phi = 0.3, zeta = 0.7
)
# Without.Balanced
batch <- rep("1", 320)
batch[81:240] <- "1'"
data$metadata$batch <- as.factor(batch)
file <- "data/simulated/microarray/dea/phi30/without-bal.rds"
# saveRDS(data, file)

# Without.Imbalanced
batch <- rep("1", 320)
batch[41:200] <- "1'"
data$metadata$batch <- as.factor(batch)
file <- "data/simulated/microarray/dea/phi30/without-imbal.rds"
# saveRDS(data, file)

# Plot
## PCA
ax_pca <- ggplot_pca(
  data$X, data$metadata, cex = 1,
  col = 'batch', pch = 'class'
)
file <- "tmp/fig/pca-dea_phi30_without_imbal.png"
ggsave(file, ax_pca, width = 4, height = 3)

## Feature
i <- 0
i <- i + 1
print(i)
feature <- data.frame(
  value = as.vector(data.matrix(data$X[i, ])),
  batch = data$metadata$batch,
  class = data$metadata$class
)

ax <- ggplot(feature) +
  geom_point(
    aes(x = batch, y = value, pch = class, col = batch),
    position = position_jitterdodge(jitter.width = .5),
    show.legend = TRUE
  ) +
  ylim(0, 15)
file <- "tmp/fig/feature1-dea_phi30_with_bal.png"
ggsave(file, ax, width = 5, height = 5)


npercond <- c(400, 1200, 1200, 400)
crosstab <- matrix(npercond, 2, 2)
data <- simulate_microarray(
  10000, crosstab, delta = 0.7,
  phi = 0.3, zeta = 0.7
)
file <- "data/simulated/scrna/with_imbal.rds"
saveRDS(data, file)

idx <- sample(1:3200, 800)
ax_pca <- ggplot_pca(
  data$X[, idx], data$metadata[idx, ], cex = 1,
  col = 'batch', pch = 'class'
)
file <- "tmp/fig/scrna-pca_without_imbal.png"
ggsave(file, ax_pca, width = 4, height = 3)