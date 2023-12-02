library(ggplot2)
library(magrittr)
theme_set(theme_bw(base_size = 8))

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


# Simulate: Different magnitudes
## Balanced
nclass <- 2
nbatch <- 2
npercond <- 20

crosstab <- matrix(npercond, nclass, nbatch)
print(crosstab)

for (delta in seq(0.1, 1.0, 0.1)) {
  simdata <- simulate_microarray(
    10000, crosstab,
    delta = delta,
    gamma = 0.5,
    phi = 0.2,
    c = 10, d = 6,
    epsilon = 0.5,
    kappa = 0.2,
    dropout = FALSE
  )
  file <- sprintf("data/simulated/microarray/magnitude/bal-%.02f.rds", delta)
  print(file)
  saveRDS(simdata, file)
  # Plot
  ## Normalised data
  ax <- 2 ^ simdata$X %>%
    scale_trimmed() %>%
    log2_transform() %>%
    ggplot_pca(simdata$metadata, col = 'batch', pch = 'class')
  file <- sprintf("tmp/fig/microarray/magnitude/pca-bal_%.02f.jpg", delta)
  ggsave(file, ax, width = 5, height = 5)
  ## Raw data
  ax <- ggplot_pca(
    simdata$X, simdata$metadata,
    col = 'batch', pch = 'class'
  )
  file <- sprintf("tmp/fig/microarray/magnitude/pca_raw-bal_%.02f.jpg", delta)
  ggsave(file, ax, width = 5, height = 5)
}

# No batch effects
crosstab <- matrix(40, 2, 1)
simdata <- simulate_microarray(
  10000, crosstab,
  delta = 0.1,
  gamma = 0.5,
  phi = 0.2,
  c = 10, d = 6,
  epsilon = 0.5,
  kappa = 0.2,
  dropout = FALSE
)
print(table(simdata$metadata$class, simdata$metadata$batch))
batch <- rep("1", 80)
batch[21:60] <- "2"
simdata$metadata$batch <- as.factor(batch)
print(table(simdata$metadata$class, simdata$metadata$batch))
file <- "data/simulated/microarray/magnitude/bal-0.00.rds"
print(file)
saveRDS(simdata, file)

## Normalised data
ax <- 2 ^ simdata$X %>%
  scale_trimmed() %>%
  log2_transform() %>%
  ggplot_pca(simdata$metadata, col = 'batch', pch = 'class')
file <- "tmp/fig/microarray/magnitude/pca-bal_0.00.jpg"
ggsave(file, ax, width = 5, height = 5)

# Measure omega terms
total_vars <- numeric()
deltas <- seq(0, 1, 0.1)
for (delta in deltas) {
  file <- sprintf("data/simulated/microarray/magnitude/bal-%.02f.rds", delta)
  print(file)
  simdata <- readRDS(file)
  total_var <- sum(apply(simdata$batch.terms, 1, var))
  total_vars <- c(total_vars, total_var)
}
batch_var <- data.frame(deltas, total_vars)

ax <- ggplot(batch_var, aes(x = deltas, y = total_vars)) +
  geom_line()
file <- "tmp/fig/microarray/magnitude/delta-total_var.jpg"
ggsave(file, ax)

## Imbalanced
nclass <- 2
nbatch <- 2
crosstab <- matrix(c(10, 30, 30, 10), nclass, nbatch)
print(crosstab)

for (delta in seq(0.1, 1.0, 0.1)) {
  simdata <- simulate_microarray(
    10000, crosstab,
    delta = delta,
    gamma = 0.5,
    phi = 0.2,
    c = 10, d = 6,
    epsilon = 0.5,
    kappa = 0.2,
    dropout = FALSE
  )
  file <- sprintf("data/simulated/microarray/magnitude/imbal-%.01f.rds", delta)
  print(file)
  saveRDS(simdata, file)
  # Plot
  ## Normalised data
  ax <- 2 ^ simdata$X %>%
    scale_trimmed() %>%
    log2_transform() %>%
    ggplot_pca(simdata$metadata, col = 'batch', pch = 'class')
  file <- sprintf("tmp/fig/microarray/magnitude/pca-imbal_%.01f.jpg", delta)
  ggsave(file, ax, width = 5, height = 5)
}

# No batch effects
crosstab <- matrix(40, 2, 1)
simdata <- simulate_microarray(
  10000, crosstab,
  delta = 0.1,
  gamma = 0.5,
  phi = 0.2,
  c = 10, d = 6,
  epsilon = 0.5,
  kappa = 0.2,
  dropout = FALSE
)
print(table(simdata$metadata$class, simdata$metadata$batch))
batch <- rep("1", 80)
batch[11:50] <- "2"
simdata$metadata$batch <- as.factor(batch)
print(table(simdata$metadata$class, simdata$metadata$batch))
file <- "data/simulated/microarray/magnitude/imbal-0.0.rds"
print(file)
saveRDS(simdata, file)

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

# Without.Balanced
crosstab <- matrix(160, 2, 1)
data <- simulate_microarray(
  10000, crosstab, delta = 0.7,
  phi = 0.3, zeta = 0.7
)
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

# Simulate: Different magnitudes
crosstab <- matrix(20, 2, 2)
data <- simulate_microarray(
  10000, crosstab, delta = 0.7,
  phi = 0.3, zeta = 0.7
)


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
