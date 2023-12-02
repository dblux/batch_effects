library(ggplot2)

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


# Simulate data set
npercond <- c(20, 20, 20, 20)
crosstab <- matrix(npercond, 2, 2)

for (delta in seq(0.2, 1.2, 0.1)) {
  print(delta)
  simdata <- simulate_microarray(
    8000, sum(npercond), crosstab,
    delta = delta, phi = 0.1, zeta = 1.5,
    dropout = FALSE
  )

  file1 <- sprintf('../data/simulated/small/balanced/bal-%.1f.pdf', delta)
  write.table(simdata$X, file1, quote = FALSE, sep = '\t', row.names = F)

  ax_pca <- ggplot_pca(simdata$X, simdata$metadata, col = 'batch', pch = 'class')
  file2 <- sprintf('~/Dropbox/tmp/pca_bal-%.1f.pdf', delta)
  ggsave(file2, ax_pca, width = 6, height = 4)

  ax_top_pc <- ggplot_top_pc(
    simdata$X, simdata$metadata,
    x_axis = 'batch', col = 'batch', pch = 'class'
  )
  file3 <- sprintf('~/Dropbox/tmp/top_pc_bal-%.1f.pdf', delta)
  ggsave(file3, ax_top_pc, width = 6, height = 4)

  file4 <- sprintf('~/Dropbox/tmp/log_beta_bal-%.1f.pdf', delta)
  pdf(file4, width = 6, height = 4)
  hist(simdata$log_beta[, 1])
  dev.off()
}

# Load simulated data
dir <- '../data/simulated/small/balanced'
files <- list.files(dir, full.names = T)
files <- files[c(9, 1:8, 10:11)]
print(files)
list_data <- lapply(files, read.table, sep = "\t", header = T)

ids <- sapply(files, function(x) substring(x, 34), USE.NAMES = F)
print(ids)
names(list_data) <- ids

# metadata
ncond <- 20
batch <- as.factor(rep(1:2, each = ncond * 2))
class <- rep(rep(LETTERS[1:2], each = ncond), 2)
metadata <- data.frame(batch, class, row.names = colnames(list_data[[1]]))

print(names(list_data))
X <- list_data[[2]]

# Dev: Data simulation
# Orthogonality between batch effects and class effects
#   - Features with both batch and class effects results in non-orthogonality
#   - Visualisation: Random unit vectors contained in a unit sphere
#   - By chance, class effect genes may overlap with genes with high batch effects, resulting in non-orthogonality

nclass <- 2
nbatch <- 2
npercond <- 20
# n_kg <- rep(npercond, nbatch * nclass)
crosstab <- matrix(npercond, nclass, nbatch)
print(crosstab)

delta <- 0.7
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

## Class log fold-change
class_fc <- as.vector(simdata$class.logfc)
nonzero_cfc <- class_fc[class_fc != 0]
head(sort(abs(nonzero_cfc), decreasing = T), 30)
head(sort(abs(nonzero_cfc), decreasing = F), 30)

# Plot
## Normalised data
ax_pca <- 2 ^ simdata$X %>%
  scale_trimmed() %>%
  log2_transform() %>%
  ggplot_pca(simdata$metadata, col = 'batch', pch = 'class')
ax_pca

file <- "tmp/pca-simulated_scaled.png"
ggsave(file, ax_pca, width = 5, height = 5)

i <- 1
feature <- data.frame(
  value = as.vector(data.matrix(X_scaled[i, ])),
  batch = simdata$metadata$batch,
  class = simdata$metadata$class
)

ax <- ggplot(feature) +
  geom_point(
    aes(x = batch, y = value, pch = class, col = batch),
    position = position_jitterdodge(jitter.width = .5),
    show.legend = FALSE
  ) +
  ylim(0, 15)
ax
file <- "tmp/feature-simulated_scaled.png"
ggsave(file, ax, width = 5, height = 5)

## Un-normalised data
ax_pca <- ggplot_pca(
  simdata$X, simdata$metadata,
  col = 'batch', pch = 'class'
)
ax_pca

file <- "tmp/pca-simulated_b10c5_delta20.png"
ggsave(file, ax_pca, width = 5, height = 5)

## Feature
i <- 0
i <- i + 1
feature <- data.frame(
  value = as.vector(data.matrix(simdata$X[i, ])),
  batch = simdata$metadata$batch,
  class = simdata$metadata$class
)

ax <- ggplot(feature) +
  geom_point(
    aes(x = batch, y = value, pch = class, col = batch),
    position = position_jitterdodge(jitter.width = .5),
    show.legend = FALSE
  ) +
  ylim(0, 15)
ax

file <- "tmp/feature-simulated.png"
ggsave(file, ax, width = 5, height = 5)


## Colsum
l1 <- colSums(simdata$X)
feature <- data.frame(
  value = l1,
  batch = simdata$metadata$batch,
  class = simdata$metadata$class
)

ax <- ggplot(feature) +
  geom_point(
    aes(x = batch, y = value, pch = class, col = batch),
    position = position_jitterdodge(jitter.width = .5),
    show.legend = FALSE
  )
ax

file <- "tmp/l1-simulated.png"
ggsave(file, ax, width = 5, height = 5)

## Histogram
hist(simdata$Z[1, ], breaks = 20)
