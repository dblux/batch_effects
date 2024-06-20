library(splatter)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
theme_set(theme_bw(base_size = 6))
source("R/plot.R")
source("R/simulate.R")
source("R/HVP.R")


# Microarray simulation
set.seed(1)
crosstab <- matrix(5, 2, 2)
simdata <- simulate_microarray(
  1000, crosstab,
  delta = 1.0,
  gamma = 0.5,
  phi = 0.2,
  c = 10, d = 6,
  epsilon = 0.5,
  kappa = 0.2,
  dropout = FALSE
)
names(simdata)

log_psi <- as.matrix(simdata$log_psi)
log_rho <- abs(simdata$class.logfc[, 2, drop = FALSE])
log_alpha <- abs(t(simdata$log_alpha))
log_alpha_mat <- abs(t(replicate(30, simdata$log_alpha)))
Z <- simdata$Z
# Make terms from one class negative so that it takes on a different color
class_col <- ifelse(simdata$metadata$class == "B", -1, 1)
Z_col <- Z %*% diag(class_col)
log_beta1 <- abs(simdata$batch.logfc[, 1, drop = FALSE])
log_beta2 <- abs(simdata$batch.logfc[, 2, drop = FALSE])
log_omega <- abs(simdata$batch.terms)
# Make terms from one batch negative so that it takes on a different color
batch_col <- ifelse(simdata$metadata$batch == 2, 1, -1)
log_omega_col <- log_omega %*% diag(batch_col)
X <- simdata$X

# Heatmap
pheatmap(
  log_psi,
  color = brewer.pal(9, "Greys"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/log_psi.pdf"
)
pheatmap(
  log_rho,
  color = brewer.pal(9, "Reds"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/log_rho.pdf"
)

rho_breaks <- seq(-11.4, 11.4, length.out = 17)
rho_palette <- colorRampPalette(brewer.pal(
  n = 9, name = "RdGy"
))(length(rho_breaks))
pheatmap(
  Z_col,
  color = rho_palette, breaks = rho_breaks, border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/Z.pdf"
)

pheatmap(
  log_alpha,
  color = brewer.pal(9, "Greens"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/log_alpha.pdf"
)
pheatmap(
  log_alpha_mat,
  color = brewer.pal(9, "Greens"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/log_alpha_mat.pdf"
)

pheatmap(
  log_beta1,
  color = brewer.pal(9, "Oranges"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/log_beta1.pdf"
)
pheatmap(
  log_beta2,
  color = brewer.pal(9, "Purples"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/log_beta2.pdf"
)

beta_breaks <- seq(-3.9, 3.9, length.out = 17)
beta_palette <- colorRampPalette(brewer.pal(
  n = 9, name = "PuOr"
))(length(beta_breaks))
pheatmap(
  log_omega_col,
  color = beta_palette, breaks = beta_breaks, border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/log_omega.pdf"
)

pheatmap(
  X,
  color = brewer.pal(9, "PuRd"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/X.pdf"
)

# Splat - Methodology:
# Simulate batch log-normal terms for each batch
# Simulate batch cell means (BatchCellMeans)
# Simulate DE genes
# Simulate different library sizes for each sample (BaseCellMeans)
# Enforce mean-variance trend (CellMeans)
# Sample counts from Poisson distribution (TrueCounts)
# Apply dropouts (counts)
# Data is normalised and log-transformed before measuring for batch effects

# Simulate data
set.seed(0)
batch_scale <- 0.3
params <- newSplatParams()
params <- setParams(
  params,
  nGenes = 30,
  batchCells = c(20, 20),
  group.prob = c(0.5, 0.5),
  batch.facLoc = 0, # log-normal
  batch.facScale = batch_scale # log-normal
)
splat <- splatSimulate(params, method = "groups")

# Heatmap: 30 x 20
sid <- c(
  c(1,2,8,10,11),
  3:7,
  c(22:25,27),
  c(26,28:30,32)
)
metadata <- colData(splat)[sid, ]
table(metadata$Batch, metadata$Group)

# Information
assays(splat)
rowData(splat)

# BaseGeneMean * OutlierFactor = GeneMean
bgm <- rowData(splat)[, "BaseGeneMean"]
of <- rowData(splat)[, "OutlierFactor"]
gm <- rowData(splat)[, "GeneMean"]
gm_mat <- replicate(20, gm)
bf1 <- rowData(splat)[, "BatchFacBatch1"]
bf2 <- rowData(splat)[, "BatchFacBatch2"]
bf_mat <- cbind(replicate(10, -bf1), replicate(10, bf2))
cf1 <- rowData(splat)[, "DEFacGroup1"]
cf2 <- rowData(splat)[, "DEFacGroup2"]
cf_mat <- cbind(
  replicate(5, cf1), replicate(5, -cf2),
  replicate(5, cf1), replicate(5, -cf2)
)
libsize <- colData(splat)$ExpLibSize[sid]
libsize_mat <- t(replicate(30, libsize))
bcm_mat <- assays(splat)$BaseCellMeans[, sid]
counts <- assays(splat)$counts[, sid]

pheatmap(
  bgm,
  color = brewer.pal(9, "Greys"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/bgm.pdf"
)
pheatmap(
  of,
  color = brewer.pal(9, "Greys"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/of.pdf"
)
pheatmap(
  gm_mat,
  color = brewer.pal(9, "Greys"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/gm_mat.pdf"
)

pheatmap(
  bf1,
  color = brewer.pal(9, "Oranges"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/bf1.pdf"
)
pheatmap(
  bf2,
  color = brewer.pal(9, "Purples"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/bf2.pdf"
)
pheatmap(
  bf_mat,
  color = brewer.pal(9, "PuOr"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/bf_mat.pdf"
)

pheatmap(
  cf1,
  color = brewer.pal(9, "Blues"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/cf1.pdf"
)
pheatmap(
  cf2,
  color = brewer.pal(9, "Reds"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/cf2.pdf"
)
pheatmap(
  cf_mat,
  color = brewer.pal(9, "RdBu"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/cf_mat.pdf"
)

pheatmap(
  t(libsize),
  color = brewer.pal(9, "Greens"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/libsize.pdf"
)
pheatmap(
  libsize_mat, 
  color = brewer.pal(9, "Greens"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/libsize_mat.pdf"
)

pheatmap(
  bcm_mat,
  color = brewer.pal(9, "PuRd"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/bcm_mat.pdf"
)

pheatmap(
  counts,
  color = brewer.pal(9, "PuRd"), border_color = NA,
  cellwidth = 3, cellheight = 3,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  scale = "none", legend = FALSE,
  file = "tmp/fig/illustrations/counts.pdf"
)


# RVP: Jitter plot
set.seed(1)
crosstab <- matrix(10, 2, 2)
simdata <- simulate_microarray(
  30, crosstab,
  delta = 1.0,
  gamma = 0.5,
  phi = 0.2,
  c = 10, d = 6,
  epsilon = 0.5,
  kappa = 0.2,
  dropout = FALSE
)
feature <- data.frame(
  value = simdata$X[6, ],
  batch = simdata$metadata$batch,
  class = simdata$metadata$class
)

batch_cols <- brewer.pal(7, "Dark2")[2:4]
class_cols <- brewer.pal(7, "Dark2")[5:7]
ax <- ggplot(
  feature,
  aes(x = class, y = value, col = batch, pch = class)
) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),
    cex = 1, show.legend = TRUE
  ) +
  geom_hline(
    yintercept = mean(feature$value),
    lty = "longdash", linewidth = 0.15
  ) +
  scale_color_manual(values = batch_cols) +
  labs(
    x = "Class", y = "Expression value", col = "Batch", pch = "Class"
  ) +
  stat_summary(
    fun = mean, geom = "point", pch = 5, cex = 1
  ) +
  stat_summary(
    fun = mean, geom = "crossbar",
    col = "darkgrey", lty = "longdash", linewidth = 0.08
  ) +
  theme(
    legend.key.size = unit(4, "mm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
file <- "tmp/fig/illustrations/feature.pdf"
ggsave(file, ax, width = 3, height = 1.5)


# Plot HVP sum of squares
res <- HVP(simdata$X, simdata$metadata$batch, simdata$metadata$class)

ax <- plot_hvp(res)
file <- "tmp/fig/hvp.pdf"
ggsave(file, ax, width = 4.5, height = 1.5)
