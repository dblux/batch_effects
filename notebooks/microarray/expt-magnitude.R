library(ggplot2)
library(RColorBrewer)
library(tidyr)
theme_set(theme_bw(base_size = 7))
src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


# results_bal <- results_imbal <- list()
# for (delta in seq(0, 1, 0.1)) {
#   file1a <- sprintf("data/simulated/microarray/magnitude/bal-%.01f.rds", delta)
#   file1b <- sprintf("tmp/microarray/magnitude/results-bal_%.01f.rds", delta)
#   print(file1a)
#   print(file1b)
# 
#   simdata <- readRDS(file1a)
#   batch_var <- sum(apply(simdata$batch.terms, 1, var))
#   X_var <- sum(apply(simdata$X, 1, var))
#   if (ncol(simdata$batch.logfc) == 2) {
#     expected_omega <- cbind(
#       matrix(simdata$batch.logfc[, 1], nrow = 10000, ncol = 40),
#       matrix(simdata$batch.logfc[, 2], nrow = 10000, ncol = 40)
#     )
#   } else if (ncol(simdata$batch.logfc) == 1) {
#     expected_omega <- matrix(simdata$batch.logfc, nrow = 10000, ncol = 80)
#   }
#   expected_batch_var <- sum(apply(expected_omega, 1, var))
# 
#   results <- readRDS(file1b)
#   rvp <- results$rvp$percent.batch
#   gpca <- results$gpca$delta
#   pvca <- results$pvca$dat[3]
#   bal_metrics <- c(delta, expected_batch_var, batch_var, X_var, rvp, gpca, pvca) 
#   names(bal_metrics) <- c(
#     "delta", "expected_batch_var", "batch_var", "X_var", "rvp", "gpca", "pvca"
#   )
#   results_bal <- append(results_bal, list(bal_metrics))
# 
#   file2a <- sprintf("data/simulated/microarray/magnitude/imbal-%.01f.rds", delta)
#   file2b <- sprintf("tmp/microarray/magnitude/results-imbal_%.01f.rds", delta)
#   print(file2a)
#   print(file2b)
# 
#   simdata <- readRDS(file2a)
#   batch_var <- sum(apply(simdata$batch.terms, 1, var))
#   X_var <- sum(apply(simdata$X, 1, var))
#   if (ncol(simdata$batch.logfc) == 2) {
#     expected_omega <- cbind(
#       matrix(simdata$batch.logfc[, 1], nrow = 10000, ncol = 40),
#       matrix(simdata$batch.logfc[, 2], nrow = 10000, ncol = 40)
#     )
#   } else if (ncol(simdata$batch.logfc) == 1) {
#     expected_omega <- matrix(simdata$batch.logfc, nrow = 10000, ncol = 80)
#   }
#   expected_batch_var <- sum(apply(expected_omega, 1, var))
# 
#   results <- readRDS(file2b)
#   rvp <- results$rvp$percent.batch
#   gpca <- results$gpca$delta
#   pvca <- results$pvca$dat[3]
#   imbal_metrics <- c(delta, expected_batch_var, batch_var, X_var, rvp, gpca, pvca) 
#   names(imbal_metrics) <- c(
#     "delta", "expected_batch_var", "batch_var", "X_var", "rvp", "gpca", "pvca"
#   )
#   results_imbal <- append(results_imbal, list(imbal_metrics))
# }
# 
# # Compile results
# short_bal <- data.frame(t(data.frame(results_bal)), row.names = NULL)
# file <- "tmp/microarray/magnitude/short_bal.csv"
# write.csv(short_bal, file, row.names = FALSE)
# 
# short_imbal <- data.frame(t(data.frame(results_imbal)), row.names = NULL)
# file <- "tmp/microarray/magnitude/short_imbal.csv"
# write.csv(short_imbal, file, row.names = FALSE)

# Plots
cex <- 0.5
metric_ord <- c("RVP", "gPCA", "PVCA", "CMS", "kBET", "LISI")
metric_cols <- ggplot_palette(6)
metric_cols[6] <- "#BF80FF"
names(metric_cols) <- metric_ord
balance_palette <- brewer.pal(8, "Set1")
theoretical_lab <- "Theoretical batch effects variance"
observed_lab <- "Observed batch effects variance"

file <- "tmp/microarray/magnitude/short_bal.csv"
short_bal <- read.csv(file)
file <- "tmp/microarray/magnitude/short_imbal.csv"
short_imbal <- read.csv(file)

ax <- ggplot(short_bal, aes(x = expected_batch_var, y = batch_var)) +
  geom_line() +
  geom_point(cex = cex) +
  labs(
    title = "Batch-class balanced",
    x = theoretical_lab,
    y = observed_lab
  )
file <- "tmp/fig/microarray/magnitude/bal-observed_theoretical.jpg"
ggsave(file, ax, width = 2.1, height = 2.1)

ax <- ggplot(short_imbal, aes(x = expected_batch_var, y = batch_var)) +
  geom_line() +
  geom_point(cex = cex) +
  labs(
    title = "Batch-class imbalanced",
    x = theoretical_lab,
    y = observed_lab
  )
file <- "tmp/fig/microarray/magnitude/imbal-observed_theoretical.jpg"
ggsave(file, ax, width = 2.1, height = 2.1)

balance_cols <- balance_palette[c(2, 8)]
names(balance_cols) <- c("Balanced", "Imbalanced")
ax <- ggplot() +
  geom_line(
    mapping = aes(x = expected_batch_var, y = batch_var, color = "Balanced"),
    data = short_bal, linetype = "twodash"
  ) +
  geom_point(
    mapping = aes(x = expected_batch_var, y = batch_var, color = "Balanced"),
    data = short_bal, cex = 1, alpha = 0.7
  ) +
  geom_line(
    mapping = aes(x = expected_batch_var, y = batch_var, color = "Imbalanced"),
    data = short_imbal, linetype = "dashed"
  ) +
  geom_point(
    mapping = aes(x = expected_batch_var, y = batch_var, color = "Imbalanced"),
    data = short_imbal, cex = 1, alpha = 0.7
  ) +
  scale_color_manual(values = balance_cols) +
  labs(
    x = theoretical_lab,
    y = observed_lab
  ) +
  theme(
    legend.position = c(0.01, 1.04),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.key = element_rect(fill='transparent'),
    legend.key.size = unit(4, "mm"),
    legend.background = element_rect(fill='transparent'),
    legend.title = element_blank()
  )
file <- "tmp/fig/microarray/magnitude/observed_theoretical.jpg"
ggsave(file, ax, width = 2.1, height = 2.1)

# # Consolidate results
# bal_var <- data.frame(
#   expected_batch_var = short_bal$expected_batch_var,
#   rvp = short_bal$rvp * short_bal$X_var,
#   gpca = short_bal$gpca * short_bal$X_var,
#   pvca = short_bal$pvca * short_bal$X_var
# )
# imbal_var <- data.frame(
#   expected_batch_var = short_imbal$expected_batch_var,
#   rvp = short_imbal$rvp * short_imbal$X_var,
#   gpca = short_imbal$gpca * short_imbal$X_var,
#   pvca = short_imbal$pvca * short_imbal$X_var
# )
# long_bal <- gather(bal_var, "metric", "value", -expected_batch_var)
# long_imbal <- gather(imbal_var, "metric", "value", -expected_batch_var)
# metrics <- rbind(long_bal, long_imbal)
# metrics$batch_class <- c(
#   rep("balanced", nrow(long_bal)),
#   rep("imbalanced", nrow(long_imbal))
# )
# file <- "tmp/microarray/magnitude/metrics-expected_var.csv"
# write.csv(metrics, file, row.names = FALSE)

# Plot: Metrics - Simulated microarray
file <- "tmp/microarray/magnitude/metrics-expected_var.csv"
metrics <- read.csv(file)
metrics$Metric <- factor(metrics$Metric, levels = c("gPCA", "PVCA", "RVP"))

imbalance_labs <- c("Batch-class balanced", "Batch-class imbalanced")
names(imbalance_labs) <- c("Balanced", "Imbalanced")
ax <- ggplot(metrics, aes(x = Expected.Batch.Var, y = Value, col = Metric)) +
  facet_wrap(
    ~ Batch.Class, nrow = 1, scales = "fixed",
    labeller = labeller(Batch.Class = imbalance_labs)
  ) +
  geom_line() +
  geom_point(cex = cex) +
  geom_abline(
    aes(slope = 1, intercept = 0), 
    color = "darkgray", linetype = "longdash"
  ) +
  labs(
    x = theoretical_lab,
    y = "Estimated batch effects variance" 
  ) +
  scale_color_manual(values = metric_cols) +
  annotate(
    geom = "text", x = 4200, y = 5000,
    label = "Perfect estimate",
    color = "darkgray", cex = 2, angle = 20
  ) +
  theme(
    legend.position = c(0.01, 1.04),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.key = element_rect(fill='transparent'),
    legend.key.size = unit(4, "mm"),
    legend.background = element_rect(fill='transparent'),
    legend.title = element_blank()
  )
file <- "tmp/fig/microarray/magnitude/metrics-expected_var.jpg"
ggsave(file, ax, width = 4, height = 2)

# Plot: PCA
batch_cols <- brewer.pal(7, "Dark2")[2:4]
class_cols <- brewer.pal(7, "Dark2")[5:7]

# Balanced
file <- "tmp/microarray/magnitude/short_bal.csv"
short_bal <- read.csv(file)
show.legend <- FALSE
width <- 1.5
batch_sizes <- seq(0, 1, 0.1)
for (i in seq_len(length(batch_sizes))) {
  if (i == 11) {
    show.legend <- TRUE
    width <- 1.9
  }
  batch_size <- batch_sizes[i]
  file1a <- sprintf(
    "data/simulated/microarray/magnitude/bal-%.01f.rds", batch_size
  )
  bvar_title <- sprintf(
    "Variance = %.00f",
    short_bal$expected_batch_var[i]
  )
  simdata <- readRDS(file1a)
  colnames(simdata$metadata) <- c("Class", "Batch")
  ax <- 2 ^ simdata$X %>%
    scale_trimmed() %>%
    log2_transform() %>%
    ggplot_pca(
      simdata$metadata, col = "Batch", pch = "Class",
      cex = 1.5, alpha = 0.7, show.legend = show.legend, plot.axis = FALSE
    )
  ax <- ax +
    labs(title = bvar_title) +
    theme(
      title = element_text(size = 6),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 5),
      axis.title.y = element_text(size = 5),
      legend.key.size = unit(4, "mm")
    ) +
    scale_color_manual(values = batch_cols)
  file <- sprintf(
    "tmp/fig/microarray/magnitude/pca-microarray_bal_%.01f.jpg", batch_size
  )
  ggsave(file, ax, width = width, height = 1.5)
  print(file)
}

# Imbalanced
file <- "tmp/microarray/magnitude/short_imbal.csv"
short_imbal <- read.csv(file)
show.legend <- FALSE
width <- 1.5
for (i in seq_len(length(batch_sizes))) {
  if (i == 11) {
    show.legend <- TRUE
    width <- 1.9
  }
  batch_size <- batch_sizes[i]
  file1a <- sprintf(
    "data/simulated/microarray/magnitude/imbal-%.01f.rds", batch_size
  )
  bvar_title <- sprintf(
    "Variance = %.00f",
    short_imbal$expected_batch_var[i]
  )
  simdata <- readRDS(file1a)
  colnames(simdata$metadata) <- c("Class", "Batch")
  ax <- 2 ^ simdata$X %>%
    scale_trimmed() %>%
    log2_transform() %>%
    ggplot_pca(
      simdata$metadata, col = "Batch", pch = "Class",
      cex = 1.5, alpha = 0.7, show.legend = show.legend, plot.axis = FALSE
    )
  ax <- ax +
    labs(title = bvar_title) +
    theme(
      title = element_text(size = 6),
      axis.title.x = element_text(size = 5),
      axis.title.y = element_text(size = 5),
      legend.key.size = unit(4, "mm")
    ) +
    scale_color_manual(values = batch_cols)
  file <- sprintf(
    "tmp/fig/microarray/magnitude/pca-microarray_imbal_%.01f.jpg", batch_size
  )
  ggsave(file, ax, width = width, height = 1.5)
  print(file)
}
