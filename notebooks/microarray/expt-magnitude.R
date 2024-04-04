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


# # Compile results
# results_bal <- results_imbal <- list()
# for (delta in seq(0, 1, 0.1)) {
#   # Balanced
#   file1a <- sprintf("data/simulated/microarray/magnitude/bal-%.01f.rds", delta)
#   file1b <- sprintf("tmp/microarray/magnitude/results-bal_%.01f.rds", delta)
#   print(file1a)
#   print(file1b)
#   simdata <- readRDS(file1a)
#   observed_batchvar <- sum(apply(simdata$batch.terms, 1, var))
#   X_norm <- log2_transform(scale_trimmed(2 ^ simdata$X)) # log-normalise
#   X_var <- sum(apply(X_norm, 1, var))
#   # Theoretical batch effects variance
#   if (ncol(simdata$batch.logfc) == 2) {
#     expected_omega <- cbind(
#       matrix(simdata$batch.logfc[, 1], nrow = 10000, ncol = 40),
#       matrix(simdata$batch.logfc[, 2], nrow = 10000, ncol = 40)
#     )
#   } else if (ncol(simdata$batch.logfc) == 1) {
#     expected_omega <- matrix(simdata$batch.logfc, nrow = 10000, ncol = 80)
#   }
#   theoretical_batchvar <- sum(apply(expected_omega, 1, var))
#   # Metric values
#   results <- readRDS(file1b)
#   rvp <- results$rvp$percent.batch
#   gpca <- results$gpca$delta
#   pvca <- results$pvca$dat[3]
#   bal_metrics <- c(
#     delta, theoretical_batchvar, observed_batchvar,
#     X_var, rvp, gpca, pvca
#   )
#   names(bal_metrics) <- c(
#     "delta", "theoretical_batchvar", "observed_batchvar",
#     "X_var", "rvp", "gpca", "pvca"
#   )
#   results_bal <- append(results_bal, list(bal_metrics))
#   # Imbalanced
#   file2a <- sprintf("data/simulated/microarray/magnitude/imbal-%.01f.rds", delta)
#   file2b <- sprintf("tmp/microarray/magnitude/results-imbal_%.01f.rds", delta)
#   print(file2a)
#   print(file2b)
#   simdata <- readRDS(file2a)
#   observed_batchvar <- sum(apply(simdata$batch.terms, 1, var))
#   X_norm <- log2_transform(scale_trimmed(2 ^ simdata$X)) # log-normalise
#   X_var <- sum(apply(X_norm, 1, var))
#   # Theoretical batch effects variance
#   if (ncol(simdata$batch.logfc) == 2) {
#     expected_omega <- cbind(
#       matrix(simdata$batch.logfc[, 1], nrow = 10000, ncol = 40),
#       matrix(simdata$batch.logfc[, 2], nrow = 10000, ncol = 40)
#     )
#   } else if (ncol(simdata$batch.logfc) == 1) {
#     expected_omega <- matrix(simdata$batch.logfc, nrow = 10000, ncol = 80)
#   }
#   theoretical_batchvar <- sum(apply(expected_omega, 1, var))
#   # Metric values
#   results <- readRDS(file2b)
#   rvp <- results$rvp$percent.batch
#   gpca <- results$gpca$delta
#   pvca <- results$pvca$dat[3]
#   imbal_metrics <- c(
#     delta, theoretical_batchvar, observed_batchvar,
#     X_var, rvp, gpca, pvca
#   )
#   names(imbal_metrics) <- c(
#     "delta", "theoretical_batchvar", "observed_batchvar",
#     "X_var", "rvp", "gpca", "pvca"
#   )
#   results_imbal <- append(results_imbal, list(imbal_metrics))
# }
# 
# short_bal <- data.frame(t(data.frame(results_bal)), row.names = NULL)
# file <- "tmp/microarray/magnitude/short_bal.csv"
# write.csv(short_bal, file, row.names = FALSE)
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

balance_cols <- balance_palette[c(2, 8)]
names(balance_cols) <- c("Balanced", "Imbalanced")
ax <- ggplot() +
  geom_line(
    aes(x = theoretical_batchvar, y = observed_batchvar, color = "Balanced"),
    data = short_bal, linetype = "twodash"
  ) +
  geom_point(
    aes(x = theoretical_batchvar, y = observed_batchvar, color = "Balanced"),
    data = short_bal, cex = 1, alpha = 0.7
  ) +
  geom_line(
    aes(x = theoretical_batchvar, y = observed_batchvar, color = "Imbalanced"),
    data = short_imbal, linetype = "dashed"
  ) +
  geom_point(
    aes(x = theoretical_batchvar, y = observed_batchvar, color = "Imbalanced"),
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
file <- "tmp/fig/microarray/magnitude/observed-theoretical.pdf"
ggsave(file, ax, width = 2.4, height = 2.4)

# # Consolidate results
# bal_var <- data.frame(
#   theoretical_batchvar = short_bal$theoretical_batchvar,
#   RVP = short_bal$rvp * short_bal$X_var,
#   gPCA = short_bal$gpca * short_bal$X_var,
#   PVCA = short_bal$pvca * short_bal$X_var
# )
# imbal_var <- data.frame(
#   theoretical_batchvar = short_imbal$theoretical_batchvar,
#   RVP = short_imbal$rvp * short_imbal$X_var,
#   gPCA = short_imbal$gpca * short_imbal$X_var,
#   PVCA = short_imbal$pvca * short_imbal$X_var
# )
# long_bal <- gather(bal_var, "metric", "value", -theoretical_batchvar)
# long_imbal <- gather(imbal_var, "metric", "value", -theoretical_batchvar)
# metrics <- rbind(long_bal, long_imbal)
# metrics$batch_class <- c(
#   rep("balanced", nrow(long_bal)),
#   rep("imbalanced", nrow(long_imbal))
# )
# file <- "tmp/microarray/magnitude/metrics-expected_var.csv"
# write.csv(metrics, file, row.names = FALSE)

# Plot: Metrics - Simulated microarray
file <- "tmp/microarray/magnitude/metrics-expected_var.csv"
batchvar <- read.csv(file)
batchvar$metric <- factor(batchvar$metric, levels = c("gPCA", "PVCA", "RVP"))

imbalance_labs <- c("Batch-class balanced", "Batch-class imbalanced")
names(imbalance_labs) <- c("balanced", "imbalanced")
ax <- ggplot(
  batchvar,
  aes(x = theoretical_batchvar, y = value, col = metric)
) +
  facet_wrap(
    ~ batch_class, nrow = 1, scales = "fixed",
    labeller = labeller(batch_class = imbalance_labs)
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
file <- "tmp/fig/microarray/magnitude/estimated-theoretical.pdf"
ggsave(file, ax, width = 4.8, height = 2.4)

# var_gpca_lab <- expression(paste("gPCA ", delta, " * ", S[bold(X)]^2))
# var_pvca_lab <- expression(paste("PVCA * ", S[bold(X)]^2))
# var_rvp_lab <- expression(paste("RVP * ", S[bold(X)]^2))
