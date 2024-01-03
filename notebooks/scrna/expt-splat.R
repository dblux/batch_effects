library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 7))
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(scater)
library(CellMixS)
library(kBET)
library(lisi)

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


# # Results
# results_bal <- results_imbal <- list()
# # Approximate magnitude from batch log-normal factors
# # Assume: Number of samples in batch 1 and 2 are equal
# niter <- 5
# batch_scales <- c(seq(0.00, 0.09, 0.01), seq(0.1, 0.3, 0.04))
# for (batch_scale in batch_scales) {
#   for (i in 1:niter) {
#     file1 <- sprintf(
#       "data/simulated/scrna/loc0/splat-scale_%.2f-%02d.rds",
#       batch_scale, i
#     )
#     print(file1)
#     splat <- readRDS(file1)
#     fac_scale <- splat@metadata$Params@batch.facScale[1]
#     splat_means <- rowData(splat)
#     magnitude <- lpnorm(
#       splat_means$BatchFacBatch1 - splat_means$BatchFacBatch2
#     )
# 
#     file2 <- sprintf(
#       "tmp/scrna/splat/results-%.2f-%02d.rds",
#       batch_scale, i
#     )
#     results <- readRDS(file2)
# 
#     obj <- results[[1]]
#     k <- obj$kbet$params$k0
#     rvp <- obj$rvp$percent.batch
#     cms <- mean(obj$cms$cms)
#     kbet <- obj$kbet$summary$kBET.observed[1]
#     lisi <- mean(obj$lisi$Batch)
#     bal_row <- c(batch_scale, magnitude, rvp, cms, kbet, lisi)
#     names(bal_row) <- c(
#       "batch_scale", "magnitude", "rvp", "cms", "kbet", "lisi"
#     )
#     results_bal <- append(results_bal, list(bal_row))
# 
#     obj <- results[[2]]
#     k <- obj$kbet$params$k0
#     rvp <- obj$rvp$percent.batch
#     cms <- mean(obj$cms$cms)
#     kbet <- obj$kbet$summary$kBET.observed[1]
#     lisi <- mean(obj$lisi$Batch)
#     imbal_row <- c(batch_scale, magnitude, rvp, cms, kbet, lisi)
#     names(imbal_row) <- c(
#       "batch_scale", "magnitude", "rvp", "cms", "kbet", "lisi"
#     )
#     results_imbal <- append(results_imbal, list(imbal_row))
#   }
# }
# 
# # Compile results
# short_bal <- data.frame(t(data.frame(results_bal)), row.names = NULL)
# long_bal <- gather(
#   short_bal, key = "metric", value = "value",
#   -batch_scale, -magnitude
# )
# short_imbal <- data.frame(t(data.frame(results_imbal)), row.names = NULL)
# long_imbal <- gather(
#   short_imbal, key = "metric", value = "value",
#   -batch_scale, -magnitude
# )
# 
# splat_metrics <- rbind(long_bal, long_imbal)
# splat_metrics$batch_class <- c(
#   rep("balanced", nrow(long_bal)),
#   rep("imbalanced", nrow(long_imbal))
# )
# file <- "tmp/scrna/splat/splat-metrics.csv"
# write.csv(splat_metrics, file, row.names = FALSE)

# Read saved results
file <- "tmp/scrna/splat/splat-metrics.csv"
splat_metrics <- read.csv(file)
splat_metrics$metric <- factor(
  splat_metrics$metric, levels = c("rvp", "kbet", "cms", "lisi")
)

splat_summary <- splat_metrics %>%
  group_by(batch_scale, metric, batch_class) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    magnitude = mean(magnitude),
  )

imbalance_labs <- c("Batch-class balanced", "Batch-class imbalanced")
names(imbalance_labs) <- c("balanced", "imbalanced")
metric_labs <- c("RVP", "CMS", "kBET", "LISI")
names(metric_labs) <- c("rvp", "cms", "kbet", "lisi")

METRIC_ORD <- c("rvp", "gpca", "pvca", "cms", "kbet", "lisi")
CEX <- 0.5
GGCOLS <- ggplot_palette(6)
GGCOLS[6] <- "#BF80FF"

names(GGCOLS) <- METRIC_ORD
ax <- ggplot(
  splat_summary,
  aes(x = batch_scale, y = mean, colour = metric)
) +
  facet_grid(
    metric ~ batch_class, scales = "free_y",
    labeller = labeller(metric = metric_labs, batch_class = imbalance_labs)
  ) +
  geom_line(show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.008, linewidth = 0.3, color = "grey55"
  ) +
  # scale_y_reverse() +
  scale_color_manual(values = GGCOLS) +
  labs(x = "Batch scale", y = "Metric value")
file <- "tmp/fig/splat/splat.jpg"
ggsave(file, ax, width = 4, height = 4)
