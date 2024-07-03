library(dplyr)
library(scater)
library(ggplot2)
theme_set(theme_bw(base_size = 7))
library(cowplot)
library(RColorBrewer)
src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


# ## Compile results: batch_scale v.s. values
# results_bal <- results_imbal <- list()
# niter <- 5
# batch_scales <- c(seq(0.00, 0.09, 0.01), seq(0.1, 0.3, 0.04))
# for (batch_scale in batch_scales) {
#   for (i in 1:niter) {
#     file <- sprintf("tmp/scrna/splat/results-%.2f-%02d.rds", batch_scale, i)
#     results1 <- readRDS(file)
#     file <- sprintf("tmp/scrna/splat/gpca_pvca-%.2f-%02d.rds", batch_scale, i)
#     results2 <- readRDS(file)
#     # Balanced: hvp, ... 
#     obj <- results1[[1]]
#     k <- obj$kbet$params$k0
#     hvp <- obj$rvp$percent.batch
#     cms <- mean(obj$cms$cms)
#     kbet <- obj$kbet$summary$kBET.observed[1]
#     lisi <- mean(obj$lisi$Batch)
#     # Balanced: gpca, pvca
#     obj2 <- results2[[1]]
#     gpca <- obj2$gpca$delta
#     pvca <- obj2$pvca$dat[3]
#     values_bal <- data.frame(
#       batch_scale = batch_scale,
#       batch_var = batch_scale ^ 2,
#       batch_class = "balanced",
#       metric = c("HVP", "CMS", "kBET", "LISI", "gPCA", "PVCA"),
#       value = c(hvp, cms, kbet, lisi, gpca, pvca)
#     )
#     results_bal <- append(results_bal, list(values_bal))
#     # Imbalanced: hvp, ... 
#     obj <- results1[[2]]
#     k <- obj$kbet$params$k0
#     hvp <- obj$rvp$percent.batch
#     cms <- mean(obj$cms$cms)
#     kbet <- obj$kbet$summary$kBET.observed[1]
#     lisi <- mean(obj$lisi$Batch)
#     # Imbalanced: gpca, pvca 
#     obj2 <- results2[[2]]
#     gpca <- obj2$gpca$delta
#     pvca <- obj2$pvca$dat[3]
#     gpca <- obj2$gpca$delta
#     pvca <- obj2$pvca$dat[3]
#     values_imbal <- data.frame(

#       batch_scale = batch_scale,
#       batch_var = batch_scale ^ 2,
#       batch_class = "imbalanced",
#       metric = c("HVP", "CMS", "kBET", "LISI", "gPCA", "PVCA"),
#       value = c(hvp, cms, kbet, lisi, gpca, pvca)
#     )
#     results_imbal <- append(results_imbal, list(values_imbal))
#   }
# }
# metrics_bal <- do.call(rbind, results_bal)
# metrics_imbal <- do.call(rbind, results_imbal)
# metrics <- rbind(metrics_bal, metrics_imbal)
# 
# file <- "tmp/scrna/splat/splat-metrics.csv"
# write.csv(metrics, file, row.names = FALSE)


## Plot: Metrics
file <- "tmp/scrna/splat/splat-metrics.csv"
metrics <- read.csv(file, stringsAsFactors = TRUE)

splat_summary <- metrics %>%
  group_by(batch_var, metric, batch_class) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )
metric_order <- c("HVP", "kBET", "CMS", "LISI", "gPCA", "PVCA")
splat_summary$metric <- factor(splat_summary$metric, levels = metric_order)
levels(splat_summary$batch_class) <- c(
  "Batch-class balanced",
  "Batch-class imbalanced"
)

GGCOLS <- ggplot_palette(6)
GGCOLS[6] <- "#BF80FF"
names(GGCOLS) <- c("HVP", "gPCA", "PVCA", "CMS", "kBET", "LISI")
ax <- splat_summary %>%
  subset(metric %in% c("HVP", "CMS", "kBET", "LISI")) %>%
  ggplot(aes(
    x = batch_var, y = mean,
    colour = factor(metric, levels = metric_order)
  )) +
  facet_grid(
   factor(metric, levels = metric_order) ~ batch_class, scales = "free_y"
  ) +
  geom_line(linewidth = 0.4, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.0017, linewidth = 0.2,
    color = "grey55"
  ) +
  scale_y_reverse() +
  scale_color_manual(values = GGCOLS) +
  labs(x = "Batch scale", y = "Metric value")
file <- "tmp/fig/splat/splat-rev.pdf"
ggsave(file, ax, width = 4.5, height = 5)

# Supplementary: gPCA, PVCA
ax <- splat_summary %>%
  subset(metric %in% c("gPCA", "PVCA")) %>%
  ggplot(aes(
    x = batch_var, y = mean,
    colour = factor(metric, levels = metric_order)
  )) +
  facet_grid(
   factor(metric, levels = metric_order) ~ batch_class, scales = "free_y"
  ) +
  geom_line(linewidth = 0.4, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.001, linewidth = 0.3,
    color = "grey55"
  ) +
  scale_color_manual(values = GGCOLS) +
  labs(x = "Batch scale", y = "Metric value")
file <- "tmp/fig/splat/splat-suppl.pdf"
ggsave(file, ax, width = 4.5, height = 2.8)

# Compile: Variance of data set and observed batch effects variance (variance
# of log batch factors matrix)
batch_vars <- total_vars_bal <- total_vars_imbal <- numeric()
niter <- 5
batch_scales <- c(seq(0, 0.09, 0.01), seq(0.1, 0.3, 0.04))
for (batch_scale in batch_scales) {
  for (i in 1:niter) {
    file1 <- sprintf(
      "data/simulated/scrna/loc0/splat-scale_%.2f-%02d.rds",
      batch_scale, i
    )
    file2 <- sprintf(
      "data/simulated/scrna/loc0/datasets-%.2f-%02d.rds",
      batch_scale, i
    )
    print(file1)
    splat <- readRDS(file1)
    if ("BatchFacBatch1" %in% colnames(rowData(splat))) {
      batch1_fac <- rowData(splat)$BatchFacBatch1
      batch2_fac <- rowData(splat)$BatchFacBatch2
      batch_fac <- cbind(
        replicate(2000, batch1_fac), replicate(2000, batch2_fac)
      )
      batch_var <- sum(rowVars(log2(batch_fac)))
    } else {
      batch_var <- 0
    }
    batch_vars <- c(batch_vars, batch_var)
    # Total variance of logcounts 
    datasets <- readRDS(file2)
    total_var_bal <- sum(rowVars(logcounts(datasets$bal)))
    total_var_imbal <- sum(rowVars(logcounts(datasets$imbal)))
    total_vars_bal <- c(total_vars_bal, total_var_bal)
    total_vars_imbal <- c(total_vars_imbal, total_var_imbal)
  }
}


# splat_vars <- data.frame(
#   batch_scale = rep(rep(batch_scales, each = 5), 2),
#   rep = rep(1:5, 32),
#   batch_class = rep(c("balanced", "imbalanced"), each = 80),
#   observed_var = rep(batch_vars, 2),
#   data_var = c(total_vars_bal, total_vars_imbal)
# )
# file <- "tmp/scrna/splat/splat-var.csv"
# write.csv(splat_vars, file, row.names = FALSE)

## Plot: Estimated v.s. observed batch effects variance
file <- "tmp/scrna/splat/splat-var.csv"
splat_vars <- read.csv(file, stringsAsFactors = TRUE)

file <- "tmp/scrna/splat/splat-metrics-old.csv"
splat_metrics <- read.csv(file)

splat_hvp <- subset(splat_metrics, metric == "RVP")
splat_vars$hvp <- splat_hvp$value
splat_vars$estimated_var <- splat_vars$hvp * splat_vars$data_var
levels(splat_vars$batch_class) <- c(
  "Batch-class balanced", "Batch-class imbalanced"
)

ax <- ggplot(splat_vars) +
  facet_wrap(~batch_class) +
  geom_point(
    aes(x = observed_var, y = estimated_var, col = batch_scale ^ 2),
    pch = 1, cex = 1
  ) +
  geom_abline(slope = 1, col = "grey", lty = "dashed") +
  labs(
    x = "Observed batch effects variance",
    y = "Estimated batch effects variance",
    col = "Batch scale"
  ) +
  annotate(
    geom = "text", x = 225, y = 265,
    label = "Ideal estimate",
    color = "darkgray", cex = 1.8, angle = 54
  ) +
  theme(
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 6)
  )
file <- "tmp/fig/splat/hvp-vars.pdf"
ggsave(file, ax, height = 2.2, width = 7.4)
