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


# Dev: Evaluate batch effects
k <- 300
results <- lapply(
  splat_objs, eval_batch,
  "Batch", "Group",
  nperm = 0, k.cms = k, k0 = k, perplexity = k,
  ret.scores = FALSE
)
saveRDS(results, file)

# files <- character()
# nsample <- 5
# for (i in 1:nsample) {
#   file <- sprintf("tmp/scrna/splat/loc0-v1/results-0.00-%02d.rds", i)
#   files <- c(files, file)
# }
# batch_scales <- c(seq(0.01, 0.09, 0.01), seq(0.1, 0.3, 0.04))
# for (batch_scale in batch_scales) {
#   file <- sprintf("tmp/scrna/splat/loc0-v1/results-%.2f.rds", batch_scale)
#   files <- c(files, file)
# }
# print(files)


# Results
results_bal <- results_imbal <- list()
# Approximate magnitude from batch log-normal factors
# Assume: Number of samples in batch 1 and 2 are equal
niter <- 5
batch_scales <- c(seq(0.00, 0.09, 0.01), seq(0.1, 0.3, 0.04))
for (batch_scale in batch_scales) {
  for (i in 1:niter) {
    file1 <- sprintf(
      "data/simulated/scrna/loc0/splat-scale_%.2f-%02d.rds",
      batch_scale, i
    )
    print(file1)
    splat <- readRDS(file1)
    fac_scale <- splat@metadata$Params@batch.facScale[1]
    splat_means <- rowData(splat)
    magnitude <- lpnorm(
      splat_means$BatchFacBatch1 - splat_means$BatchFacBatch2
    )

    file2 <- sprintf(
      "tmp/scrna/splat/results-%.2f-%02d.rds",
      batch_scale, i
    )
    results <- readRDS(file2)

    obj <- results[[1]]
    k <- obj$kbet$params$k0
    rvp <- obj$rvp$percent.batch
    cms <- mean(obj$cms$cms)
    kbet <- obj$kbet$summary$kBET.observed[1]
    lisi <- mean(obj$lisi$Batch)
    bal_row <- c(batch_scale, magnitude, rvp, cms, kbet, lisi)
    names(bal_row) <- c(
      "batch_scale", "magnitude", "rvp", "cms", "kbet", "lisi"
    )
    results_bal <- append(results_bal, list(bal_row))

    obj <- results[[2]]
    k <- obj$kbet$params$k0
    rvp <- obj$rvp$percent.batch
    cms <- mean(obj$cms$cms)
    kbet <- obj$kbet$summary$kBET.observed[1]
    lisi <- mean(obj$lisi$Batch)
    imbal_row <- c(batch_scale, magnitude, rvp, cms, kbet, lisi)
    names(imbal_row) <- c(
      "batch_scale", "magnitude", "rvp", "cms", "kbet", "lisi"
    )
    results_imbal <- append(results_imbal, list(imbal_row))
  }
}

# Compile results
short_bal <- data.frame(t(data.frame(results_bal)), row.names = NULL)
long_bal <- gather(
  short_bal, key = "metric", value = "value",
  -batch_scale, -magnitude
)
short_imbal <- data.frame(t(data.frame(results_imbal)), row.names = NULL)
long_imbal <- gather(
  short_imbal, key = "metric", value = "value",
  -batch_scale, -magnitude
)

splat_metrics <- rbind(long_bal, long_imbal)
splat_metrics$batch_class <- c(
  rep("balanced", nrow(long_bal)),
  rep("imbalanced", nrow(long_imbal))
)
file <- "tmp/scrna/splat/splat-metrics.csv"
write.csv(splat_metrics, file, row.names = FALSE)

# Read saved results
file <- "tmp/scrna/splat/splat-metrics.csv"
splat_metrics <- read.csv(file)

ax <- splat_metrics %>%
  subset(subset = metric == "rvp" & batch_class == "balanced") %>%
  ggplot(aes(x = batch_scale, y = magnitude)) +
    geom_point()
file <- "tmp/fig/splat/scale-magnitude.jpg"
ggsave(file, ax, width = 7, height = 7)

splat_summary <- splat_metrics %>%
  group_by(batch_scale, metric, batch_class) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    magnitude = mean(magnitude),
  )

lw <- 0.3
ax <- ggplot(
  splat_summary,
  aes(x = batch_scale ^ 2, y = mean, colour = metric)
) +
  facet_grid(metric ~ batch_class, scales = "free_y") +
  geom_line(linewidth = lw) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.002, linewidth = lw, color = "grey55"
  ) +
  theme(text = element_text(size = 8))
file <- "tmp/fig/splat/splat.jpg"
ggsave(file, ax, width = 7, height = 7)


# # Plots
# file <- "tmp/fig/splat/heatmap-bal.png"
# pheatmap(
#   short1_bal,
#   filename = file,
#   scale = "row",
#   color = colorRampPalette(brewer.pal(9, "Blues"))(100),
#   show_colnames = TRUE,
#   cluster_rows = FALSE, cluster_cols = FALSE
# )
# file <- "tmp/fig/splat/heatmap-imbal.png"
# pheatmap(
#   short1_imbal,
#   filename = file,
#   scale = "row",
#   color = colorRampPalette(brewer.pal(9, "Blues"))(100),
#   show_colnames = TRUE,
#   cluster_rows = FALSE, cluster_cols = FALSE
# )

# Investigate
file1 <- sprintf(
  "data/simulated/scrna/loc0/splat-scale_%.2f-%02d.rds",
  0.01, 1
)
splat <- readRDS(file1)
hist(log(rowData(splat)$BatchFacBatch1))

# ax <- ggplot(bal_summary, aes(x = batch_scale, y = mean, colour = metric)) +
#   facet_wrap(~ metric, ncol = 1, scales = "free_y") +
#   geom_line(linewidth = lw, show.legend = FALSE) +
#   geom_errorbar(
#     aes(ymin = mean - sd, ymax = mean + sd),
#     width = 0.006, linewidth = lw, color =  "grey55"
#   )
# file <- "tmp/fig/splat/splat-bal.jpg"
# ggsave(file, ax, width = 3.5, height = 7)
# 
# ax <- ggplot(imbal_summary, aes(x = batch_scale, y = mean, colour = metric)) +
#   facet_wrap(~ metric, ncol = 1, scales = "free_y") +
#   geom_line(linewidth = lw, show.legend = FALSE) +
#   geom_errorbar(
#     aes(ymin = mean - sd, ymax = mean + sd),
#     width = 0.006, linewidth = lw, color = "grey55" 
#   )
# file <- "tmp/fig/splat/splat-imbal.jpg"
# ggsave(file, ax, width = 3.5, height = 7)

# ax <- ggplot(long_imbal, aes(x = batch_scale, y = value, colour = metric)) +
#   facet_wrap(~ metric, ncol = 1, scales = "free_y") +
#   geom_point() +
#   stat_summary(fun.data = "mean_sdl", geom = "errorbar")
# file <- "tmp/fig/splat/splat-imbal1.jpg"
# ggsave(file, ax)
