library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw(base_size = 7))


ratios <- read.csv("tmp/benchmark/tables.csv")
head(ratios)

# Horizontal bar chart

datanames <- c("villani", "halfmix", "cellbench", "maqc", "yeoh", "westlake")
datanames <- "cellbench"
for (dataname in datanames) {
  result_sub <- subset(result, subset = dataset == dataname)
  ax_snr <- ggplot(result_sub) +
    geom_bar(
      stat = "identity",
      aes(x = log(snr), y = method, fill = metric),
      show.legend = FALSE
    ) +
    labs(y = dataname)
  ax_psr <- ggplot(result_sub) +
    geom_bar(
      stat = "identity",
      aes(x = psr, y = method, fill = metric),
      show.legend = FALSE
    ) +
    labs(y = dataname)
  ax <- plot_grid(ax_snr, ax_psr)

  file <- sprintf("tmp/fig/%s-snr_psr.png", dataname)
  ggsave(file, ax, width = 8, height = 2)
}


left_ylab <- "CMS, kBET"
villani <- read.csv("tmp/villani-k.csv")
villani_k <- gather(villani, key = "condition", "value", -Metric, -Param)
ordered_lvls <- colnames(villani)[-c(1, 2)]
villani_k$condition <- factor(villani_k$condition, levels = ordered_lvls)

villani_ratio <- subset(
  villani_k, subset = condition %in% c("SNR", "PSR")
)
villani_cond <- subset(
  villani_k, subset = !(condition %in% c("SNR", "PSR"))
)
villani_cond[villani$Metric == "LISI", "value"] <-
  villani_cond[villani$Metric == "LISI", "value"] - 1

ax_cond <- ggplot(villani_cond, aes(x = Param, y = value, col = Metric)) +
  facet_wrap(~condition, nrow = 1, scales = "fixed") +
  geom_point() +
  geom_line(aes(y = value)) +
  scale_y_continuous(
    name = left_ylab, limits = c(0, 1),
    sec.axis = sec_axis(trans = ~ . + 1, name = "LISI")
  )
file <- "tmp/fig/villani-cond.png"
ggsave(file, ax_cond, width = 10, height = 2.5)

ax_ratio <- ggplot(villani_ratio, aes(x = Param, y = value, col = Metric)) +
  facet_wrap(~condition, nrow = 1, scales = "free") +
  geom_point() +
  geom_line(aes(y = value)) +
  scale_y_continuous(
    name = left_ylab,
    sec.axis = sec_axis(trans = ~ . + 1, name = "LISI")
  )
file <- "tmp/fig/villani-ratio.png"
ggsave(file, ax_ratio, width = 6, height = 2.5)


halfmix <- read.csv("tmp/halfmix-k.csv")
halfmix_k <- gather(halfmix, key = "condition", "value", -Metric, -Param)
ordered_lvls <- colnames(halfmix)[-c(1, 2)]
halfmix_k$condition <- factor(halfmix_k$condition, levels = ordered_lvls)

halfmix_ratio <- subset(
  halfmix_k, subset = condition %in% c("SNR", "PSR")
)
halfmix_cond <- subset(
  halfmix_k, subset = !(condition %in% c("SNR", "PSR"))
)
y_offset <- -1.2
halfmix_cond[halfmix$Metric == "LISI", "value"] <-
  halfmix_cond[halfmix$Metric == "LISI", "value"] + y_offset

ax_cond <- ggplot(halfmix_cond, aes(x = Param, y = value, col = Metric)) +
  facet_wrap(~condition, nrow = 1, scales = "fixed") +
  geom_point() +
  geom_line(aes(y = value)) +
  scale_y_continuous(
    name = left_ylab, limits = c(0, 1),
    sec.axis = sec_axis(trans = ~ . - y_offset, name = "LISI")
  )
file <- "tmp/fig/halfmix-cond.png"
ggsave(file, ax_cond, width = 12, height = 2.5)

ax_ratio <- ggplot(halfmix_ratio, aes(x = Param, y = value, col = Metric)) +
  facet_wrap(~condition, nrow = 1, scales = "free") +
  geom_point() +
  geom_line(aes(y = value)) +
  scale_y_continuous(
    name = left_ylab,
    sec.axis = sec_axis(trans = ~ . + 1, name = "LISI")
  )
file <- "tmp/fig/halfmix-ratio.png"
ggsave(file, ax_ratio, width = 6, height = 2.5)


cellbench <- read.csv("tmp/cellbench-k.csv")
cellbench <- cellbench[-1, ] # remove RVP
cellbench_k <- gather(cellbench, key = "condition", "value", -Metric, -Param)
ordered_lvls <- colnames(cellbench)[-c(1, 2)]
cellbench_k$condition <- factor(cellbench_k$condition, levels = ordered_lvls)

cellbench_ratio <- subset(
  cellbench_k, subset = condition %in% c("SNR", "PSR")
)
cellbench_cond <- subset(
  cellbench_k, subset = !(condition %in% c("SNR", "PSR"))
)
y_offset <- -1
y_scale <- 0.5
cellbench_cond[cellbench$Metric == "LISI", "value"] <-
  (cellbench_cond[cellbench$Metric == "LISI", "value"] + y_offset) * y_scale

ax_cond <- ggplot(cellbench_cond, aes(x = Param, y = value, col = Metric)) +
  facet_wrap(~condition, nrow = 1, scales = "fixed") +
  geom_point() +
  geom_line(aes(y = value)) +
  scale_y_continuous(
    name = left_ylab, limits = c(0, 2),
    sec.axis = sec_axis(trans = ~ (. / y_scale - y_offset), name = "LISI")
  )
file <- "tmp/fig/cellbench-cond.png"
ggsave(file, ax_cond, width = 10, height = 2.5)

ax_ratio <- ggplot(cellbench_ratio, aes(x = Param, y = value, col = Metric)) +
  facet_wrap(~condition, nrow = 1, scales = "free") +
  geom_point() +
  geom_line(aes(y = value)) +
  scale_y_continuous(
    name = left_ylab,
    sec.axis = sec_axis(trans = ~ . + 1, name = "LISI")
  )
file <- "tmp/fig/cellbench-ratio.png"
ggsave(file, ax_ratio, width = 6, height = 2.5)


file <- "tmp/scrna/villani-results_k20.rds"
results <- readRDS(file)
str(results)
results$negctrl_bal$cms$cms
results$negctrl_imbal$cms$cms
results$bal$cms$cms
results$bal$cms$cms

