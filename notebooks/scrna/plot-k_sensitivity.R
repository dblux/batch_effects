library(magrittr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_bw(base_size = 7))
src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


METRIC_ORD <- c("HVP", "gPCA", "PVCA", "CMS", "kBET", "LISI")
GGCOLS <- ggplot_palette(6)
GGCOLS[6] <- "#BF80FF"
names(GGCOLS) <- METRIC_ORD


# k-sensitivity
beffs_labs <- c("With batch effects", "Without batch effects")
names(beffs_labs) <- c("With", "Without")
imbalance_labs <- c("Batch-class balanced", "Batch-class imbalanced")
names(imbalance_labs) <- c("Balanced", "Imbalanced")
custom_theme <- theme(
  legend.key.size = unit(4, "mm"),
  title = element_text(size = 6),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_text(size = 6)
)

# villani
file <- "tmp/scrna/villani/villani-k.csv"
villani <- read.csv(file)
villani_k <- villani %>% 
  gather(key = "Condition", "Value", -Metric, -Param) %>%
  separate_wider_delim(
    cols = "Condition", delim = ".",
    names = c("Batch.Effects", "Imbalance")
  )
y_offset <- -1
villani_k[villani_k$Metric == "LISI", "Value"] <-
  villani_k[villani_k$Metric == "LISI", "Value"] + y_offset
left_ylab <- "CMS, kBET"
# ordered_lvls <- colnames(villani)[-c(1, 2)]
# villani_k$condition <- factor(villani_k$condition, levels = ordered_lvls)
ax <- ggplot(villani_k, aes(x = Param, y = Value, col = Metric)) +
  facet_grid(
    Imbalance ~ Batch.Effects, scales = "fixed",
    labeller = labeller(Imbalance = imbalance_labs, Batch.Effects = beffs_labs)
  ) +
  geom_line(show.legend = FALSE) +
  geom_point(cex = 0.5, show.legend = FALSE) +
  scale_y_continuous(
    name = left_ylab, limits = c(0, 1),
    sec.axis = sec_axis(trans = ~ . - y_offset, name = "LISI")
  ) +
  scale_color_manual(values = GGCOLS) +
  labs(title = "Villani et al.", x = "k") +
  custom_theme
file <- "tmp/fig/real/villani-k.pdf"
ggsave(file, ax, width = 2.8, height = 2.6)

# halfmix
file <- "tmp/scrna/halfmix/halfmix-k.csv"
halfmix <- read.csv(file)
halfmix_k <- halfmix %>% 
  gather(key = "Condition", "Value", -Metric, -Param) %>%
  separate_wider_delim(
    cols = "Condition", delim = ".",
    names = c("Batch.Effects", "Imbalance")
  )
y_offset <- -1.2
halfmix_k[halfmix_k$Metric == "LISI", "Value"] <-
  halfmix_k[halfmix_k$Metric == "LISI", "Value"] + y_offset
left_ylab <- "CMS, kBET"
# ordered_lvls <- colnames(halfmix)[-c(1, 2)]
# halfmix_k$condition <- factor(halfmix_k$condition, levels = ordered_lvls)
ax <- ggplot(halfmix_k, aes(x = Param, y = Value, col = Metric)) +
  facet_grid(
    Imbalance ~ Batch.Effects, scales = "fixed",
    labeller = labeller(Imbalance = imbalance_labs, Batch.Effects = beffs_labs)
  ) +
  geom_line(show.legend = FALSE) +
  geom_point(cex = 0.5, show.legend = FALSE) +
  scale_y_continuous(
    name = left_ylab, limits = c(0, 1),
    sec.axis = sec_axis(trans = ~ . - y_offset, name = "LISI")
  ) +
  scale_color_manual(values = GGCOLS) +
  labs(title = "Zheng et al.", x = "k") +
  custom_theme
file <- "tmp/fig/real/halfmix-k.pdf"
ggsave(file, ax, width = 2.8, height = 2.6)

# cellbench
file <- "tmp/scrna/cellbench/cellbench-k.csv"
cellbench <- read.csv(file)
cellbench_k <- cellbench %>% 
  gather(key = "Condition", "Value", -Metric, -Param) %>%
  separate_wider_delim(
    cols = "Condition", delim = ".",
    names = c("Batch.Effects", "Imbalance")
  )
y_offset <- -1
y_scale <- 0.5
cellbench_k[cellbench_k$Metric == "LISI", "Value"] <-
  (cellbench_k[cellbench_k$Metric == "LISI", "Value"] + y_offset) * y_scale
left_ylab <- "CMS, kBET"
# ordered_lvls <- colnames(cellbench)[-c(1, 2)]
# cellbench_k$condition <- factor(cellbench_k$condition, levels = ordered_lvls)
ax <- ggplot(cellbench_k, aes(x = Param, y = Value, col = Metric)) +
  facet_grid(
    Imbalance ~ Batch.Effects, scales = "fixed",
    labeller = labeller(Imbalance = imbalance_labs, Batch.Effects = beffs_labs)
  ) +
  geom_line(show.legend = TRUE) +
  geom_point(cex = 0.5, show.legend = FALSE) +
  scale_y_continuous(
    name = left_ylab, limits = c(0, 1),
    sec.axis = sec_axis(trans = ~ . / y_scale - y_offset, name = "LISI")
  ) +
  scale_color_manual(values = GGCOLS) +
  labs(title = "Tian et al.", x = "k") +
  custom_theme
file <- "tmp/fig/real/cellbench-k.pdf"
ggsave(file, ax, width = 3.3, height = 2.6)
