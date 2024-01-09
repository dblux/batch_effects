library(ggplot2)
library(cowplot)
theme_set(theme_bw(base_size = 7))

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


cex <- 0.5
metric_ord <- c("RVP", "gPCA", "PVCA", "CMS", "kBET", "LISI")
metric_cols <- ggplot_palette(6)
metric_cols[6] <- "#BF80FF"
names(metric_cols) <- metric_ord

# CPU time
file <- "tmp/benchmark/time.txt" 
time <- read.table(file, header = T, stringsAsFactors = T)
time$cpu.time <- time$user + time$system
time$metric <- factor(time$metric, levels = metric_ord)

xlab <- "Number of samples"
ylab <- "CPU time (s)"
ax <- ggplot(time, aes(x = n, y = cpu.time, col = metric)) +
  geom_point(cex = CEX, show.legend = FALSE) +
  geom_line(linewidth = CEX) +
  labs(x = xlab, y = ylab, col = "Metric") +
  theme(
    # legend.spacing.y = unit(1, "pt"),
    legend.key.size = unit(4, "mm"),
    legend.position = c(0.05, 1),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.key = element_rect(fill="transparent"),
    legend.background = element_rect(fill="transparent"),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = metric_cols)
#   scale_y_continuous(trans = "log10")

file <- "tmp/fig/time.jpg"
ggsave(file, ax, width = 4, height = 2.5)

# Peak memory
file <- "tmp/benchmark/memory.txt"
memory <- read.table(file, header = T, stringsAsFactors = T)
memory$maxsize <- memory$maxsize / 1e6
memory$metric <- factor(memory$metric, levels = metric_ord)

xlab <- "Number of samples"
ylab <- "Peak memory usage (GB)"
ax <- ggplot(memory, aes(x = n, y = maxsize, col = metric)) +
  geom_point(cex = CEX, show.legend = F) +
  geom_line(linewidth = CEX) +
  labs(x = xlab, y = ylab, col = "Metric") +
  theme(
    legend.key.size = unit(4, "mm"),
    legend.position = c(0.05, 1),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.key = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = metric_cols)
#   scale_y_continuous(trans = "log10")

file <- "tmp/fig/memory.jpg"
ggsave(file, ax, width = 4, height = 2.5)
