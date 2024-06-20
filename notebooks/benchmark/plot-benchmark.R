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
metric_ord <- c("HVP", "gPCA", "PVCA", "CMS", "kBET", "LISI")
metric_cols <- ggplot_palette(6)
metric_cols[6] <- "#BF80FF"
names(metric_cols) <- metric_ord

# CPU time
file <- "tmp/benchmark/hvp/time3.txt"
time <- read.table(file, header = T, stringsAsFactors = T)
# time <- subset(time, metric != "HVPS")
time$cpu.time <- time$user + time$system

time[order(time$n, time$metric), ]

time$metric <- factor(time$metric, levels = metric_ord)

xlab <- "Number of samples"
ylab <- "CPU time (s)"
ax <- ggplot(time, aes(x = n, y = cpu.time, col = metric)) +
  geom_point(cex = cex, show.legend = FALSE) +
  geom_line(linewidth = cex) +
  labs(x = xlab, y = ylab, col = "Metric") +
  theme(
    # legend.spacing.y = unit(1, "pt"),
    legend.key.size = unit(5, "mm"),
    # legend.position = c(0.05, 1),
    # legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.key = element_rect(fill="transparent"),
    legend.background = element_rect(fill="transparent"),
  ) +
  scale_color_manual(values = metric_cols) +
  scale_y_continuous(trans = "log10")
file <- "tmp/fig/log-cputime.pdf"
ggsave(file, ax, width = 3.5, height = 2)

# Peak memory
# file <- "tmp/benchmark/rvp_sparse/rerun/memory.txt"
file <- "tmp/benchmark/hvp/memory3.txt"
memory <- read.table(file, header = T, stringsAsFactors = T)
memory <- subset(memory, metric != "HVPS")
memory$maxsize <- memory$maxsize / 1e6

memory$metric <- factor(memory$metric, levels = metric_ord)

memory[order(memory$n, memory$metric), ]

xlab <- "Number of samples"
ylab <- "Peak memory usage (GB)"
ax <- ggplot(memory, aes(x = n, y = maxsize, col = metric)) +
  geom_point(cex = cex, show.legend = F) +
  geom_line(linewidth = cex) +
  labs(x = xlab, y = ylab, col = "Metric") +
  theme(
    legend.key.size = unit(5, "mm"),
    # legend.position = c(0.05, 1),
    # legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.key = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),
  ) +
  scale_color_manual(values = metric_cols)
file <- "tmp/fig/memory.pdf"
ggsave(file, ax, width = 3.5, height = 2)
