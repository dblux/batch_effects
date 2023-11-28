library(ggplot2)
library(cowplot)
theme_set(theme_bw(base_size = 7))


SIZE <- 0.2
KEYSIZE <- 9

# CPU time
file <- "tmp/benchmark/time.txt" 
time <- read.table(file, header = T, stringsAsFactors = T)
time$cpu.time <- time$user + time$system

xlab <- "Number of samples"
ylab <- "CPU time (s)"

ax <- ggplot(time) +
  geom_point(
    aes(x = n, y = cpu.time, col = metric),
    cex = SIZE, show.legend = FALSE
  ) +
  geom_line(
    aes(x = n, y = cpu.time, col = metric), linewidth = SIZE
  ) +
  labs(x = xlab, y = ylab, col = "Metric") +
  theme(
    # legend.spacing.y = unit(1, "pt"),
    legend.key.size = unit(KEYSIZE, "pt"),
    legend.position = c(.25, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.key = element_rect(fill="transparent"),
    legend.background = element_rect(fill="transparent"),
    legend.title = element_blank()
  )
#   scale_y_continuous(trans = "log10")

file <- "tmp/fig/time.png"
ggsave(file, ax, width = 2.5, height = 2)

perf_rvp <- subset(performance, metric == "RVP")

# CPU time: RVP
xlab <- "Number of samples"
ylab <- "Runtime (s)"

ax <- ggplot(perf_rvp, aes(x = samples, y = runtime, col = metric)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  labs(x = xlab, y = ylab, col = "Metric") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = c(.25, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.key = element_rect(fill="transparent"),
    legend.background = element_rect(fill="transparent"),
    legend.title = element_blank()
  )
#   scale_y_continuous(trans = "log10")


# Peak memory
file <- "tmp/benchmark/memory.txt"
memory <- read.table(file, header = T, stringsAsFactors = T)
memory$maxsize <- memory$maxsize / 1e6

xlab <- "Number of samples"
ylab <- "Peak memory usage (GB)"

ax <- ggplot(memory) +
  geom_point(aes(x = n, y = maxsize, col = metric), cex = SIZE, show.legend = F) +
  geom_line(aes(x = n, y = maxsize, col = metric), linewidth = SIZE) +
  labs(x = xlab, y = ylab, col = "Metric") +
  theme(
    legend.key.size = unit(KEYSIZE, "pt"),
    legend.position = c(.25, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.key = element_rect(fill="transparent"),
    legend.background = element_rect(fill="transparent"),
    legend.title = element_blank()
  )
#   scale_y_continuous(trans = "log10")

file <- "tmp/fig/memory.png"
ggsave(file, ax, width = 2.5, height = 2)
