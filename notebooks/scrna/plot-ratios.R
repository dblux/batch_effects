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


# Horizontal bar chart
file <- "tmp/snr_psr.csv"
ratios <- read.csv(file)
ratios["log2(SNR)"] <- log2(ratios$SNR)
ratios$Metric <- factor(ratios$Metric, levels = METRIC_ORD)

datasets <- list(
  c("villani", 40, "Villani et al. (scRNA-seq)"),
  c("halfmix", 200, "Zheng et al. (scRNA-seq)"),
  c("cellbench", 40, "Tian et al. (scRNA-seq)"),
  c("maqc", NA, "MAQC (microarray)"),
  c("yeoh", NA, "Yeoh et al. (microarray)"),
  c("westlake", NA, "Wang et al. (proteomics)")
)

for (params in datasets) {
  dataset <- params[1]
  k <- params[2]
  author <- params[3]
  ratios_dataset <- ratios %>%
    subset(
      subset = Dataset == dataset,
      select = c("Metric", "log2(SNR)", "PSR")
    ) %>% 
    gather(key = "Ratio", value = "Value", -Metric) %>%
    na.omit()
  ax <-  ggplot(ratios_dataset) +
      facet_wrap(~Ratio, nrow = 1, scales = "free", drop = TRUE) +
      geom_bar(
        stat = "identity",
        aes(x = Value, y = Metric, fill = Metric),
        show.legend = FALSE
      ) +
      theme(
        title = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      ) +
      scale_fill_manual(values = GGCOLS) +
      labs(title = author) +
      scale_y_discrete(limits = rev)
  file <- sprintf("tmp/fig/real/%s-ratios.pdf", dataset)
  print(file)
  ggsave(file, ax, width = 3.6, height = 0.9)
}

