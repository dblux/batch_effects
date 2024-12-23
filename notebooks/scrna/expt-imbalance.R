src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


file <- "data/simulated/imbalance/datasets-be.rds"
data_be <- readRDS(file)

names(data_be)
data <- data_be[[6]]
keep_features <- remove_sparse(
  assay(data), 0.8, class = data$Group,
  ret.features = TRUE, keep = TRUE 
)
length(keep_features)
data <- data[keep_features, ]
res <- HVP(data, "Batch", "Group")
print(res$HVP)

hvp <- c(0.00756, 0.00731, 0.00581, 0.00503, 0.00352, 0.00305)
file <- "tmp/fig/imbalance-hvp.pdf"
pdf(file)
plot(as.numeric(names(data_be)), hvp, ylim = c(0, 0.008))
dev.off()

dim(data)
pct_zero <- rowSums(assay(data) == 0) / ncol(data)

file <- "tmp/fig/imbalance-pct_zeros.pdf"
pdf(file)
hist(pct_zero)
dev.off()

