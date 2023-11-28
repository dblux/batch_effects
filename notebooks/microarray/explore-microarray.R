library(magrittr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
theme_set(theme_bw())

src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


# MAQC
file_maqc <- "data/MAQC-I/processed/mas5_original-ref.tsv"
raw_maqc <- read.table(file_maqc, sep = "\t", header = T, row.names = 1)

# MAQC metadata
# Class A - Universal Human Reference RNA (UHRR)
# Class B - Human Brain Reference RNA (HBRR)
batch <- as.factor(rep(1:6, each = 10))
class <- rep(rep(LETTERS[1:2], each = 5), 6)
metadata_maqc <- data.frame(
  batch, class,
  row.names = colnames(raw_maqc)
)

# normalise -> remove -> filter -> log
scaled_maqc <- raw_maqc %>%
  normaliseMeanScaling() %>%
  remove_affymetrix()
full_maqc <- log2_transform(scaled_maqc)

log_maqc <- scaled_maqc %>%
  remove_sparse(0.3, metadata_maqc$class) %>%
  log2_transform()

# MAR only present in data sets with severe batch effects (e.g. diff tech)
# Block drop outs due to tech that detects unique features that others cannot
# MNAR seems to be the main cause of missing values
# Higher detection limit in some batches are an example of MNAR
# Non-BEAMs may have batch effects due to remaining values with batch effects
# Missing values in scRNA-Seq are mostly MNAR and block drop outs (MAR).
# Both MNAR and MAR are BEAMs. Most missing values are BEAMs
a_maqc <- subset_cols(full_maqc, metadata_maqc, class == "A")

pct_zero <- rowSums(a_maqc == 0) / ncol(a_maqc)
nonzero_mean <- apply(a_maqc, 1, function(x) mean(x[x > 0]))
dropout <- na.omit(data.frame(nonzero_mean, pct_zero))
dropout_missing <- subset(dropout, subset = pct_zero > 0)

ax <- ggplot(dropout_missing, aes(x = nonzero_mean, y = pct_zero)) +
  geom_hex()
file <- "tmp/fig/maqc_a-dropout_nonzero_mean1.png"
ggsave(file, ax)

# inspection
funny_features <- rownames(subset(
  dropout_missing,
  subset = nonzero_mean > 8 & pct_zero > 0.8
))
length(funny_features)

# marginal
llim <- 9
rlim <- llim + 0.5
file <- sprintf("tmp/fig/maqc_a-marginal_%.1f.png", llim)
pct_zero_marginal <- pct_zero[which(nonzero_mean > llim & nonzero_mean < rlim)]
png(file)
hist(
  pct_zero_marginal,
  xlab = "Percentage of zeros",
  main = sprintf(
    "Features with mean in (%.1f, %.1f); n = %d",
    llim, rlim, length(pct_zero_marginal)
  )
)
dev.off()

png("tmp/fig/maqc_a-dropout.png")
hist(pct_zero)
dev.off()

a_maqc_fltr <- a_maqc[-which(rowSums(a_maqc == 0) %in% c(0, 30)), ]

pheatmap(
  a_maqc_fltr, color = colorRampPalette(c("white", "navy"))(50),
  cluster_rows = F, cluster_cols = F
)

# total percentage zero
sum(a_maqc_fltr[, 1:5] == 0) / length(as.matrix(a_maqc_fltr[, 1:5]))
sum(a_maqc_fltr[, 6:10] == 0) / length(as.matrix(a_maqc_fltr[, 6:10]))
sum(a_maqc_fltr[, 11:15] == 0) / length(as.matrix(a_maqc_fltr[, 11:15]))
sum(a_maqc_fltr[, 16:20] == 0) / length(as.matrix(a_maqc_fltr[, 16:20]))
sum(a_maqc_fltr[, 21:25] == 0) / length(as.matrix(a_maqc_fltr[, 21:25]))

# BEAMs
x1 <- rowSums(a_maqc_fltr[, 1:5] == 0) / 5 
x2 <- rowSums(a_maqc_fltr[, 6:10] == 0) / 5 
x3 <- rowSums(a_maqc_fltr[, 11:15] == 0) / 5 
x4 <- rowSums(a_maqc_fltr[, 16:20] == 0) / 5
x5 <- rowSums(a_maqc_fltr[, 21:25] == 0) / 5 
batch_pct_zero <- data.frame(x1, x2, x3, x4, x5)
batch_maxmin <- apply(batch_pct_zero, 1, function(x) max(x) - min(x))
features <- which(batch_maxmin >= 0.8)


# 2000 BEAMs, 12000
pheatmap(
  a_maqc_fltr[-features, ], color = colorRampPalette(c("white", "navy"))(50),
  cluster_rows = F, cluster_cols = F
)

ggplot_pca(
  a_maqc_fltr[-features, ], metadata_maqc,
  col = 'batch', pch = 'class'
)

nonbeam <- a_maqc_fltr[-features, ]
beam <- a_maqc_fltr[features, ]
# single feature
i <- 0
i <- i + 1
print(i)
feature <- data.frame(
  value = as.vector(data.matrix(a_maqc_fltr[funny_features[i], ])),
  batch = metadata_maqc$batch,
  class = metadata_maqc$class
)
ggplot(feature) +
  geom_point(
    aes(x = batch, y = value, pch = class, col = batch),
    position = position_jitterdodge(jitter.width = 2),
    show.legend = FALSE
  )

# boxplot
long_a <- gather(a_maqc_fltr, key = "sid")
long_a$batch <- as.factor(substring(long_a$sid, 5, 5))
long_a_fltr <- subset(long_a, subset = value > 0)

ggplot(long_a_fltr, aes(x = sid, y = value, col = batch)) +
  geom_boxplot()

# # Measuring batch effects
# - Different batch-class imbalance
# - No batch effects
# - Different magnitude of batch effects when there is batch-class imbalance
# - Measure when there is different number of features
# - Different batch sizes

# ### Comparisons
# - gPCA delta: 0-1 (proportion of variance)
#     - Weak in quantifying small amounts of batch effects?
#     - Problems when there is no batch effects but there is class imbalance

# subsetting maqc
maqc_bal <- log_maqc[, c(1:4, 6:9, 11:14, 16:19)]
maqc_imbal <- log_maqc[, c(1:3, 6:10, 11:15, 16:18)]
maqc_fakebatch <- log_maqc[, c(1:4, 6:9)]

fake_metadata <- metadata_maqc[colnames(maqc_fakebatch), ]
fake_metadata[3, 1] <- 2
fake_metadata[6, 1] <- 1

# create expressionset from biobase package
meta_metadata <- data.frame(labelDescription = colnames(metadata))
pheno_data <- new("AnnotatedDataFrame", data = metadata, varMetadata = meta_metadata)
maqc_eset <- ExpressionSet(assayData = as.matrix(log_maqc), phenoData = pheno_data)
eset <- maqc_eset[, c(1:4, 6:9, 11:14, 16:19)]
pvca_obj <- pvcaBatchAssess(eset, c('batch', 'class'), 0.6)
var_prop <- as.vector(pvca_obj$dat)
factors <- as.vector(pvca_obj$label)
print(factors)
rounded_prop <- sapply(var_prop, round, digits = 4)
do.call(paste, as.list(c(rounded_prop, sep = '; ')))

# Ma-Spore ALL

# TODO: Include rest of the ipynb notebook
