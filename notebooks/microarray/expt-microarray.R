library(dplyr)
library(tibble)
library(tidyr)

library(Biobase)
library(pvca)

library(ggplot2)
library(cowplot)
theme_set(theme_bw())


quantify_batch <- function(X, metadata) {
  batch <- metadata[colnames(X), 'batch']
  class <- metadata[colnames(X), 'class']

  rvp <- RVP(t(X), batch, class)
  gpca <- gPCA.batchdetect(t(X), batch)
    
  # create expressionset from biobase package
  meta_metadata <- data.frame(labelDescription = colnames(metadata))
  pheno_data <- new(
    "AnnotatedDataFrame",
    data = metadata, varMetadata = meta_metadata
  )
  eset <- ExpressionSet(assayData = as.matrix(X), phenoData = pheno_data)
  pvca_obj <- pvcaBatchAssess(eset, c('batch', 'class'), 0.6)
  var_pcts <- as.vector(pvca_obj$dat)
  names(var_pcts) <- as.vector(pvca_obj$label)

  c(
    gpca = gpca$delta,
    pvca = var_pcts['batch'],
    rvp = rvp
  )
}



# MAQC
file_maqc <- "../data/MAQC-I/processed/mas5_original-ref.tsv"
raw_maqc <- read.table(file_maqc, sep = "\t", header = T, row.names = 1)

# MAQC metadata
# Class A - Universal Human Reference RNA (UHRR)
# Class B - Human Brain Reference RNA (HBRR)
batch_info <- as.factor(rep(1:6, each = 10))
class_info <- rep(rep(LETTERS[1:2], each = 5), 6)
metadata_maqc <- data.frame(
  batch_info, class_info,
  row.names = colnames(raw_maqc)
)

# SCALE->REMOVE->FILTER->LOG
scaled_maqc <- raw_maqc %>%
  normaliseMeanScaling()

log_maqc <- scaled_maqc %>%
  removeProbesets() %>%
  filterProbesets(0.7, metadata_maqc) %>%
  log2_transform()

ggplot_pca(
  maqc_bal, metadata_maqc,
  col = 'batch_info', pch = 'class_info'
)

feat_vars <- apply(log_maqc, 1, var)
i <- 1
idx <- order(-feat_vars)
print(feat_vars[idx[i]])
i <- i + 1
print(i)
feature <- data.frame(
  value = as.vector(data.matrix(log_maqc[idx[i], ])),
  batch = metadata_maqc$batch_info,
  class = metadata_maqc$class_info
)
ggplot(feature) +
  geom_point(
    aes(x = batch, y = value, pch = class, col = batch),
    position = position_jitterdodge(jitter.width = 2),
    show.legend = FALSE
  )

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
pvca_obj <- pvcaBatchAssess(eset, c('batch_info', 'class_info'), 0.6)
var_prop <- as.vector(pvca_obj$dat)
factors <- as.vector(pvca_obj$label)
print(factors)
rounded_prop <- sapply(var_prop, round, digits = 4)
do.call(paste, as.list(c(rounded_prop, sep = '; ')))

# Ma-Spore ALL

# TODO: Include rest of the ipynb notebook
