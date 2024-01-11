library(magrittr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(Biobase)
library(pvca)
theme_set(theme_bw(base_size = 7))
# source files
src_files <- list.files("R", full.names = TRUE)
for (f in src_files) {
  source(f)
  cat(sprintf('Sourced file: %s\n', f))
}

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

# log_nozero <- log_maqc[rowSums(log_maqc == 0) == 0, ]

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

# create expressionset from biobase package
print(fake_metadata)
fake_pheno_data <- new("AnnotatedDataFrame", data = fake_metadata, varMetadata = meta_metadata)
fakebatch_eset <- ExpressionSet(assayData = as.matrix(maqc_fakebatch), phenoData = fake_pheno_data)

pvca_obj <- pvcaBatchAssess(eset, c('batch_info', 'class_info'), 0.6)
var_prop <- as.vector(pvca_obj$dat)
factors <- as.vector(pvca_obj$label)
print(factors)
rounded_prop <- sapply(var_prop, round, digits = 4)
do.call(paste, as.list(c(rounded_prop, sep = '; ')))

# Metadata
METADATA_SID <- "../../relapse_prediction/data/GSE67684/processed/metadata/sid-metadata_v2.tsv"
METADATA_PID <- "../../relapse_prediction/data/GSE67684/processed/metadata/pid-metadata_v7.tsv"
metadata_sid <- read.table(METADATA_SID, sep = "\t")
metadata_pid <- read.table(METADATA_PID, sep = "\t", row.names = 1, quote = '"')

## Data
# Removed outliers, patients with timepoints from different batches and batch 5
SUBSET_RPATH <- "../../relapse_prediction/data/GSE67684/processed/subset_yeoh.tsv"
raw_yeoh <- read.table(SUBSET_RPATH, sep = "\t")

# Metadata
metadata_sid$label <- as.factor(metadata_sid$label)
levels(metadata_sid$label) <- c('Remission', 'Relapse')
metadata_sid$batch_info <- as.factor(metadata_sid$batch_info) 
metadata_pid$label <- as.factor(metadata_pid$label)
levels(metadata_pid$label) <- c('Remission', 'Relapse')

# SCALE->REMOVE->FILTER->LOG
scaled_yeoh <- normaliseMeanScaling(raw_yeoh)
selected_yeoh <- removeProbesets(scaled_yeoh)
yeoh <- log2_transform(filterProbesets(selected_yeoh, 0.7, metadata_sid))

# All features
yeoh_allps <- log2_transform(scaled_yeoh)
yeoh_unfltr <- log2_transform(selected_yeoh)

metadata_yeoh <- metadata_sid[colnames(yeoh), ]
metadata_telaml1 <- subset(
  metadata_sid,
  subtype == 'TEL-AML1' &
  label == 'Remission' &
  batch_info %in% c(2, 9)
)

rvp_sel <- RVP(t(yeoh), metadata_yeoh$batch_info, metadata_yeoh$class_info)
rvp_all <- RVP(t(yeoh_allps), metadata_yeoh$batch_info, metadata_yeoh$class_info)

del_feat <- setdiff(rownames(yeoh_allps), rownames(yeoh))
yeoh_del <- yeoh_allps[del_feat, ]

bal_sids <- get_sid(c(
  'P022', 'P023', 'P024', 'P025',
  'P099', 'P106', 'P120', 'P121'
))
imbal_sids <- c(
  'P022_D0', 'P023_D0', 'P024_D0', 'P025_D0', 'P026_D0', 'P022_D8', 'P023_D8', 'P024_D8',
  'P099_D0', 'P106_D0', 'P120_D0', 'P099_D8', 'P106_D8', 'P120_D8', 'P121_D8', 'P127_D8'
)

yeoh_bal <- yeoh[, bal_sids]
yeoh_imbal <- yeoh[, imbal_sids]

fakebatch_sids <- get_sid(c(
  'P022', 'P023', 'P024', 'P025',
  'P026', 'P027', 'P035', 'P036'
))
yeoh_fakebatch <- yeoh[, fakebatch_sids]

fake_metadata <- metadata_sid[fakebatch_sids, ]
fake_metadata1 <- fake_metadata
fake_metadata2 <- fake_metadata
fake_metadata1[c(1:4, 9:12), 1] <- 1
fake_metadata2[c(1:6, 11:16), 1] <- 1

table(metadata_sid$batch_info)

metadata_yeoh1 <- subset(metadata_yeoh, batch_info %in% c(2, 9))
pidx <- rownames(metadata_yeoh1)

head(metadata_sid)

hist(yeoh[,1], breaks = 30)

ggplot_pca(yeoh[, pidx], metadata_yeoh1, col = 'batch_info', pch = 'class_info')

i <- i + 1
print(i)
feature <- data.frame(
  value = as.vector(data.matrix(yeoh[i, pidx])),
  batch = metadata_yeoh1$batch_info,
  class = metadata_yeoh1$class_info,
  subtype = metadata_yeoh1$subtype
)
ggplot(feature) +
  geom_point(
    aes(x = batch, y = value, col = subtype, pch = class),
    position = position_jitterdodge(jitter.width = .1),
    show.legend = FALSE
  ) +
  ylim(0, 15)



table(metadata_sid$batch_info, metadata_sid$class_info)

ggplot_pca(yeoh_bal, metadata_sid, col = 'batch_info', pch = 'class_info')

rbinom(100, 1, .8)

metrics1 <- quantify_batch(yeoh_bal, metadata_sid)

metrics2 <- quantify_batch(yeoh_imbal, metadata_sid)

metrics3 <- quantify_batch(yeoh_fakebatch, fake_metadata1)

metrics4 <- quantify_batch(yeoh_fakebatch, fake_metadata2)

metrics4

metadata_yeoh <- metadata_sid[colnames(yeoh), ]

# create expressionset from biobase package
meta_metadata <- data.frame(labelDescription = colnames(metadata_sid))
pheno_data <- new("AnnotatedDataFrame", data = metadata_yeoh, varMetadata = meta_metadata)
yeoh_eset <- ExpressionSet(assayData = as.matrix(yeoh), phenoData = pheno_data)

# TODO
eset <- yeoh_eset[, imbal_sids]

pvca_obj <- pvcaBatchAssess(eset, c('batch_info', 'class_info'), 0.6)
var_prop <- as.vector(pvca_obj$dat)
factors <- as.vector(pvca_obj$label)
print(factors)
rounded_prop <- sapply(var_prop, round, digits = 4)
do.call(paste, as.list(c(rounded_prop, sep = '; ')))

# create expressionset from biobase package
print(fake_metadata2)
fake_pheno_data <- new("AnnotatedDataFrame", data = fake_metadata2, varMetadata = meta_metadata)
fakebatch_eset <- ExpressionSet(assayData = as.matrix(yeoh_fakebatch), phenoData = fake_pheno_data)

pvca_obj <- pvcaBatchAssess(fakebatch_eset, c('batch_info', 'class_info'), 0.6)
var_prop <- as.vector(pvca_obj$dat)
factors <- as.vector(pvca_obj$label)
print(factors)
rounded_prop <- sapply(var_prop, round, digits = 4)
do.call(paste, as.list(c(rounded_prop, sep = '; ')))

# WESTLAKE
# creating new subsets with more features (DDA)
# feature selection: all features with missing values were removed
# normalisation? check DDA pipeline
file <- "data/westlake/raw/DDA/A549.csv"
a549 <- read.csv(file, row.names = 1)
file <- "data/westlake/raw/DDA/K562.csv"
k562 <- read.csv(file, row.names = 1)

symbols <- intersect(rownames(a549), rownames(k562))
dda <- cbind(k562[symbols, ], a549[symbols, ])
dda[is.na(dda)] <- 0
dda_class <- sapply(
  strsplit(colnames(dda), "_"),
  function(l) substring(l[[2]], 1, 4)
)
dda_sel <- remove_sparse(dda, 0.3, dda_class)

file <- "data/westlake/processed/metadata/with-balanced.csv"
metadata_with_bal <- read.csv(file, row.names = 1)
file <- "data/westlake/processed/metadata/with-severe1.csv"
metadata_with_imbal <- read.csv(file, row.names = 1)
file <- "data/westlake/processed/metadata/without-balanced1.csv"
metadata_without_bal <- read.csv(file, row.names = 1)
file <- "data/westlake/processed/metadata/without-severe1.csv"
metadata_without_imbal <- read.csv(file, row.names = 1)

metadata_with_bal$machine <- as.factor(metadata_with_bal$machine)
metadata_with_imbal$machine <- as.factor(metadata_with_imbal$machine)
metadata_without_bal$machine <- as.factor(metadata_without_bal$machine)
metadata_without_imbal$machine <- as.factor(metadata_without_imbal$machine)

with_bal <- dda_sel[, rownames(metadata_with_bal)]
file <- "data/westlake/processed/sparse/with-bal.csv"
write.csv(with_bal, file)
with_imbal <- dda_sel[, rownames(metadata_with_imbal)]
file <- "data/westlake/processed/sparse/with-imbal.csv"
write.csv(with_imbal, file)
without_bal <- dda_sel[, rownames(metadata_without_bal)]
file <- "data/westlake/processed/sparse/without-bal.csv"
write.csv(without_bal, file)
without_imbal <- dda_sel[, rownames(metadata_without_imbal)]
file <- "data/westlake/processed/sparse/without-imbal.csv"
write.csv(without_imbal, file)

datasets <- list(
  with_bal = list(X = with_bal, metadata = metadata_with_bal),
  with_imbal = list(X = with_imbal, metadata = metadata_with_imbal),
  without_bal = list(X = without_bal, metadata = metadata_without_bal),
  without_imbal = list(X = without_imbal, metadata = metadata_without_imbal)
)
file <- "data/westlake/processed/westlake-datasets.rds"
saveRDS(datasets, file)

# view results
file <- "tmp/microarray/westlake/results.rds"
results <- readRDS(file)
str(results)
