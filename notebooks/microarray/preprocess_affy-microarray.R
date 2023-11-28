library(affy)

path <- '../data/GSE24061/raw/GPL1261/'
files <- list.files(path, full.names = T)
affybatch <- ReadAffy(filenames = files)

# Scan dates of microarrays
scan_dates <- affybatch@protocolData@data
accessions <- sapply(
  strsplit(rownames(scan_dates), '_'),
  function(x) x[1]
)
rownames(scan_dates) <- accessions
file <- '../data/GSE24061/processed/scan_dates.tsv'
write.table(scan_dates, file, quote = F, sep = '\t')

# Processing CEL
# MAS5 without trimmed mean scaling
data_eset <- mas5(affybatch, normalize = F)
# MAS5 detection call object
calls_eset <- mas5calls(affybatch)

# # RMA without quantile normalisation
# # Use expresso function to mix and match algorithm
# rma_data <- expresso(raw_data,
#                     bgcorrect.method = "rma",
#                     normalize = F,
#                     pmcorrect.method = "pmonly",
#                     summary.method = "medianpolish")

# Extract expression data
raw_exprs <- exprs(data_eset)
colnames(raw_exprs) <- accessions

# Assigns detection calls based on default threshold
# M call: 0.04 < p-value <= 0.06
mas5_calls <- exprs(calls_eset)
colnames(mas5_calls) <- accessions

# If call == "P" cell -> 1
indicator_mat <- (mas5_calls == "P") * 1
# mas5 data preserving only "P" calls
filtered_exprs <- raw_exprs * indicator_mat

# Extracts pvalues from MAS5 call object
mas5_pvalues <- assayData(calls_eset)[["se.exprs"]]
colnames(mas5_pvalues) <- accessions

file1 <- '../data/GSE24061/processed/GSE24061-mas5_exprs.tsv'
write.table(filtered_exprs, file1, sep = "\t", quote = F)
file2 <- '../data/GSE24061/processed/GSE24061-mas5_calls.tsv'
write.table(mas5_calls, file2, sep = "\t", quote = F)
file3 <- '../data/GSE24061/processed/GSE24061-mas5_pvalues.tsv'
write.table(mas5_pvalues, file3, sep = "\t", quote = F)

path <- '../data/GSE6116/raw/GPL1261/'
files <- list.files(path, full.names = T)
affybatch <- ReadAffy(filenames = files)

# Scan dates of microarrays
scan_dates <- affybatch@protocolData@data
accessions <- substring(rownames(scan_dates), 0, 9)
rownames(scan_dates) <- accessions

file <- '../data/GSE6116/processed/scan_dates.tsv'
write.table(scan_dates, file, quote = F, sep = '\t')

# Processing CEL
# MAS5 without trimmed mean scaling
data_eset <- mas5(affybatch, normalize = F)
# MAS5 detection call object
calls_eset <- mas5calls(affybatch)

# # RMA without quantile normalisation
# # Use expresso function to mix and match algorithm
# rma_data <- expresso(raw_data,
#                     bgcorrect.method = "rma",
#                     normalize = F,
#                     pmcorrect.method = "pmonly",
#                     summary.method = "medianpolish")

# Extract expression data
raw_exprs <- exprs(data_eset)
colnames(raw_exprs) <- accessions

# Assigns detection calls based on default threshold
# M call: 0.04 < p-value <= 0.06
mas5_calls <- exprs(calls_eset)
colnames(mas5_calls) <- accessions

# If call == "P" cell -> 1
indicator_mat <- (mas5_calls == "P") * 1
# mas5 data preserving only "P" calls
filtered_exprs <- raw_exprs * indicator_mat

# Extracts pvalues from MAS5 call object
mas5_pvalues <- assayData(calls_eset)[["se.exprs"]]
colnames(mas5_pvalues) <- accessions

file1 <- '../data/GSE6116/processed/GSE6116-mas5_exprs.tsv'
write.table(filtered_exprs, file1, sep = "\t", quote = F)
file2 <- '../data/GSE6116/processed/GSE6116-mas5_calls.tsv'
write.table(mas5_calls, file2, sep = "\t", quote = F)
file3 <- '../data/GSE6116/processed/GSE6116-mas5_pvalues.tsv'
write.table(mas5_pvalues, file3, sep = "\t", quote = F)

# # Visualise each microarray
# # Remove chips with too many artefacts
# DIR_WPATH <- sub("raw", "img", CEL_DIRPATH)
# print(DIR_WPATH)

# # Creates new directory
# dir.create(DIR_WPATH)

# # Saves each microarray image in folder
# for (i in 1:length(fpaths)) {
#   wpath <- sprintf("%s%03d.jpg", DIR_WPATH, i)
#   jpeg(wpath)
#   image(raw_data[, i])
#   dev.off()
# }

# # Plot boxplots of probe intensities for each microarray
# par(mar = c(10,4,4,4),
#     cex.axis = 0.8)
# boxplot(raw_data, las = 2)

# # Plot density curves of microarrays
# hist(raw_data)
