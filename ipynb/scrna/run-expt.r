library(optparse)


option_list <- list(
  make_option(
    c("-f", "--file"), action = "store",
    help = "File name of datasets.rds file"
  ),
  make_option(
    c("-o", "--outfile"), action = "store",
    help = "File name of output file"
  ),
  make_option(
    c("-k", "--knn"), action = "store", default = NULL,
    help = paste0(
      "k-nn param for kBET and perplexity param for LISI.",
      "k = 0 resorts to default parameters."
    )
  )
)
opt <- parse_args(OptionParser(option_list = option_list))


library(Seurat)
library(kBET)
library(lisi)
src_files <- list.files("R", full.names = TRUE)
cat("Sourcing files:", fill = TRUE)
for (f in src_files) {
  source(f)
  cat(f, fill = TRUE)
}


cat(sprintf("Reading file: %s", opt$file), fill = TRUE)
datasets <- readRDS(opt$file)

if (is.null(opt$k)) {
  cat("Performing RVP...", fill = TRUE)
  results <- lapply(
    datasets, rvp,
    "batch", "celltype", nperm = 100
  )
} else if (as.numeric(opt$k) == 0) {
  k <- as.numeric(opt$k)
  cat(sprintf("Performing kBET and LISI (k = %d)...", k), fill = TRUE)
  results <- lapply(
    datasets, eval_batch,
    "batch", "celltype",
    ret.scores = FALSE, do.rvp = FALSE
  )
} else {
  k <- as.numeric(opt$k)
  cat(sprintf("Performing kBET and LISI (k = %d)...", k), fill = TRUE)
  results <- lapply(
    datasets, eval_batch,
    "batch", "celltype",
    k0 = k, perplexity = k,
    ret.scores = FALSE, do.rvp = FALSE
  )
}

cat(sprintf("Writing file: %s", opt$outfile), fill = TRUE)
saveRDS(results, opt$outfile)
