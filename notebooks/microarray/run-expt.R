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
    c("-m", "--metrics"), action = "store", default = "all",
    help = "Metric to evaluate. Defaults to all."
  ),
  make_option(
    c("-k", "--knn"), action = "store", default = 0,
    help = paste(
      "k-nn param for kBET and perplexity param for LISI.",
      "k = 0 resorts to default parameters."
    )
  ),
  make_option(
    c("-p", "--nperm"), action = "store", default = 0,
    help = "Number of permutations for RVP."
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

nperm <- as.numeric(opt$nperm)
k <- as.numeric(opt$k)
if (k == 0) {
  k.cms <- 50
  k0 <- NULL
  perplexity <- 30
} else {
  k.cms <- k0 <- perplexity <- k
}
# workout which metric to run
if (opt$metrics %in% c("cms", "kbet", "lisi")) {
  metrics <- opt$metrics
} else if (opt$metrics == "all") {
  metrics <- c("rvp", "cms", "kbet", "lisi", "gpca", "pvca")
} else if (opt$metrics == "ckl") {
  metrics <- c("cms", "kbet", "lisi")
} else if (opt$metrics == "rgp") {
  metrics <- c("rvp", "gpca", "pvca")
} else {
  stop("flag:metric is not in acceptable list of values.")
}
cat(sprintf(
  "Params: nperm = %d, k = %d, metrics = %s",
  nperm, k, capture.output(cat(metrics))
), fill = TRUE)



# Imports
files <- list.files("R", full.names = TRUE)
for (file in files) {
  source(file)
  cat(sprintf("Sourced %s", file), fill = TRUE)
}


# Evaluate batch effects
cat(sprintf("Reading file: %s", opt$file), fill = TRUE)
datasets <- readRDS(opt$file)

results <- lapply(
  datasets, eval_batch, "machine", "class",
  metrics = metrics,
  nperm = nperm, k.cms = k.cms, k0 = k0, perplexity = perplexity
)

cat(sprintf("Writing file: %s", opt$outfile), fill = TRUE)
saveRDS(results, opt$outfile)
