library(profmem)
library(magrittr)
library(dplyr)
library(scater)
library(pryr)

source("R/HVP.R")
source("R/utils.R")


data <- readRDS("data/simulated/profiling-data.rds")
data <- logNormCounts(data)
dim(data)
object_size(logcounts(data))

profile <- profmem(
  .HVP(
    logcounts(data),
    batch = data$Batch,
    cls = data$Group
  )
)

filter_profile <- profile %>%
  filter(bytes > 2e6)
  # arrange(desc(bytes))
print(filter_profile)

for (i in 1:nrow(filter_profile)) {
  filter_profile[[i, 3]] %>%
    rev() %>%
    paste(collapse = '() -> ') %>%
    paste0('()') %>%
    print()
  print(sprintf("%.1fMB", filter_profile[i, 2] / 1e6))
}

### Investigate ###

x <- logcounts(data)

profmem({
  split_cols(x, data$Batch)
  # ss <- rowSums((x - rowMeans(x)) ^ 2)
}) %>%
  filter(bytes > 2e6) %>%
  print()

obj <- HVP(data, "Batch", "Group")

library(Matrix)

obj <- .HVP_sparseMatrix(
  logcounts(data),
  batch = data$Batch,
  cls = data$Group
)
