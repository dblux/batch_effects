library(testthat)

library(SingleCellExperiment)
library(Matrix)

source("../R/rvp.R")

# SETUP
X <- readRDS("../data/batchqc/sizes/batchqc-2000.rds")
n <- ncol(X)
ncond <- n / 4
batch <- as.factor(rep(1:2, each = ncond * 2))
class <- rep(rep(LETTERS[1:2], each = ncond), 2)
metadata <- data.frame(batch, class, row.names = colnames(X))

villani <- readRDS("../data/scrna/villani/processed/villani.rds")
assay(villani, "tpm_sparse") <- Matrix(tpm(villani))


test_that("RVP.matrix gives correct value.", {
  expect_equal(
    RVP(t(X), batch, class),
    0.0442705
  )
})

test_that("RVP.data.frame gives correct value.", {
  expect_equal(
    RVP(data.frame(t(X)), batch, class),
    0.0442705
  )
})

test_that("RVP.SingleCellExperiment on dGCMatrix gives correct value.", {
  expect_equal(
    RVP.SingleCellExperiment(villani, "batch", "cell_type1"),
    RVP.SingleCellExperiment(villani, "batch", "cell_type1", assayname ="tpm_sparse")
  )
})
