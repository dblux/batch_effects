library(testthat)
library(Seurat)
library(SingleCellExperiment)
library(Matrix)
source("R/rvp.R")


# SETUP
X <- readRDS("data/batchqc/sizes/batchqc-2000.rds")
n <- ncol(X)
ncond <- n / 4
batch <- as.factor(rep(1:2, each = ncond * 2))
class <- rep(rep(LETTERS[1:2], each = ncond), 2)
metadata <- data.frame(batch, class, row.names = colnames(X))

villani <- readRDS("data/scrna/villani/processed/villani.rds")
assay(villani, "tpm_sparse") <- Matrix(tpm(villani))


# UNIT TESTS
test_that("rvp.matrix gives correct value.", {
  expect_equal(rvp(t(X), batch, class), 0.0442705)
})

test_that("rvp.data.frame gives correct value.", {
  expect_equal(rvp(data.frame(t(X)), batch, class), 0.0442705)
})

test_that("rvp.SingleCellExperiment on dGCMatrix gives correct value.", {
  expect_equal(
    rvp.SingleCellExperiment(villani, "batch", "cell_type1"),
    rvp.SingleCellExperiment(
      villani, "batch", "cell_type1", assayname = "tpm_sparse"
    )
  )
})