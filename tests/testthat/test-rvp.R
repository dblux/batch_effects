library(testthat)
library(SingleCellExperiment)
library(Matrix)
source("R/rvp.R")


### SETUP ###
X <- readRDS("data/batchqc/sizes/batchqc-2000.rds")
n <- ncol(X)
ncond <- n / 4
batch <- rep(1:2, each = ncond * 2)
class1 <- rep(rep(LETTERS[1:2], each = ncond), 2)
class1_numeric <- ifelse(class1 == "A", 1, 2)
class2 <- rep(LETTERS[3:4], ncond * 2)
class3 <- rep(LETTERS[3:4], each = ncond * 2)
metadata <- data.frame(batch, class1, class2, row.names = colnames(X))

villani <- readRDS("data/villani/processed/villani-sce.rds")
# Convert tpm assay from a dense to sparse matrix
tpm(villani) <- Matrix(tpm(villani))


## UNIT TESTS ###
# BatchQC dataset
test_that("split_cols works when arg: f is a list of vectors", {
  Xs <- split_cols(X, list(class1, class2), drop = TRUE)

  expect_named(Xs, c("A.C", "B.C", "A.D", "B.D"))
  expect_identical(unname(sapply(Xs, nrow)), rep(nrow(X), 4))
})

RVP_BATCHQC <- 0.0442705
test_that("rvp works on a dense matrix", {
  # test that arg: batch, cls can be of type {numeric, character, factor}
  expect_equal(rvp(X, batch, class1)$RVP, RVP_BATCHQC)
  expect_equal(rvp(X, as.factor(batch), class1)$RVP, RVP_BATCHQC)
  expect_equal(rvp(X, as.character(batch), class1)$RVP, RVP_BATCHQC)
  expect_equal(rvp(X, batch, as.factor(class1))$RVP, RVP_BATCHQC)
  expect_equal(rvp(X, batch, class1_numeric)$RVP, RVP_BATCHQC)
})

test_that("rvp works when arg: f is a list of vectors", {
  expect_no_error(
    rvp(X, batch, list(class1, class2))
  )
})

test_that("rvp throws an error when batch and cls are completely confounded", {
  expect_error(
    rvp(X, batch, list(class1, class3)),
    "batch and class are completely confounded"
  )
})

test_that("rvp works on a data frame", {
  expect_equal(rvp(data.frame(X), batch, class1)$RVP, RVP_BATCHQC)
})

# Villani dataset
RVP_VILLANI <- 0.08885395
test_that("rvp works on a sparse matrix", {
  expect_equal(
    rvp(tpm(villani), villani$batch, villani$cell_type1)$RVP,
    RVP_VILLANI 
  )
})

test_that("rvp.SingleCellExperiment S3 method works", {
  expect_equal(
    rvp.SingleCellExperiment(villani, "batch", "cell_type1")$RVP,
    RVP_VILLANI
  )
})
