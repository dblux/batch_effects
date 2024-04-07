set.seed(0)
library(testthat)
library(Matrix)
library(SingleCellExperiment)
library(SeuratObject)
source("R/rvp.R")

########## SETUP ########## 

# Data
n <- 200
m <- 1000
sids <-  sprintf("SID%03d", seq_len(n)) # n is 3 digits
x <- rgamma(m * n, 40, 5)
missing <- rbinom(m * n, 1, 0.8) == 1
x[missing] <- 0
dense <- matrix(
  x, m, n,
  dimnames = list(seq_len(m), sids)
)
sparse <- as(dense, "sparseMatrix")
data <- data.frame(dense)

# Metadata
ncond <- n / 4
batch <- rep(1:2, each = ncond * 2)
class1 <- rep(rep(LETTERS[1:2], each = ncond), 2)
class1_int <- ifelse(class1 == "A", 1, 2)
class2 <- rep(LETTERS[3:4], ncond * 2)
class3 <- rep(LETTERS[3:4], each = ncond * 2)
metadata <- data.frame(batch, class1, class2, row.names = colnames(dense))

sce <- SingleCellExperiment(
  assays = list(logcounts = sparse),
  colData = metadata
)
seu <- CreateSeuratObject(
  counts = sparse,
  meta.data = metadata
)

########## UNIT TESTS ##########

# Test: Utility functions
test_that("split_cols works when arg: f is a list of vectors", {
  list_dense <- split_cols(dense, list(class1, class2), drop = TRUE)

  expect_named(list_dense, c("A.C", "B.C", "A.D", "B.D"))
  expect_identical(
    unname(sapply(list_dense, nrow)),
    rep(nrow(dense), 4)
  )
})

# Test: RVP
RVP_TRUE <- 0.01011387
test_that("RVP works on a dense matrix", {
  # test that arg: batch, cls can be of type {numeric, character, factor}
  expect_equal(RVP(dense, batch, class1)$RVP, RVP_TRUE)
  expect_equal(RVP(dense, as.factor(batch), class1)$RVP, RVP_TRUE)
  expect_equal(RVP(dense, as.character(batch), class1)$RVP, RVP_TRUE)
  expect_equal(RVP(dense, batch, as.factor(class1))$RVP, RVP_TRUE)
  expect_equal(RVP(dense, batch, class1_int)$RVP, RVP_TRUE)
})

test_that("RVP works when arg: f is a list of vectors", {
  expect_no_error(RVP(dense, batch, list(class1, class2)))
})

test_that("RVP throws error when batch and class are completely confounded", {
  expect_error(
    RVP(dense, batch, list(class1, class3)),
    "batch and class are completely confounded"
  )
})

test_that("RVP works on a sparse matrix", {
  expect_equal(RVP(sparse, batch, class1)$RVP, RVP_TRUE)
})

test_that("RVP works on a data frame", {
  expect_equal(RVP(data, batch, class1)$RVP, RVP_TRUE)
})

test_that("RVP works on a SingleCellExperiment object", {
  expect_equal(RVP(sce, "batch", "class1")$RVP, RVP_TRUE)
})

test_that("RVP works on a Seurat object", {
  expect_equal(RVP(seu, "batch", "class1")$RVP, RVP_TRUE)
})

PVALUE <- 0.404
cat("Running permutation test...", fill = TRUE)
test_that("Permutation test returns p.value and null.distribution", {
  expect_no_error({
    obj <- RVP(dense, batch, class1, nperm = 1000)
    expect_equal(obj$p.value, PVALUE)
  })
})

test_that("Helper RVP function that uses sparse matrices works", {
  expect_equal(RVP(sparse, batch, class1, use.sparse = TRUE)$RVP, RVP_TRUE)
})
