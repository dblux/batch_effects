library(testthat)
source("R/subset.R")


# SETUP
X <- matrix(1, 10, 10)
X[1, ] <- rep(0, 10)
X[2, ] <- c(1, rep(0, 9))
X[3, ] <- rep(c(1, rep(0, 4)), 2)
X[4, ] <- rep(c(0, 1), each = 5)
data <- data.frame(X)

class <- rep(1:2, each = 5)


# TESTS
test_that("remove_sparse works for matrix.", {
  expect_equal(nrow(remove_sparse(X, 0.9)), 8) # row 1, 2 removed
  expect_equal(nrow(remove_sparse(X, 0.9, class)), 9) # row 1 removed
  expect_equal(nrow(remove_sparse(X, 0.9, class, any)), 7) # row 1, 2, 4 removed
})

test_that("remove_sparse works for data.frame.", {
  expect_equal(nrow(remove_sparse(data, 0.9)), 8) # row 1, 2 removed
  expect_equal(nrow(remove_sparse(data, 0.9, class)), 9) # row 1 removed
  expect_equal(nrow(remove_sparse(data, 0.9, class, any)), 7) # 1, 2, 4 removed
})
