#' Recursive variance partitioning (RVP)
#'
#' Calculates percentage of variance in data due to batch effects. It is a
#' generic function that invokes convenience S3 methods for the
#' SingleCellExperiment and Seurat classes.
#'
#' @param X Dataframe or matrix with dim (n_samples, n_features).
#' @param batch Vector containing batch labels of samples.
#' @param cls Vector or list of vectors containing class labels of samples.
#' @param ret.percent Logical indicating whether to only return percentage of
#'   variance due to batch.
#' @return numeric indicating total percentage of variance in data due to batch effects.
rvp <- function(x, ...) UseMethod("rvp", x)


#' Recursive variance partitioning (RVP)
#'
#' Calculates percentage of variance in data due to batch effects. Default S3
#' method of the generic rvp function for data.frame or matrix classes.
#'
#' @param X Dataframe or matrix with dim (n_samples, n_features).
#' @param batch Vector containing batch labels of samples.
#' @param cls Vector or list of vectors containing class labels of samples.
#' @param ret.percent Logical indicating whether to only return percentage of
#'   variance due to batch.
#' @return numeric indicating total percentage of variance in data due to batch effects.
rvp.default <- function(X, batch, cls = NULL, ret.percent = TRUE) {
  if (nrow(X) != length(batch))
    stop("Length of batch does not match number of rows in X!")
  if (is.vector(cls) || is.factor(cls))
    cls <- as.character(cls)
  if (length(unique(batch)) == 1) {
    # Use NA as is.na works on lists
    return(list(percent.batch = 0, sum.squares = NA)) # only one batch is present
  }

  X[is.na(X)] <- 0
  batch <- as.character(batch)
  
  # COMPUTE RVP 
  if (is.null(cls)) {
    feature_means <- colMeans(X)
    ss_total <- colSums(sweep(X, 2, feature_means, `-`) ^ 2)
    X_batches <- split.data.frame(X, batch)
    # rm(X)
    batch_means <- sapply(X_batches, function(X) colMeans(X))
    nperbatches <- sapply(X_batches, nrow)
    squares <- (batch_means - feature_means) ^ 2
    ss_batch <- rowSums(sweep(squares, 2, nperbatches, `*`))

    stopifnot(length(ss_batch) == ncol(X))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    if (ret.percent) {
      return(pct_batch)
    } else {
      return(list(
        percent.batch = pct_batch,
        sum.squares = data.frame(ss_batch, ss_total, row.names = colnames(X))
      ))
    }
  } else {
    ss_total <- colSums(sweep(X, 2, colMeans(X), `-`) ^ 2)
    # If cls is list, splits by all vectors in list
    X_classes <- split.data.frame(X, cls)
    # rm(X)
    X_classes <- Filter(function(X) nrow(X) != 0, X_classes)
    classes_string <- do.call(paste, as.list(names(X_classes)))
    # message(sprintf("Split into classes: %s", classes_string))
    batch_classes <- split(batch, cls)
    batch_classes <- Filter(function(x) length(x) != 0, batch_classes)
    # Warning: recursive call
    objs <- mapply(
      rvp.default, X_classes, batch_classes,
      MoreArgs = list(cls = NULL, ret.percent = FALSE),
      SIMPLIFY = FALSE
    )
    sumsquares_classes <- lapply(objs, function(obj) obj$sum.squares)
    # Filters out obj$sum.squares == NA
    sumsquares_classes <- sumsquares_classes[!is.na(sumsquares_classes)]
    ss_batch_classes <- lapply(sumsquares_classes, function(X) X$ss_batch)
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    if (ret.percent) {
      return(pct_batch)
    } else {
      return(list(
        percent.batch = pct_batch,
        sum.squares = data.frame(ss_batch, ss_total, row.names = colnames(X))
      ))
    }
  }
}


#' Calculates percentage of variance in data due to batch effects
#'
#' @param sce SingleCellExperiment object
#' @param batchname Character vector of column name of colData representing batch information.
#' @param classname Character vector of column name/s of colData representing class information.
#' @param dataname Character vector of assay name of SCE object. By default
#'   the first assay is used.
#' @param ret.percent Logical indicating whether to only return percentage of
#'   variance due to batch.
#' @return numeric indicating total percentage of variance in data due to batch effects.
#' @import SingleCellExperiment, Matrix
rvp.SingleCellExperiment <- function(
  sce, batchname, classname,
  dataname = NULL, ret.percent = TRUE
) {
  X <- if (is.null(dataname)) {
    Matrix::t(assay(sce))
  } else {
    Matrix::t(assay(sce, dataname))
  }
  batch <- sce[[batchname]]
  # TODO: Handle classname == NULL
  cls <- if (length(classname) > 1) {
    lapply(classname, function(name) sce[[name]])
  } else {
    sce[[classname]]
  }
  return(rvp.default(X, batch, cls, ret.percent))
}


#' Calculates percentage of variance in data due to batch effects
#'
#' @param obj Seurat object
#' @param batchname Character vector of column name of colData representing batch information.
#' @param classname Character vector of column name/s of colData representing class information.
#' @param slot Character indicating slot in assay to use. E.g. counts
#' @param nperm Number of permutations to simulate in the Monte Carlo
#'   permutation test. A mininum value of 100 is required. By default, no
#'   permutation testing is performed.
#' @param ret.percent Logical indicating whether to only return percentage of
#'   variance due to batch.
#' @return numeric indicating total percentage of variance in data due to batch effects.
#' @import Seurat, Matrix
#' @export
rvp.Seurat <- function(
  obj, batchname, classname,
  slot = "data", nperm = NULL, ret.percent = TRUE
) {
  X <- if (length(find("LayerData")) == 0) {
    Matrix::t(Seurat::GetAssayData(obj, slot))
  } else {
    Matrix::t(Seurat::LayerData(obj, slot)) # Seurat v5
  }
  batch <- obj@meta.data[[batchname]]
  # TODO: Handle classname == NULL
  cls <- if (length(classname) > 1) {
    # obj[[name]] returns dataframe instead of vector!
    lapply(classname, function(name) obj@meta.data[[name]])
  } else {
    obj@meta.data[[classname]]
  }
  rvp_obj <- rvp.default(X, batch, cls, ret.percent = FALSE)
  if (is.numeric(nperm)) {
    stopifnot(nperm >= 100)
    null_pcts <- numeric()
    for (i in seq_len(nperm)) {
      shuffled_batch <- sample(batch)
      null_pct <- rvp.default(X, shuffled_batch, cls, ret.percent = TRUE)
      null_pcts <- c(null_pcts, null_pct)
    }
    pct_batch <- rvp_obj$percent.batch
    pvalue <- sum(null_pcts > pct_batch) / nperm
    rvp_obj$null.percentages <- null_pcts
    rvp_obj$p.value <- pvalue
  }
  if (ret.percent) {
    return(rvp_obj$percent.batch)
  } else {
    return(rvp_obj)
  }
}


#' Plots graph of truncated RVP at different feature lengths
#'
#' @param sum_sq Dataframe of sum of squares with dim (n_features, 2).
#' @param m Number of features to retain from sum_sq (from the top).
plot.rvp <- function(sum_sq, m = NULL, cex = 1) {
  xlab <- "Feature index"
  if (is.numeric(m))
    sum_sq <- sum_sq[seq_len(m), , drop = FALSE]

  ax1 <- ggplot(sum_sq) +
    geom_line(
      aes(x = seq_len(nrow(sum_sq)), y = cumsum(ss_batch)),
      color = 'blue'
    ) +
    geom_line(
      aes(x = seq_len(nrow(sum_sq)), y = cumsum(ss_total))
    ) +
    labs(x = xlab, y = "Cumulative sum of squares")
  
  ax2 <- ggplot(sum_sq) +
    geom_line(
      aes(x = seq_len(nrow(sum_sq)), y = cumsum(ss_total))
    ) +
    labs(x = xlab, y = "Cumulative sum of squares (total)")

  ax3 <- ggplot(sum_sq) +
    geom_line(
      aes(
        x = seq_len(nrow(sum_sq)),
        y = cumsum(ss_batch) / cumsum(ss_total)
      )
    ) +
    labs(x = xlab, y = "RVP")

  # ax4 <- ggplot(sum_sq) +
  #   geom_point(
  #     aes(x = seq_len(nrow(sum_sq)), y = ss_batch),
  #     color = 'blue', cex = cex
  #   ) + 
  #   geom_point(
  #     aes(x = seq_len(nrow(sum_sq)), y = ss_total),
  #     cex = cex
  #   ) +
  #   labs(x = xlab, y = "Sum of squares")

  cowplot::plot_grid(ax1, ax3, nrow = 1)
}
