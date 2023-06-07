#' Calculates percentage of variance in data due to batch effects
#' 
#' @param X Dataframe or matrix with dim (n_samples, n_features).
#' @param batch Vector containing batch labels of samples.
#' @param class Vector or list of vectors containing class labels of samples.
#' @param ret.obj Logical indicating whether to return object or percentage of variance.
#' @return numeric indicating total percentage of variance in data due to batch effects.
RVP <- function(X, batch, class = NULL, ret.obj = FALSE) {
  if (nrow(X) != length(batch))
    stop("Length of batch does not match number of rows in X!")
  if (is.vector(class) || is.factor(class))
    class <- as.character(class)
  if (length(unique(batch)) == 1) {
    # Use NA as is.na works on lists
    return(list(percentage = 0, sum_squares = NA)) # only one batch is present
  }

  X[is.na(X)] <- 0
  batch <- as.character(batch)
  
  # COMPUTE RVP 
  if (is.null(class)) {
    feature_means <- colMeans(X)
    ss_total <- colSums(sweep(X, 2, feature_means, `-`) ^ 2)
    X_batches <- split.data.frame(X, batch)
    # rm(X)
    batch_means <- sapply(X_batches, function(X) colMeans(X))
    nperbatches <- sapply(X_batches, nrow)
    squares <- (batch_means - feature_means) ^ 2
    ss_batch <- rowSums(sweep(squares, 2, nperbatches, `*`))

    stopifnot(length(ss_batch) == ncol(X))
    total_percentage <- sum(ss_batch) / sum(ss_total) 
    if (ret.obj) {
      return(list(
        percentage = total_percentage,
        sum_squares = data.frame(ss_batch, ss_total, row.names = colnames(X))
      ))
    } else {
      return(total_percentage)
    }
  } else {
    ss_total <- colSums(sweep(X, 2, colMeans(X), `-`) ^ 2)
    # If class is list, splits by all vectors in list 
    X_classes <- split.data.frame(X, class)
    # rm(X)
    X_classes <- Filter(function(X) nrow(X) != 0, X_classes)
    classes_string <- do.call(paste, as.list(names(X_classes)))
    # message(sprintf("Split into classes: %s", classes_string))
    batch_classes <- split(batch, class)
    batch_classes <- Filter(function(x) length(x) != 0, batch_classes)
    # Warning: recursive call
    objs <- mapply(
      RVP, X_classes, batch_classes,
      MoreArgs = list(class = NULL, ret.obj = TRUE),
      SIMPLIFY = FALSE
    )
    sumsquares_classes <- lapply(objs, function(obj) obj$sum_squares)
    # Filters out obj$sum_squares == NA
    sumsquares_classes <- sumsquares_classes[!is.na(sumsquares_classes)]
    ss_batch_classes <- lapply(sumsquares_classes, function(X) X$ss_batch)
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    stopifnot(length(ss_batch) == length(ss_total))
    total_percentage <- sum(ss_batch) / sum(ss_total)
    if (ret.obj) {
      return(list(
        percentage = total_percentage,
        percentage_classes = sapply(objs, function(obj) obj$percentage),
        sum_squares = data.frame(ss_batch, ss_total, row.names = colnames(X))
      ))
    } else {
      return(total_percentage)
    }
  }
}


#' Calculates percentage of variance in data due to batch effects
#' 
#' @param X Dataframe with dim (n_samples, n_features).
#' @param batchname Character vector of column name of colData representing batch information.
#' @param classname Character vector of column name/s of colData representing class information.
#' @param assayname Character vector of assay name of SCE object. By default
#'   the first assay is used.
#' @param ret.obj Logical indicating whether to return object or percentage of variance.
#' @return numeric indicating total percentage of variance in data due to batch effects.
RVP.SingleCellExperiment <- function(
  sce, batchname, classname, assayname = NULL, ret.obj = FALSE
) {
  X <- if (is.null(assayname)) {
    t(assay(sce))
  } else {
    t(assay(sce, assayname))
  }
  batch <- colData(sce)[[batchname]]
  class <- if (length(classname) > 1) {
    lapply(classname, function(name) colData(sce)[[name]])
  } else {
    colData(sce)[[classname]]
  }
  
  return(RVP(X, batch, class, ret.obj))
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
    labs(x = xlab, y = "Cumulative sum")
  
  ax2 <- ggplot(sum_sq) +
    geom_line(
      aes(x = seq_len(nrow(sum_sq)), y = cumsum(ss_total))
    ) +
    labs(x = xlab, y = "Cumulative sum")

  ax3 <- ggplot(sum_sq) +
    geom_line(
      aes(
        x = seq_len(nrow(sum_sq)),
        y = cumsum(ss_batch) / cumsum(ss_total)
      )
    ) +
    labs(x = xlab, y = "RVP")

  ax4 <- ggplot(sum_sq) +
    geom_point(
      aes(x = seq_len(nrow(sum_sq)), y = ss_batch),
      color = 'blue', cex = cex
    ) + 
    geom_point(
      aes(x = seq_len(nrow(sum_sq)), y = ss_total),
      cex = cex
    ) +
    labs(x = xlab, y = "Sum of squares")

  cowplot::plot_grid(ax1, ax2, ax3, ax4, nrow = 2)
}
