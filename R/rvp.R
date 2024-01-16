#' Recursive variance partitioning (RVP)
#'
#' Calculates percentage of variance in data due to batch effects (S3 generic
#' function).
#'
#' @param X Dataframe or matrix with dim (n_samples, n_features).
#' @param batch Vector containing batch labels of samples.
#' @param cls Vector or list of vectors containing class labels of samples.
#' @return numeric indicating total percentage of variance in data due to batch effects.
rvp <- function(x, ...) UseMethod("rvp", x)


# TODO: Refactor permutation test as a standalone function
# TODO: Refactor rvp.default as an inner function .rvp.default
# TODO: Create private S3 generic .rvp which dispatches to .rvp.default and.rvp.sparseMatrix
# TODO: Change rvp.default to include permutation tests and to call .rvp

#' Recursive variance partitioning (RVP)
#'
#' Calculates percentage of variance in data due to batch effects. Default S3
#' method of the generic rvp function for data.frame or matrix classes.
#'
#' @param X Dataframe or matrix with dim (n_features, n_samples).
#' @param batch Vector containing batch labels of samples.
#' @param cls Vector or list of vectors containing class labels of samples.
#' @return numeric indicating total percentage of variance in data due to batch effects.
rvp.default <- function(X, batch, cls = NULL) {
  if (ncol(X) != length(batch)) {
    stop("Length of batch does not match number of columns in X!")
  }
  if (length(unique(batch)) == 1L) {
    # Use NA as is.na works on lists
    message("Only one batch present!")
    return(list(RVP = 0, sum.squares = NA)) # only one batch is present
  }
  X[is.na(X)] <- 0
  
  # COMPUTE RVP 
  if (is.null(cls)) {
    feature_means <- rowMeans(X)
    ss_total <- rowSums((X - feature_means) ^ 2)
    X_batches <- split_cols(X, batch, drop = TRUE)
    batch_means <- lapply(X_batches, rowMeans)
    batch_means <- do.call(cbind, batch_means)
    nperbatches <- sapply(X_batches, ncol)
    squares <- (batch_means - feature_means) ^ 2
    # Multiplying a matrix by a vector (broadcasted across rows) is equivalent
    # to multiplying a matrix by a diagonal matrix = diag(vector)
    ss_batch <- rowSums(squares %*% diag(nperbatches))
    stopifnot(length(ss_batch) == nrow(X))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- cbind(ss_batch, ss_total)
    rownames(SS) <- rownames(X)
    colnames(SS) <- c("ss_batch", "ss_total")
    # print("rvp.default: vars")
    # print(sapply(ls(), function(x) object_size(mget(x, inherits = TRUE))))
    return(list(RVP = pct_batch, sum.squares = SS))
  } else {
    feature_means <- rowMeans(X)
    ss_total <- rowSums((X - feature_means) ^ 2)
    X_classes <- split_cols(X, cls, drop = TRUE)
    batch_classes <- split(batch, cls, drop = TRUE)
    # Warning: recursive call
    objs <- mapply(
      rvp.default, X_classes, batch_classes,
      MoreArgs = list(cls = NULL),
      SIMPLIFY = FALSE
    )
    SS_classes <- lapply(objs, function(obj) obj$sum.squares)
    # Filters out obj$sum.squares == NA
    SS_classes <- SS_classes[!is.na(SS_classes)]
    if (length(SS_classes) == 0L) {
      confound_message <- paste(
        "RVP is unable to quantify batch effects as batch and class",
        "are completely confounded!"
      )
      stop(confound_message)
    }
    ss_batch_classes <- lapply(SS_classes, function(X) X[, "ss_batch"])
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- cbind(ss_batch, ss_total)
    rownames(SS) <- rownames(X)
    colnames(SS) <- c("ss_batch", "ss_total")

    return(list(RVP = pct_batch, sum.squares = SS))
  }
}


#' Recursive variance partitioning (RVP)
#'
#' Calculates percentage of variance in data due to batch effects. Default S3
#' method of the generic rvp function for data.frame or matrix classes.
#'
#' @param X Dataframe or matrix with dim (n_features, n_samples).
#' @param batch Vector containing batch labels of samples.
#' @param cls Vector or list of vectors containing class labels of samples.
#' @return numeric indicating total percentage of variance in data due to batch effects.
#' @import Matrix
rvp.sparseMatrix <- function(X, batch, cls = NULL) {
  if (ncol(X) != length(batch)) {
    stop("Length of batch does not match number of columns in X!")
  }
  if (length(unique(batch)) == 1L) {
    # Use NA as is.na works on lists
    message("Only one batch present!")
    return(list(RVP = 0, sum.squares = NA)) # only one batch is present
  }
  X[is.na(X)] <- 0
  
  # COMPUTE RVP 
  if (is.null(cls)) {
    feature_means <- rowMeans(X, sparseResult = TRUE)
    ss_total <- rowSums((X - feature_means) ^ 2, sparseResult = TRUE)
    X_batches <- split_cols(X, batch, drop = TRUE)
    batch_means <- lapply(
      X_batches,
      function(X) as(rowMeans(X, sparseResult = TRUE), "sparseMatrix")
    )
    # N.B. cbind only works on sparseMatrix and not sparseVector
    batch_means <- do.call(cbind, batch_means)
    nperbatches <- sapply(X_batches, ncol)
    squares <- (batch_means - feature_means) ^ 2
    # Multiplying a matrix by a vector (broadcasted across rows) is equivalent
    # to multiplying a matrix by a diagonal matrix = diag(vector)
    ss_batch <- rowSums(
      squares %*% .sparseDiagonal(x = nperbatches),
      sparseResult = TRUE
    )
    stopifnot(length(ss_batch) == nrow(X))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- Matrix(
      cbind(as(ss_batch, "sparseMatrix"), as(ss_total, "sparseMatrix")),
      dimnames = list(rownames(X), c("ss_batch", "ss_total")),
      sparse = TRUE
    )
    # print("rvp.sparseMatrix: vars")
    # print(sapply(ls(), function(x) object_size(mget(x, inherits = TRUE))))
    return(list(RVP = pct_batch, sum.squares = SS))
  } else {
    feature_means <- rowMeans(X, sparseResult = TRUE)
    ss_total <- rowSums((X - feature_means) ^ 2, sparseResult = TRUE)
    X_classes <- split_cols(X, cls, drop = TRUE)
    batch_classes <- split(batch, cls, drop = TRUE)
    # Warning: recursive call
    objs <- mapply(
      rvp.sparseMatrix, X_classes, batch_classes,
      MoreArgs = list(cls = NULL),
      SIMPLIFY = FALSE
    )
    SS_classes <- lapply(objs, function(obj) obj$sum.squares)
    # Filters out obj$sum.squares == NA
    SS_classes <- SS_classes[!is.na(SS_classes)]
    if (length(SS_classes) == 0L) {
      confound_message <- paste(
        "RVP is unable to quantify batch effects as batch and class",
        "are completely confounded!"
      )
      stop(confound_message)
    }
    ss_batch_classes <- lapply(
      SS_classes,
      # N.B. drop = False to return sparseMatrix instead of dense vector
      function(X) X[, "ss_batch", drop = FALSE]
    )
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    # TODO: Check if there is a problem with ss_batch being a sparseMatrix 
    # and ss_total being a sparseVector
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- Matrix(
      cbind(as(ss_batch, "sparseMatrix"), as(ss_total, "sparseMatrix")),
      dimnames = list(rownames(X), c("ss_batch", "ss_total")),
      sparse = TRUE
    )

    return(list(RVP = pct_batch, sum.squares = SS))
  }
}


#' Calculates percentage of variance in data due to batch effects
#'
#' @param sce SingleCellExperiment object
#' @param batchname Character vector of column name of colData representing batch information.
#' @param classname Character vector of column name/s of colData representing class information.
#' @param dataname Character vector of assay name of SCE object. By default
#'   the first assay is used.
#' @return numeric indicating total percentage of variance in data due to batch effects.
#' @import SingleCellExperiment, Matrix
rvp.SingleCellExperiment <- function(
  sce, batchname, classname, dataname = NULL
) {
  X <- if (is.null(dataname)) {
    assay(sce)
  } else {
    assay(sce, dataname)
  }
  batch <- sce[[batchname]]
  # TODO: Handle classname == NULL
  cls <- if (length(classname) > 1) {
    lapply(classname, function(name) sce[[name]])
  } else {
    sce[[classname]]
  }
  return(rvp(X, batch, cls))
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
#' @return numeric indicating total percentage of variance in data due to batch effects.
#' @export
rvp.Seurat <- function(
  obj, batchname, classname,
  layer = "data", nperm = 0
) {
  # Enhances: Matrix, Seurat/SeuratObject
  # GetAssayData is for Seurat assay v3/v4
  # GetAssayData is either from the Seurat or SeuratObject package
  X <- GetAssayData(obj, layer = layer)
  batch <- obj@meta.data[[batchname]]
  # TODO: Handle classname == NULL
  cls <- if (length(classname) > 1) {
    # obj[[name]] returns dataframe instead of vector!
    lapply(classname, function(name) obj@meta.data[[name]])
  } else {
    obj@meta.data[[classname]]
  }
  rvp_obj <- rvp(X, batch, cls)
  if (nperm > 0) {
    # TODO: Multiprocessing for permtest
    stopifnot(nperm >= 100)
    pb <- progress::progress_bar$new(
      format = "Permutation tests: [:bar] :current/:total in :elapsed.",
      total = nperm, clear = FALSE, width = 75
    )
    capture.output(pb$tick(0), file = nullfile())
    null_distr <- numeric()
    for (i in seq_len(nperm)) {
      shuffled_batch <- sample(batch)
      null_pct <- rvp(X, shuffled_batch, cls)
      null_distr <- c(null_distr, null_pct)
      pb$tick()
    }
    pct_batch <- rvp_obj$RVP
    pvalue <- sum(null_distr > pct_batch) / nperm
    rvp_obj$null.distribution <- null_distr
    rvp_obj$p.value <- pvalue
  }
  if (nperm == 0) {
    return(rvp_obj$RVP)
  } else {
    return(rvp_obj)
  }
}


#' Plots graph of truncated RVP at different feature lengths
#'
#' @param sum_sq Dataframe of sum of squares with dim (n_features, 2).
#' @param m Number of features to retain from sum_sq (from the top).
plot_rvp <- function(rvp_obj, m = NULL, cex = 1) {
  xlab <- "Feature index"
  sum_sq <- rvp_obj$sum.squares
  sum_sq <- sum_sq[rev(order(sum_sq$ss_total)), ]
  RVP <- rvp_obj$RVP

  if (is.numeric(m))
    sum_sq <- sum_sq[seq_len(m), , drop = FALSE]

  ax_ssb <- ggplot(sum_sq) +
    geom_line(
      aes(x = seq_len(nrow(sum_sq)), y = cumsum(ss_batch)),
      color = "blue"
    ) +
    labs(x = xlab, y = "Cumulative sum of squares batch")

  ax_sst <- ggplot(sum_sq) +
    geom_line(
      aes(x = seq_len(nrow(sum_sq)), y = cumsum(ss_total)),
      color = "chartreuse3"
    ) +
    labs(x = xlab, y = "Cumulative sum of squares total")

  ax_rvp <- ggplot(sum_sq) +
    geom_line(
      aes(x = seq_len(nrow(sum_sq)), y = cumsum(ss_batch) / cumsum(ss_total)),
      col = "brown1"
    ) +
    geom_hline(yintercept = RVP, col = "darkgray", linetype = "dashed") +
    labs(x = xlab, y = "Cumulative RVP")

  full_rvp <- ax_rvp +
    ylim(c(0, 1)) +
    theme(
      axis.title.x = element_blank(), 
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA)
    )

  ax <- ggdraw() +
    draw_plot(ax_rvp) +
    draw_plot(
      full_rvp,
      x = 0.5, y = 0.35,
      width = 0.45, height = 0.4
    )

  plot_grid(ax_ssb, ax_sst, ax, nrow = 1)
}


#' Divides array-like objects according to their columns
#'
#' @param x array-like object to be divided
#' @param f vector or list of vectors indicating the grouping of columns
#' @param drop logical indicating if levels that do not occur should be dropped
split_cols <- function(x, f, drop = FALSE, ...) {
  if (is.list(f)) {
    stopifnot(all(sapply(f, length) == ncol(x)))
  } else {
    stopifnot(length(f) == ncol(x))
  }
  lapply(
    split(seq_len(ncol(x)), f, drop, ...),
    function(ind) x[, ind, drop = FALSE]
  )
}
