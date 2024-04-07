#' Recursive variance partitioning (RVP)
#'
#' `RVP` calculates the proportion of variance associated with batch effects
#' in a data set ("RVP" value of a data set). To determine whether batch
#' effects are statistically significant in a data set, a permutation test can
#' be performed by setting `nperm` to a number above 100. `RVP` is
#' an S3 generic function; methods can be added for new classes. S3 methods
#' for class: array-like objects (default), `SummarizedExperiment`, 
#' `SingleCellExperiment` and `Seurat` are provided.
#'
#' @param x object to calculate RVP for.
#' @param ... additional arguments to pass to S3 methods.
#'
#' @returns List containing the following components:
#'   \describe{
#'     \item{RVP}{the proportion of variance associated with batch effects.}
#'     \item{`sum.squares`}{matrix of sum of squares between batch and total
#'       sum of squares for all features.}
#'     \item{`p.value`}{p-value of permutation test}
#'     \item{`null.distribution`}{numeric, null distribution of RVP values.}
#'   }
#'
#' @export
#'
RVP <- function(x, ...) UseMethod("RVP", x)


#' @details
#' `RVP.default()` assumes that `X` is of class array-like. It is
#' invoked by other S3 methods and contains logic that performs the
#' permutation test if necessary and allows users to specify which helper
#' RVP function to use.
#'
#' @param X data frame or matrix with dim (nfeatures, nsamples).
#' @param batch vector, indicating the batch information of samples.
#' @param cls vector or list of vectors with class information of samples.
#' @param nperm numeric indicating number of permutations to simulate in the
#'   Monte Carlo permutation test. We recommend a value no less than 1000.
#'   By default, no permutation test is performed.
#' @param use.sparse logical indicating whether to use sparse matrices when
#'   computing RVP. N.B. Using sparse matrices may lead to slight increase
#'  in run time. 
#'
#' @rdname RVP
#' @export
#'
RVP.default <- function(
  X, batch, cls = NULL,
  nperm = 0, use.sparse = FALSE
) {
  if (!use.sparse) {
    obj <- .RVP(X, batch, cls)
  } else {
    obj <- .RVP_sparseMatrix(X, batch, cls)
  }

  # Permutation test
  if (nperm > 0L) {
    # TODO: Multiprocessing for permtest
    if (nperm < 100)
      stop("nperm has to be above 100!")

    pb <- progress::progress_bar$new(
      format = "Permutations: [:bar] :current/:total in :elapsed.",
      total = nperm, clear = FALSE, width = 75
    )
    capture.output(pb$tick(0), file = nullfile())

    null_distr <- numeric()
    for (i in seq_len(nperm)) {
      shuffled_batch <- sample(batch)
      if (!use.sparse) {
        null_rvp <- .RVP(X, shuffled_batch, cls)$RVP
      } else {
        null_rvp <- .RVP_sparseMatrix(X, shuffled_batch, cls)$RVP
      }
      null_distr <- c(null_distr, null_rvp)
      pb$tick()
    }
    obj$p.value <- sum(null_distr > obj$RVP) / nperm
    obj$null.distribution <- null_distr
  }
  obj
}


#' @param seu Seurat object.
#' @param batchname character, name of column in metadata indicating batch.
#' @param classname character, name of column/s in metadata indicating class. 
#' @param slot character, name of slot in assay to use. E.g. counts.
#' @param ... optional arguments to pass to `RVP.default`
#'
#' @rdname RVP
#' @export
#'
RVP.Seurat <- function(
  seu, batchname, classname = NULL,
  slot = "data",
  ...
) {
  # Suggests: SeuratObject
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("Please install SeuratObject package!")

  # GetAssayData is for Seurat assay v3/v4
  X <- SeuratObject::GetAssayData(seu, slot = slot)
  batch <- seu@meta.data[[batchname]]
  cls <- if (is.null(classname) || is.na(classname)) {
    NULL
  } else if (length(classname) == 1L) {
    # seu[[name]] returns dataframe instead of vector!
    seu@meta.data[[classname]]
  } else {
    lapply(classname, function(name) seu@meta.data[[name]])
  }
  RVP.default(X, batch, cls, ...)
}


#' @param se `SummarizedExperiment`/`SingleCellExperiment` object.
#'   `SingleExperiment` class inherits from the `SummarizedExperiment` class. 
#' @param assayname character, name of assay to use. By default the first
#'   assay is used.
#'
#' @rdname RVP
#' @export
#'
RVP.SummarizedExperiment <- function(
  se, batchname, classname = NULL,
  assayname = NULL,
  ...
) {
  # Suggests: SummarizedExperiment 
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("Please install SummarizedExperiment package!")

  X <- if (is.null(assayname) || is.na(classname)) {
    SummarizedExperiment::assay(se)
  } else {
    SummarizedExperiment::assay(se, assayname)
  }
  batch <- se[[batchname]]
  cls <- if (is.null(classname) || is.na(classname)) {
    NULL
  } else if (length(classname) == 1L) {
    se[[classname]]
  } else {
    lapply(classname, function(name) se[[name]])
  }
  RVP.default(X, batch, cls, ...)
}


#' Helper function that calculates RVP
#'
#' @param X data frame or matrix with dim (nfeatures, nsamples).
#' @param batch vector, indicating the batch information of samples.
#' @param cls vector or list of vectors with class information of samples.
#'
.RVP <- function(X, batch, cls) {
  ### CHECK ARGS ###
  if (ncol(X) != length(batch))
    stop("Length of batch does not match number of columns in X!")

  if (length(unique(batch)) == 1L) {
    # Use NA as is.na works on lists
    message("Only one batch present!")
    return(list(RVP = 0, sum.squares = NA)) # only one batch is present
  }
  X[is.na(X)] <- 0

  ### COMPUTE RVP ###
  if (is.null(cls)) {
    feature_names <- rownames(X)
    feature_means <- rowMeans(X)
    ss_total <- rowSums((X - feature_means) ^ 2)
    X_batches <- split_cols(X, batch, drop = TRUE)
    rm(X)
    batch_means <- lapply(X_batches, rowMeans)
    batch_means <- do.call(cbind, batch_means)
    nperbatches <- sapply(X_batches, ncol)
    rm(X_batches)
    # Multiplying a matrix by a vector (broadcasted across rows) is equivalent
    # to multiplying a matrix by a diagonal matrix = diag(vector)
    ss_batch <- rowSums(
      (batch_means - feature_means) ^ 2 %*% diag(nperbatches)
    )
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- cbind(ss_batch, ss_total)
    rownames(SS) <- feature_names 
    colnames(SS) <- c("ss_batch", "ss_total")

    return(list(RVP = pct_batch, sum.squares = SS))
  } else {
    feature_names <- rownames(X)
    ss_total <- rowSums((X - rowMeans(X)) ^ 2)
    X_classes <- split_cols(X, cls, drop = TRUE)
    rm(X)
    batch_classes <- split(batch, cls, drop = TRUE)
    # Warning: Recursive call
    SS_classes <- lapply(
      mapply(
        .RVP, X_classes, batch_classes,
        MoreArgs = list(cls = NULL),
        SIMPLIFY = FALSE
      ),
      function(obj) obj$sum.squares
    )
    rm(X_classes)
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
    rm(SS_classes)
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- cbind(ss_batch, ss_total)
    rownames(SS) <- feature_names
    colnames(SS) <- c("ss_batch", "ss_total")

    return(list(RVP = pct_batch, sum.squares = SS))
  }
}


#' Helper function that calculates RVP using sparse matrices
#'
#' @param X data frame or matrix with dim (nfeatures, nsamples).
#' @param batch vector, indicating the batch information of samples.
#' @param cls vector or list of vectors with class information of samples.
#'
#' @import Matrix
#'
.RVP_sparseMatrix <- function(X, batch, cls = NULL) {
  ### CHECK ARGS ###
  if (ncol(X) != length(batch))
    stop("Length of batch does not match number of columns in X!")
  
  if (length(unique(batch)) == 1L) {
    # Use NA as is.na works on lists
    message("Only one batch present!")
    return(list(RVP = 0, sum.squares = NA)) # only one batch is present
  }
  X[is.na(X)] <- 0
  
  ### COMPUTE RVP ###
  if (is.null(cls)) {
    feature_names <- rownames(X)
    feature_means <- rowMeans(X, sparseResult = TRUE)
    ss_total <- rowSums((X - feature_means) ^ 2, sparseResult = TRUE)
    X_batches <- split_cols(X, batch, drop = TRUE)
    rm(X)
    batch_means <- lapply(
      X_batches,
      function(X) as(rowMeans(X, sparseResult = TRUE), "sparseMatrix")
    )
    # N.B. cbind only works on sparseMatrix and not sparseVector
    batch_means <- do.call(cbind, batch_means)
    nperbatches <- sapply(X_batches, ncol)
    rm(X_batches)
    # Multiplying a matrix by a vector (broadcasted across rows) is equivalent
    # to multiplying a matrix by a diagonal matrix = diag(vector)
    ss_batch <- rowSums(
      (batch_means - feature_means) ^ 2 %*% .sparseDiagonal(x = nperbatches),
      sparseResult = TRUE
    )
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- Matrix(
      cbind(as(ss_batch, "sparseMatrix"), as(ss_total, "sparseMatrix")),
      dimnames = list(feature_names, c("ss_batch", "ss_total")),
      sparse = TRUE
    )
    # print("RVP.sparseMatrix: vars")
    # print(sapply(ls(), function(x) object_size(mget(x, inherits = TRUE))))
    return(list(RVP = pct_batch, sum.squares = SS))
  } else {
    feature_names <- rownames(X)
    ss_total <- rowSums(
      (X - rowMeans(X, sparseResult = TRUE)) ^ 2,
      sparseResult = TRUE
    )
    X_classes <- split_cols(X, cls, drop = TRUE)
    rm(X)
    batch_classes <- split(batch, cls, drop = TRUE)
    # Warning: Recursive call
    SS_classes <- lapply(
      mapply(
        .RVP_sparseMatrix, X_classes, batch_classes,
        MoreArgs = list(cls = NULL),
        SIMPLIFY = FALSE
      ),
      function(obj) obj$sum.squares
    )
    rm(X_classes)
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
    rm(SS_classes)
    stopifnot(is.list(ss_batch_classes))
    ss_batch <- Reduce(`+`, ss_batch_classes)
    # TODO: Check if there is a problem with ss_batch being a sparseMatrix 
    # and ss_total being a sparseVector
    stopifnot(length(ss_batch) == length(ss_total))
    pct_batch <- sum(ss_batch) / sum(ss_total)
    SS <- Matrix(
      cbind(as(ss_batch, "sparseMatrix"), as(ss_total, "sparseMatrix")),
      dimnames = list(feature_names, c("ss_batch", "ss_total")),
      sparse = TRUE
    )

    return(list(RVP = pct_batch, sum.squares = SS))
  }
}


#' Plots graph of truncated RVP at different feature lengths
#'
#' @param sum_sq Dataframe of sum of squares with dim (n_features, 2).
#' @param m Number of features to retain from sum_sq (from the top).
#'
#' @returns Returns ggplot of cumulative sum of squares
#'
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


#' Splits subsettable objects according to their columns
#'
#' @param x subsettable object to be split
#' @param f vector or list of vectors indicating the grouping of columns
#' @param drop logical indicating if levels that do not occur should be dropped
#' @param ... optional arguments to \code{split}
#'
#' @returns List of objects split by columns
#'
#' @export
#'
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
