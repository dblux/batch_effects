#' Subset columns according to features in annotation
#'
#' Colnames of X has to match with rownames of annot.
#'
#' @example subset_cols(X, annot, subtype == 'BCR-ABL' & class_info == 'D0')
subset_cols <- function(X, metadata, ...) {
  X[, rownames(subset(metadata[colnames(X), , drop = FALSE], ...))]
}


#' Select top n highly variable features
#' @param X matrix or dataframe with rownames
select_hvg <- function(X, n, return.features = FALSE) {
  if (is.null(rownames(X)))
    stop("X does not have rownames")
  
  row_var <- apply(X, 1, var)
  sorted_var <- sort(row_var, decreasing = TRUE)
  hvg <- names(sorted_var)[1:n]
  
  if (return.features)
    return(hvg)
  
  X[hvg, ]
}


#' Select top n highest expressed features
select_highexpr <- function(X, n, return.features = FALSE) {
  stopifnot(!is.na(sum(X))) # assert no NA values
  feat_expr <- rowMeans(X)
  idx_ranked <- order(feat_expr, decreasing = TRUE)

  if (return.features)
    return(rownames(X)[idx_ranked[seq_len(n)]])

  X[idx_ranked[seq_len(n)], , drop = FALSE]
}


#' Remove sparse features with high percentage of zeros
#'
#' @param X Dataframe with dim (nfeature, nsample)
#' @param pct_zero Percentage threshold of zeros. Features with percentage of
#'   zeros above and equal to this threshold are removed.
#' @param class Class labels of samples.
#' @param func Function that is either "all" or "any". Features with
#'   all/any classes having a larger percentage of zeros than "pct_zero" are
#'   removed.
#' @param ret.features Logical indicating whether to return sparse features
#' @return Data frame without sparse features or vector of sparse features.
remove_sparse <- function(
  X, pct_zero, class = NULL, func = all, ret.features = FALSE
) {
  is_sparse <- function(X, pct_zero) rowSums(X == 0) / ncol(X) >= pct_zero
  if (is.null(class)) {
    idx <- is_sparse(X, pct_zero)
  } else {
    stopifnot(length(class) == ncol(X))
    X_classes <- split.default(data.frame(X), class)
    list_idx <- lapply(X_classes, is_sparse, pct_zero)
    idx <- do.call(mapply, c(func, list_idx))
  }
  if (ret.features) {
    return(rownames(X)[idx])
  } else {
    return(X[!idx, ])
  }
}


#' Remove rows of dataframe which meet the condition
#'
#' Special syntax: 'row' indicates row of dataframe `X`
#'
#' @examples
#' remove_rows(X, sum(row) == 0)
#' remove_rows(X, var(row) == 0)
remove_rows <- function(X, condition) {
  idx <- vector()
  for (i in seq_len(nrow(X))) {
    row <- as.numeric(X[i, ])
    # eval env defaults to calling env
    # i.e. row is automatically used for expression
    is_cond <- eval(substitute(condition))
    idx <- c(idx, is_cond)
  }
  
  X[!idx, ]
}


#' Removes ambiguous and AFFY probesets from dataframe
#' Rowname of affymetrix probesets
remove_affymetrix <- function(df) {
  logical_vec <- grepl("[0-9]_at", rownames(df)) &
    !startsWith(rownames(df), "AFFX")
  print(paste0("No. of ambiguous and AFFY probesets removed: ",
               nrow(df) - sum(logical_vec)))
  return(df[logical_vec, , drop=F])
}


#' Find similar genes to vector of gene symbols
#' @return named list with gene symbols in x as names and
#' matching gene symbols in y as values
get_similar_genes <- function(x, y) {
  # Default value of gene with no match is NULL
  matched_genes <- vector("list", length(x))
  names(matched_genes) <- x
  
  idx_notequiv <- !(x %in% y)
  not_equiv <- x[idx_notequiv]
  idx_mult <- grepl("///", not_equiv)
  
  matched_singles <- lapply(
    not_equiv[!idx_mult],
    grep, x = y, value = T
  )
  names(matched_singles) <- not_equiv[!idx_mult]
  matched_singles_fltr <- Filter(
    function(x) ifelse(length(x) == 0, F, T),
    matched_singles
  )

  mult_symbols <- lapply(
    not_equiv[idx_mult],
    function(x) unlist(strsplit(x, " /// "))
  )
  names(mult_symbols) <- not_equiv[idx_mult]
  matched_mults <- lapply(
    mult_symbols,
    function(x) unlist(lapply(x, grep, x = y, value = T))
  )
  matched_mults_fltr <- Filter(
    function(x) ifelse(length(x) == 0, F, T),
    matched_mults
  )
  
  # initialise list with equivalent genes
  matched_genes[!idx_notequiv] <- x[!idx_notequiv]
  # replace similar singles and multiples
  matched_genes[names(matched_singles_fltr)] <- matched_singles_fltr
  matched_genes[names(matched_mults_fltr)] <- matched_mults_fltr

  matched_genes
}