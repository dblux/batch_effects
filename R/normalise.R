# Cosine normalisation / L2-norm normalisation
normalise_cosine <- function(df1) {
  l2norm_vec <- apply(df1, 2, lpnorm) 
  return(mapply(`/`, df1, l2norm_vec))
}


# Min-max scaling function
# Returns: Scaled vector with range [0,1]
normaliseMinmax <- function(vec) {(vec-min(vec))/(max(vec)-min(vec))}


#' Scale (using trimmed mean) values that have not been logged
#'
#' Trimmed mean scaling does not remove all tied values
#' @param X Data frame of numbers with dim features, samples 
#' @param target_mean Numeric indicating trimmed mean after normalisation
#' @param trim Upper and lower percentage of numbers to be trimmed
scale_trimmed <- function(X, target_mean = 500, trim = 0.02) {
  trimmed_mean <- apply(X, 2, mean, trim = trim)
  scaling_factor <- target_mean / trimmed_mean

  sweep(X, 2, scaling_factor, `*`)
}


# x <- data.frame(a = c(1:4, rep(0,6)),
#                 b = c(1:5, rep(0,5)),
#                 c = c(1:6, rep(0,4)))
# Quantile normalisation: 0 values are assigned 0 automatically
# Takes in df where columns are samples and rows are genes
normaliseQuantile <- function(df) {
  zero_filter <- df == 0
  sort_arr <- apply(df, 2, sort)
  # Creates reference distribution
  ref_distr <- apply(sort_arr, 1, mean)
  rank_arr <- apply(df, 2, rank, ties.method = "min")
  # print(rank_arr)
  qnorm_arr <- apply(rank_arr, c(1,2), function(x) ref_distr[x])
  rownames(qnorm_arr) <- rownames(df)
  qnorm_df <- as.data.frame(qnorm_arr)
  qnorm_df[zero_filter] <- 0
  return(qnorm_df)
}


# Problem: When there are too many zeros and fewer values are assigned to be 0
# .. No zeroes will be assigned
# Gene Fuzzy Scoring function transforms gene expression values
# Wilson Goh's paper
# Dense rank is used
normaliseGFS <- function(A, upper = 0.05, lower = 0.15, num_intervals = 0) {
  # Bins score with range [0,1] into intervals
  # E.g. 4 Intervals: Binned into 0.2, 0.4, 0.6, 0.8
  bin <- function(score, num_intervals) {
    for (i in 1:num_intervals) {
      if (score <= i/num_intervals) {
        return (i/(num_intervals+1))
      }
    }
  }
  cat(sprintf("Top %.2f of expressed genes are assigned GFS scores of 1\n", upper))
  cat(sprintf("Genes below the top %.2f of expressed genes are assigned GFS scores of 0\n", lower))
  # Rank function ranks largest value as 1 [-A is used]
  # Handle NaN?
  ranked_A <- apply(-A, 2, dplyr::dense_rank)
  rownames(ranked_A) <- rownames(A)
  # Returns [1,] = upper, [2,] = lower
  qtile <- apply(ranked_A, 2, quantile, probs=c(upper, lower), names=F)
  
  if (num_intervals <= 0) {
    for (c in 1:ncol(ranked_A)) {
      # Calculate qtile range
      q_range <- qtile[2,c] - qtile[1,c]
      for (r in 1:nrow(ranked_A)) {
        if (ranked_A[r,c] <= qtile[1,c]) {
          # Assign 1s
          ranked_A[r,c] <- 1
        } else if (ranked_A[r,c] > qtile[2,c]){
          # Assign 0s
          ranked_A[r,c] <- 0
        } else {
          # Assign score
          score <- (qtile[2,c] - ranked_A[r,c]) / q_range
          ranked_A[r,c] <- score
        }
      }
    }
  } else {
    # Discrete intervals
    for (c in 1:ncol(ranked_A)) {
      # Calculate qtile range
      q_range <- qtile[2,c] - qtile[1,c]
      for (r in 1:nrow(ranked_A)) {
        if (ranked_A[r,c] <= qtile[1,c]) {
          # Assign 1s
          ranked_A[r,c] <- 1
        } else if (ranked_A[r,c] > qtile[2,c]){
          # Assign 0s
          ranked_A[r,c] <- 0
        } else {
          # Assign score
          score <- (qtile[2,c] - ranked_A[r,c]) / q_range
          # Round off score
          ranked_A[r,c] <- bin(score, num_intervals)
        }
      }
    }
  }
  return (as.data.frame(ranked_A))
}


normaliseCDF <- function(df) {
  for (c in 1:ncol(df)) {
    notzero <- df[,c] != 0
    df[,c][notzero] <- rank(df[,c][notzero])
    df[,c] <- df[,c]/sum(notzero)
  }
  return(df)
}
