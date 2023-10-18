 # Log2 transforms data and handles -Inf values
log2_transform <- function(df) {
  log2_df <- log2(df)
  logical_df <- is.infinite(data.matrix(log2_df))
  log2_df[logical_df] <- 0
  return(log2_df)
}


#' Converts radians to degrees
rad2degree <- function(rad) rad/pi * 180


#' Calculates angle between two vectors (in degrees)
calcAngleVectors <- function(x, y) {
  # Prevents error when number is incorrectly calculated..
  # to be slightly above 1
  cosine_sim <- calcCosineSim(x, y)
  if (cosine_sim > 1) {
    print(sprintf("Cosine similarity: %.20f -> Rounded off!", cosine_sim))
    cosine_sim <- 1
  }
  return(rad2degree(acos(cosine_sim)))
}


#' @param rad rotation angle (counter-clockwise) in radians
#' @return rotation matrix for 2D vector
calcRotationMatrix <- function(rad) {
  matrix(
    c(cos(rad), -sin(rad), sin(rad), cos(rad)),
    2, 2, byrow = T
  )
}

calcTPM <- function(X, feature_length) {
  stopifnot(nrow(X) == length(feature_length))
  # Divide by feature length
  X_length <- X/feature_length
  tpm <- sweep(X_length, 2, colSums(X_length), "/") * 1e+6
  return(tpm)
}


norm <- function(x, p) sum(x^p) ^ (1/p)


#' Calculates l2-norm of vector
l2norm <- function(x) sqrt(sum(x ^ 2))


#' Calculates cosine similarity between two vectors
#' @return Scalar cos(theta)
calcCosineSim <- function(x, y) {
  stopifnot(length(x) == length(y))
  
  if (all(x == 0) | all(y == 0))
    stop("Unable to calculate cosine similarity on zero vector..")
  
  return(sum(x*y)/sqrt(sum(x^2)*sum(y^2)))
}

# Function that calculates the matrix of mean differences from the patient and control matrix
# Arguments: a <- larger matrix, b <- smaller matrix
# CLT: As n -> infinity, approximates more accurately a normal distribution.
# If underlying distribution is normal, then no need for n to be large
# Returns: mean_diff matrix with dimension of smaller matrix
old_calc_diff <- function(a, b) {
  a_colsum <- apply(a, 1, sum)
  a_size <- ncol(a)
  mean_diff <- array(numeric(), dim(b))
  rownames(mean_diff) <- rownames(b)
  
  for (r in 1:nrow(b)) {
    for (c in 1:ncol(b)) {
      mean_diff[r,c] <- (a_colsum[r]/a_size) - b[r,c]
    }
  }
  return (mean_diff)
}

# Function that calculates the matrix of mean differences from the patient and control matrix
# Arguments: a <- larger matrix, b <- smaller matrix
# CLT: As n -> infinity, approximates more accurately a normal distribution.
# If underlying distribution is normal, then no need for n to be large
# Returns: mean_diff matrix with dimension of smaller matrix
# pct= 0.75 smallest=4. pct=0.8 smallest=5
calc_diff <- function(a, b, sample_pct) {
  a_size <- ncol(a)
  sample_size <- floor(sample_pct * a_size)
  if (sample_size == a_size) {print("Sample size too small..")}
  
  mean_diff <- array(numeric(), dim(b))
  rownames(mean_diff) <- rownames(b)
  
  for (r in 1:nrow(b)) {
    for (c in 1:ncol(b)) {
      subsample <- unlist(sample(a[r,], sample_size, replace = F))
      mean_diff[r,c] <- mean(subsample) - b[r,c]
    }
  }
  return (mean_diff)
}

all_diff <- function(a, b) {
  all_diff <- array(numeric(), c(nrow(a), ncol(a)*ncol(b)))
  rownames(all_diff) <- rownames(a)
  for (r in 1:nrow(a)) {
    print(r)
    for (col_a in 1:ncol(a)) {
      for (col_b in 1:ncol(b)) {
        c <- col_b + (col_a - 1) * ncol(b)
        all_diff[r,c] <- b[r,col_b] - a[r,col_a]
      }
    }
  }
  return (all_diff)
}

# Calculate one-sample t-test statistic
# Using sample sd as estimator of population sd
# Use t-distribution to calculate p value
ttest_onesample <- function(vector, mu) {
  n <- length(vector)
  t_stat <- (mean(vector) - mu) / (sd(vector)/sqrt(n))
  if (is.infinite(t_stat)|is.na(t_stat)) {
    return (NaN)
  }
  if (t_stat < 0) {
    p_value <- pt(t_stat, n-1) * 2
  } else if (t_stat > 0) {
    p_value <- (1-pt(t_stat, n-1)) * 2
  } else {
    # t_stat == 0
    p_value <- 1
  }
  return (p_value)
}

# Naive row-wise two-sample t-test for every probe
# Does a t-test between every row of matrices a and b
# Arguments: Dataframe with both clases, size of class A, ...
calc_ttest <- function(df, size_a, flag = "pvalue", is_paired = F) {
  row_pvalue <- function(row_vec) {
    return(t.test(row_vec[1:size_a],
                  row_vec[-(1:size_a)],
                  paired = is_paired)$p.value)
  }
  row_tstat <- function(row_vec) {
    return(t.test(row_vec[1:size_a],
                  row_vec[-(1:size_a)],
                  paired = is_paired)$statistic)
  }
    
  row_both <- function(row_vec) {
    ttest_obj <- t.test(row_vec[1:size_a],
                        row_vec[-(1:size_a)],
                        paired = is_paired)
    return(c(pvalue = ttest_obj$p.value,
             ttest_obj$statistic))
  }
  
  if (flag == "pvalue") {
    ttest_vec <- apply(df, 1, function(row) tryCatch(row_pvalue(row),
                                                     error = function(e) NA))
  } else if (flag == "tstat") {
    ttest_vec <- apply(df, 1, function(row) tryCatch(row_tstat(row),
                                                     error = function(e) NA))
  } else if (flag == "both") {
    ttest_vec <- apply(df, 1, function(row) tryCatch(row_both(row),
                                                     error = function(e) NA))
  } else {
    stop("Flag not in options pvalue or tstat..")
  }
  return(ttest_vec)
}

# Arguments: 2 dataframes that are not log-transformed
# Log-fold change (class1/class2)
calc_logfc <- function(df1, df2, func = mean, logged = T) {
  if (logged) {
    ## Mean of logged values: Geometric mean
    return(apply(df1, 1, func) - apply(df2, 1, func))
  } else {
    # Minimum value of both df besides 0 chosen as prior value
    prior_value <- min(c(df1[df1 != 0], df2[df2 != 0]))
    print(paste("Prior value:", prior_value))
    vec1 <- apply(df1, 1, func)
    vec2 <- apply(df2, 1, func)
    # log2(0) = -Inf; log2(Inf) = -Inf; log2(0/0) = NaN
    # Reassigns 0s with prior_expr
    ## Replace means that are zero with prior value
    vec1[vec1 == 0] <- prior_value
    vec2[vec2 == 0] <- prior_value
    fc <- vec1/vec2
    return(log2(fc))
  }
}


#' Naive row-wise two-sample t-test for every probe
#' Does a t-test between every row of matrices a and b
#' Arguments: Dataframe with both clases, size of class A, ...
#' @param FUN Hypothesis test function. E.g. t.test
#' @param X Data frame with dim (n_samples, n_features)
#' @param Y Data frame with dim (n_samples, n_features)
#' @param paired logical indicating whether X contains paired data or not
calc_univariate <- function(
  FUN, X, Y = NULL, n_split = NULL,
  paired = FALSE, flag = "p.value"
) {
  if (paired) {
    # Y and n can be both NULL    
    n_split <- ncol(X) / 2
  }
  
  if (!is.null(Y)){
    if (!is.data.frame(X) | !is.data.frame(Y))
      stop("Type error: Both X and Y have to be of class data.frame.")
    X1 <- data.frame(t(X))
    Y1 <- data.frame(t(Y))
    list_htest <- mapply(FUN, X1, Y1, paired = paired, SIMPLIFY = F)
  } else if(!is.null(n_split)) {
    list_htest <- apply(X, 1, function(x) {
      idx <- seq(1, n_split)
      FUN(x[idx], x[-idx], paired = paired)
    })
  } else {
    stop("If paired is FALSE either Y or n_split must be provided.")   
  }
  
  if (flag == "p.value") {
    return(sapply(list_htest, function(obj) obj$p.value))
  } else if (flag == "statistic") {
    return(sapply(list_htest, function(obj) obj$statistic))
  } else {
    return(list_htest)
  }
}
