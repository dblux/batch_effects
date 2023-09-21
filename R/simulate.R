#' Simulate log-transformed gene expression microarray data
#'
#' @param m Number of features.
#' @param crosstab Matrix of contingency table specifying number of samples in
#'   each class-batch condition, with classes as rows and batches as columns.
#' @param delta Magnitude of additive batch effects.
#' @param gamma Magnitude of multiplicative batch effects. Variance parameter
#'   of Normal distribution modelling batch effects across samples in a batch.
#' @param phi Percentage of differentially expressed features.
#' @param zeta Magnitude of class effects. Variance parameter of Normal
#'   distribution modelling log fold change.
#' @param epsilon Magnitude of feature-wise variation.
#' @param kappa Magnitude of sample-specific variation (scaling factor).
#' @param a Shape parameter of Gamma distribution modelling basal expression.
#' @param b Scale parameter of Gamma distribution modelling basal expression.
#' @param dropout Logical indicating whether to perform dropout
#' @param c Inverse scale parameter of the sigmoid function
#' @param d Midpoint parameter of the sigmoid function
#' @param seed Numeric specifying random seed. Defaults to no seed.
#' @return Matrix of dim (m, n).
#'
#' @export
simulate_microarray <- function(
  m,
  crosstab,
  delta = 1,
  gamma = 1,
  phi = 0.1,
  zeta = 1.5,
  epsilon = 0.5, # limit = (, 1)
  kappa = 0.2, # limit = (, 0.3)
  a = 40,
  b = 0.2,
  dropout = FALSE,
  c = 2,
  d = -6,
  seed = NULL
) {
  if (!is.null(seed))
    set.seed(seed)

  n <- sum(crosstab)
  n_class <- nrow(crosstab)
  n_batch <- ncol(crosstab)
  gs <- rep(rep(seq_len(n_class), n_batch), crosstab) # class encoding
  ks <- rep(rep(seq_len(n_batch), each = n_class), crosstab) # batch encoding
  # metadata
  gs_alphabet <- LETTERS[gs]
  sid <- paste(paste0("ID", seq_len(n)), gs_alphabet, ks, sep = "_")
  metadata <- data.frame(
    class = gs_alphabet,
    batch = as.factor(ks),
    row.names = sid
  )

  log_psi <- rgamma(m, a, scale = b)
  log_rho <- cbind(
    rep(0, m), # class 1 has zero log fold change w.r.t. itself
    matrix(rnorm(m * (n_class - 1), 0, zeta), m, n_class - 1)
  )
  n_diffexpr <- round(phi * m, 0)
  # features not in the top phi percent according to log-fc are set to zero
  rank_rho <- apply(log_rho, 2, function(x) rank(-abs(x)))
  log_rho[rank_rho > n_diffexpr] <- 0

  Z <- matrix(0, m, n)
  colnames(Z) <- sid
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      g <- gs[j]
      Z[i, j] <- rnorm(1, log_psi[i] + log_rho[i, g], epsilon)
    }
  }
  log_beta <- matrix(rnorm(m * n_batch, 0, delta), m, n_batch)
  omega <- matrix(0, m, n) # batch effect terms
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      k <- ks[j]
      omega[i, j] <- rnorm(1, log_beta[i, k], gamma)
    }
  }
  log_alpha <- rnorm(n, 0, kappa) # log of sample specific scaling factor
  Z <- sweep(Z, 2, log_alpha, `+`)
  X <- Z + omega
  X[X < 0] <- 0 # set negative values to zero

  if (dropout) {
    P <- sigmoid(X, c, d)
    indicator <- matrix(0, m, n)
    for (i in seq_len(m)) {
      for (j in seq_len(n)) {
        indicator[i, j] <- rbinom(1, 1, P[i, j])
      }
    }
    X <- X * indicator
  }

  list(
    X = X, Z = Z, metadata = metadata,
    n_diffexpr = n_diffexpr, omega = omega,
    log_psi = log_psi, log_rho = log_rho, log_beta = log_beta
  )
}


#' Sigmoid function
#'
#' @param numeric scalar/vector/matrix
#' @param c Inverse scale parameter of the sigmoid function
#' @param d Midpoint parameter of the sigmoid function
sigmoid <- function(x, c = 1, d = 0) 1 / (1 + exp(-(c * (x + d))))