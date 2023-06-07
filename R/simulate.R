#' Simulate log-transformed gene expression microarray data
#'
#' @param m Number of features.
#' @param n Number of samples.
#' @param crosstab Matrix of contingency table specifying number of samples in
#'   each class and batch condition, with classes as rows and batches as columns.
#' @param delta Magnitude of additive batch effects.
#' @param gamma Magnitude of multiplicative batch effects. Variance parameter
#'   of Normal distribution modelling batch effects across samples in a batch.
#' @param epsilon Magnitude of random noise.
#' @param phi Percentage of differentially expressed features.
#' @param zeta Magnitude of class effects. Variance parameter of Normal
#'   distribution modelling log fold change.
#' @param a Shape parameter of Gamma distribution modelling basal expression.
#' @param b Scale parameter of Gamma distribution modelling basal expression.
#' @param c Inverse scale parameter of the sigmoid function
#' @param d Midpoint parameter of the sigmoid function
#' @param dropout Logical indicating whether to perform dropout
#' @return Matrix of dim (m, n).
#'
#' @importFrom extraDistr rlaplace
#' @export
simulate_microarray <- function(
  m, n, crosstab = NULL,
  delta = 1, gamma = 1,
  phi = 0.1, zeta = 1.5,
  a = 40, b = 0.2,
  c = 2, d = -6, dropout = TRUE,
  epsilon = 0.5,
  seed = 0
) {
  # checks
  stopifnot(n == sum(crosstab))
  # TODO: default is random 
  set.seed(seed)

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
  # TODO: provide option to specify batch effects magnitude of different batches
  log_beta <- matrix(extraDistr::rlaplace(m * n_batch, 0, delta), m, n_batch)
  
  alpha <- matrix(0, m, n)
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      k <- ks[j]
      alpha[i, j] <- rnorm(1, log_beta[i, k], gamma)
    }
  }
  X <- Z + alpha
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
    n_diffexpr = n_diffexpr, alpha = alpha,
    log_psi = log_psi, log_rho = log_rho, log_beta = log_beta
  )
}


#' Sigmoid function
#'
#' @param numeric scalar/vector/matrix
sigmoid <- function(x, c = 1, d = 0) 1 / (1 + exp(-(c * (x + d))))
