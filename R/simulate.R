#' Simulate log-transformed gene expression microarray data
#'
#' @param m Number of features.
#' @param crosstab Matrix of contingency table specifying number of samples in
#'   each class-batch condition, with classes as rows and batches as columns.
#' @param delta Magnitude of additive batch effects.
#' @param gamma Magnitude of multiplicative batch effects. Variance parameter
#'   of Normal distribution modelling batch effects across samples in a batch.
#' @param phi Percentage of differentially expressed features.
#' @param c Shape parameter of Gamma distribution modelling log fold-change.
#' @param d Rate parameter of Gamma distribution modelling log fold-change.
#' @param epsilon Magnitude of feature-wise variation.
#' @param kappa Magnitude of sample-specific variation (scaling factor).
#' @param a Shape parameter of Gamma distribution modelling basal expression.
#' @param b Rate parameter of Gamma distribution modelling basal expression.
#' @param dropout Logical indicating whether to perform dropout
#' @param r Inverse scale parameter of the sigmoid function.
#' @param s Midpoint parameter of the sigmoid function.
#' @param seed Numeric specifying random seed. Defaults to no seed.
#' @return Matrix of dim (m, n).
#'
#' @export
simulate_microarray <- function(
  m,
  crosstab,
  delta = 1,
  gamma = 0.5,
  phi = 0.2,
  c = 10,
  d = 6,
  epsilon = 0.5, # limit = (, 1)
  kappa = 0.2, # limit = (, 0.3)
  a = 40,
  b = 5,
  dropout = FALSE,
  r = 2,
  s = -6,
  seed = NA
) {
  # Record parameters
  params <- c(
    crosstab = crosstab,
    delta = delta, gamma = gamma,
    phi = phi, c = c, d = d,
    epsilon = epsilon, kappa = kappa,
    a = a, b = b, 
    dropout = dropout, r = r, s = s,
    seed = seed
  )

  if (!is.na(seed))
    set.seed(seed)

  n <- sum(crosstab)
  n_class <- nrow(crosstab)
  n_batch <- ncol(crosstab)
  gs <- rep(rep(seq_len(n_class), n_batch), crosstab) # class encoding
  ks <- rep(rep(seq_len(n_batch), each = n_class), crosstab) # batch encoding

  # Metadata
  gs_alphabet <- LETTERS[gs]
  sid <- paste(paste0("ID", seq_len(n)), gs_alphabet, ks, sep = "_")
  metadata <- data.frame(
    class = gs_alphabet,
    batch = as.factor(ks),
    row.names = sid
  )

  log_psi <- rgamma(m, a, rate = b)
  # Log fold-change factors for each class
  # Class A has zero log fold change w.r.t. itself
  log_rho <- matrix(0, m, n_class)
  colnames(log_rho) <- LETTERS[seq_len(n_class)]
  diff.features <- NULL
  if (n_class > 1) {
    n_upreg <- n_downreg <- round(phi * m / 2, 0)
    n_diffexpr <- n_upreg + n_downreg
    for (g in seq(2, n_class)) {
      diff.features <- sort(sample(seq_len(m), n_diffexpr))
      upreg.features <- sort(sample(diff.features, n_upreg))
      downreg.features <- setdiff(diff.features, upreg.features)
      for (i in upreg.features) {
        log_rho[i, g] <- rgamma(1, c, rate = d)
      }
      for (i in downreg.features) {
        log_rho[i, g] <- -rgamma(1, c, rate = d)
      }
    }
  }

  # Base expression values with class effects
  Z <- matrix(0, m, n)
  colnames(Z) <- sid
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      g <- gs[j]
      Z[i, j] <- rnorm(1, log_psi[i] + log_rho[i, g], epsilon)
    }
  }

  # Sample specific scaling term (in log space)
  log_alpha <- rnorm(n, 0, kappa)
  W <- sweep(Z, 2, log_alpha, `+`)

  # Batch effects
  log_beta <- matrix(rnorm(m * n_batch, 0, delta), m, n_batch)
  omega <- matrix(0, m, n)
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      k <- ks[j]
      omega[i, j] <- rnorm(1, log_beta[i, k], gamma)
    }
  }

  X <- W + omega
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
    X = X, metadata = metadata,
    diff.features = diff.features,
    W = W, batch.terms = omega,
    class.logfc = log_rho, batch.logfc = log_beta,
    log_psi = log_psi, Z = Z, log_alpha = log_alpha,
    params = params
  )
}


#' Sigmoid function
#'
#' @param numeric scalar/vector/matrix
#' @param r Inverse scale parameter of the sigmoid function
#' @param s Midpoint parameter of the sigmoid function
sigmoid <- function(x, r = 1, s = 0) 1 / (1 + exp(-(r * (x + s))))
