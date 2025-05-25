#' EM algorithm to estimate parameters of the NBS distribution
#'
#' Estimates the parameters \code{mu}, \code{sigma}, and \code{eta} of the
#' Normal-Birnbaum-Saunders (NBS) distribution from a sample \code{y} using the EM algorithm.
#'
#' @param y Numeric matrix of dimension \eqn{n \times p}, the observed sample.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{mu_est}}{Estimated mean vector \eqn{\mu}.}
#'   \item{\code{sigma_est}}{Estimated covariance matrix \eqn{\Sigma}.}
#'   \item{\code{eta_est}}{Estimated shape parameter \eqn{\eta}.}
#'   \item{\code{loglikelihood}}{Log-likelihood of the final estimates.}
#'   \item{\code{convergence}}{Logical, TRUE if converged before max iterations.}
#'   \item{\code{iterations}}{Number of iterations performed.}
#' }
#'
#' @details
#' The function uses the Expectation-Maximization (EM) algorithm to maximize the likelihood of the NBS model.
#' The Cholesky decomposition is used to compute matrix inverses and determinants efficiently.
#' Vectorized calculations improve computational speed for large datasets.
#'
#' Convergence is declared when both the change in log-likelihood and the change in parameters
#' are below specified tolerances or the maximum number of iterations is reached.
#'
#' @importFrom matrixcalc is.positive.definite
#' @importFrom stats cov
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' y <- matrix(rnorm(100 * 3), 100, 3)  # Example data
#' fit <- emNBS(y)
#' print(fit$mu_est)
#' }
emNBS <- function(y) {
  # Input validation
  if (!is.matrix(y)) stop("'y' must be a numeric matrix.")
  if (!isSymmetric(cov(y))) stop("Covariance matrix of 'y' must be symmetric.")
  if (!is.positive.definite(cov(y))) stop("Covariance matrix of 'y' must be positive definite.")

  n <- nrow(y)
  p <- ncol(y)
  a <- (p - 1) / 2

  # Initialization
  mu_est <- colMeans(y)
  eta_est <- 1
  sigma_est <- cov(y) / (1 + 1 / (2 * eta_est))

  # Initial log-likelihood
  loglik_est <- sum(apply(y, 1, function(yi) dNBS(yi, mu_est, sigma_est, eta_est, log_density = TRUE)))

  # EM control parameters
  iter <- 0
  maxiter <- 2000
  tol_loglik <- 1e-4
  tol_param <- 1e-4
  converged <- FALSE

  while (iter < maxiter && !converged) {
    mu <- mu_est
    sigma <- sigma_est
    eta <- eta_est

    # Cholesky decomposition for sigma
    chol_sigma <- tryCatch(chol(sigma), error = function(e) stop("Sigma is not positive definite during iteration."))
    sigma_inv <- chol2inv(chol_sigma)
    log_det_sigma <- 2 * sum(log(diag(chol_sigma)))

    # Compute delta for all observations: vectorized quadratic forms
    diffs <- sweep(y, 2, mu, FUN = "-")
    delta <- rowSums((diffs %*% sigma_inv) * diffs)

    w <- sqrt(eta * (eta + delta))

    # Log-density components for p_y calculation (vectorized)
    logd1 <- log(besselK(w, -a)) -
      (p / 2) * log(2 * pi) - 0.5 * log_det_sigma -
      log(besselK(eta, 0.5)) +
      a * log(eta / w)

    logd2 <- log(besselK(w, a - p)) -
      (p / 2) * log(2 * pi) - 0.5 * log_det_sigma -
      log(besselK(eta, -0.5)) -
      (a - p) * log(eta / w)

    py <- exp(logd1) / (exp(logd1) + exp(logd2))

    # Calculate v1 and v2 vectorized
    bessel_w_1_a <- besselK(w, 1 - a)
    bessel_w_m1_a <- besselK(w, -1 - a)
    bessel_w_a_p <- besselK(w, a - p)
    bessel_w_1_a_p <- besselK(w, 1 + a - p)
    bessel_w_m1_a_p <- besselK(w, -1 + a - p)

    term1_v1 <- py * (bessel_w_1_a / besselK(w, -a))
    term2_v1 <- (1 - py) * (bessel_w_1_a_p / bessel_w_a_p)
    v1 <- sqrt(1 + delta / eta) * (term1_v1 + term2_v1)

    term1_v2 <- py * (bessel_w_m1_a / besselK(w, -a))
    term2_v2 <- (1 - py) * (bessel_w_m1_a_p / bessel_w_a_p)
    v2 <- 1 / sqrt(1 + delta / eta) * (term1_v2 + term2_v2)

    # Update mu
    mu_est <- colSums(v2 * y) / sum(v2)

    # Update sigma
    diffs_new <- sweep(y, 2, mu_est, FUN = "-")
    sigma_est <- matrix(0, p, p)
    for (i in 1:n) {
      sigma_est <- sigma_est + v2[i] * tcrossprod(diffs_new[i, ])
    }
    sigma_est <- sigma_est / n
    sigma_est <- (sigma_est + t(sigma_est)) / 2

    # Update eta
    denom <- mean(v1) + mean(v2) - 2
    if (!is.finite(denom) || denom <= 0) {
      stop("Invalid value encountered in eta estimation (division by zero or negative).")
    }
    eta_est <- 1 / denom

    # Compute new log-likelihood
    loglik_new <- sum(apply(y, 1, function(yi) dNBS(yi, mu_est, sigma_est, eta_est, log_density = TRUE)))

    # Parameter change norm (Euclidean for mu and eta; Frobenius for sigma)
    param_change <- sqrt(
      sum((mu_est - mu)^2) +
        sum((sigma_est - sigma)^2) +
        (eta_est - eta)^2
    )
    # Check convergence
    if (abs(loglik_new - loglik_est) < tol_loglik && param_change < tol_param) {
      converged <- TRUE
    }

    loglik_est <- loglik_new
    iter <- iter + 1
  }

  list(
    mu_est = mu_est,
    sigma_est = sigma_est,
    eta_est = eta_est,
    loglikelihood = loglik_est,
    convergence = converged,
    iterations = iter
  )
}
