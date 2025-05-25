#' EM Algorithm for Multivariate Normal-Birnbaum-Saunders with Equicorrelation Structure
#'
#' Estimates the parameters of the multivariate Normal-Birnbaum-Saunders (NBS) distribution
#' under the null hypothesis that the covariance matrix has an equicorrelation structure,
#' i.e., \eqn{\Sigma = \sigma^2 R} with \eqn{R = (1 - \rho)I + \rho J}.
#'
#' This function uses the Expectation-Maximization (EM) algorithm to obtain maximum likelihood
#' estimates of the mean vector \eqn{\mu}, the scale parameter \eqn{\eta}, the common variance
#' \eqn{\sigma^2}, and the equicorrelation parameter \eqn{\rho}.
#'
#' @param y A numeric matrix of size \eqn{n \times p}, where each row is an observation from
#' the multivariate NBS distribution.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{mu_est}}{Estimated mean vector \eqn{\mu}.}
#'   \item{\code{sigma_est}}{Estimated covariance matrix \eqn{\Sigma = \sigma^2 R}.}
#'   \item{\code{eta_est}}{Estimated shape parameter \eqn{\eta}.}
#'   \item{\code{rho_est}}{Estimated equicorrelation parameter \eqn{\rho}.}
#'   \item{\code{sigma2_est}}{Estimated common variance component \eqn{\sigma^2}.}
#'   \item{\code{loglikelihood}}{Log-likelihood value at convergence.}
#'   \item{\code{convergence}}{Logical indicating whether the algorithm converged.}
#'   \item{\code{iterations}}{Number of iterations until convergence.}
#' }
#'
#' @details
#' This function assumes the covariance structure is of the form \eqn{\Sigma = \sigma^2 R}, where
#' \eqn{R = (1 - \rho)I + \rho J} is the equicorrelation matrix. This is useful for hypothesis testing
#' involving equicorrelation versus a general covariance structure. The function implements safeguards
#' to ensure numerical stability and positive-definiteness of the estimated matrices.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' p <- 3
#' mu <- rep(0, p)
#' eta <- 1
#' sigma2 <- 2
#' rho <- 0.3
#' R <- (1 - rho) * diag(p) + rho * matrix(1, p, p)
#' Sigma <- sigma2 * R
#' y <- rNBS(n, mu, Sigma, eta)
#' fit <- emNBS_equicorrelation(y)
#' }
#'
#' @importFrom matrixcalc is.positive.definite
#'
#' @seealso \code{\link{emNBS}} for unrestricted NBS estimation.
#'
#' @export
emNBS_equicorrelation <- function(y) {
  # Initial validation
  if (!is.matrix(y)) stop("'y' must be a numeric matrix.")
  if (!isSymmetric(cov(y))) stop("The covariance matrix of 'y' is not symmetric.")
  if (!is.positive.definite(cov(y))) stop("The covariance matrix of 'y' is not positive definite.")

  n <- nrow(y)
  p <- ncol(y)
  a <- (p - 1) / 2
  I <- diag(p)
  J <- matrix(1, p, p)

  # Initialization
  mu_est <- colMeans(y)
  eta_est <- 1
  S <- cov(y)
  sigma2_est <- mean(diag(S))
  rho_est <- mean(S[lower.tri(S)]) / sigma2_est
  R_est <- (1 - rho_est) * I + rho_est * J
  sigma_est <- sigma2_est * R_est

  loglik_est <- sum(apply(y, 1, function(yi) dNBS(yi, mu_est, sigma_est, eta_est, log_density = TRUE)))

  # EM control settings
  iter <- 0
  maxiter <- 2000
  tol_loglik <- 1e-4
  tol_param <- 1e-4
  converged <- FALSE

  while (iter < maxiter && !converged) {
    mu <- mu_est
    sigma <- sigma_est
    eta <- eta_est

    chol_sigma <- tryCatch(chol(sigma), error = function(e) stop("Sigma is not positive definite during iteration."))
    sigma_inv <- chol2inv(chol_sigma)
    log_det_sigma <- 2 * sum(log(diag(chol_sigma)))

    # Vectorization
    diffs <- sweep(y, 2, mu, FUN = "-")
    delta <- rowSums((diffs %*% sigma_inv) * diffs)
    w <- sqrt(eta * (eta + delta))

    logd1 <- log(besselK(w, -a)) - (p / 2) * log(2 * pi) - 0.5 * log_det_sigma -
      log(besselK(eta, 0.5)) + a * log(eta / w)

    logd2 <- log(besselK(w, a - p)) - (p / 2) * log(2 * pi) - 0.5 * log_det_sigma -
      log(besselK(eta, -0.5)) - (a - p) * log(eta / w)

    # Numerical stability
    max_log <- pmax(logd1, logd2)
    denom <- exp(logd1 - max_log) + exp(logd2 - max_log)
    py <- exp(logd1 - max_log) / denom

    # Components for v1 and v2
    b1 <- besselK(w, 1 - a)
    b2 <- besselK(w, -1 - a)
    b3 <- besselK(w, a - p)
    b4 <- besselK(w, 1 + a - p)
    b5 <- besselK(w, -1 + a - p)

    v1 <- sqrt(1 + delta / eta) * (py * b1 / besselK(w, -a) + (1 - py) * b4 / b3)
    v2 <- 1 / sqrt(1 + delta / eta) * (py * b2 / besselK(w, -a) + (1 - py) * b5 / b3)

    # Update mu
    mu_est <- colSums(v2 * y) / sum(v2)

    # Update sigma2 and rho
    diffs <- sweep(y, 2, mu_est)
    norm2 <- rowSums(diffs^2)
    sum_vec <- rowSums(diffs)

    cumsum1 <- sum(v2 * norm2)
    cumsum2 <- sum(v2 * sum_vec^2)

    rho_est <- (cumsum2 - cumsum1) / ((p - 1) * cumsum1)
    rho_est <- max(min(rho_est, 0.99), -1 / (p - 1) + 1e-6)  # avoid out-of-range values

    denom_rho <- (1 - rho_est) * (1 - rho_est + rho_est * p)
    if (denom_rho <= 0) stop("Invalid value in sigma2 estimation.")

    sigma2_est <- ((1 - rho_est + rho_est * p) * cumsum1 - rho_est * cumsum2) / (n * p * denom_rho)
    if (sigma2_est <= 0) stop("Sigma2 must be positive.")

    R_est <- (1 - rho_est) * I + rho_est * J
    sigma_est <- sigma2_est * R_est

    # Update eta
    denom_eta <- mean(v1) + mean(v2) - 2
    if (!is.finite(denom_eta) || denom_eta <= 0) stop("Error in eta estimation (invalid division).")
    eta_est <- 1 / denom_eta

    # Log-likelihood and convergence
    loglik_new <- sum(apply(y, 1, function(yi) dNBS(yi, mu_est, sigma_est, eta_est, log_density = TRUE)))

    param_change <- sqrt(sum((mu_est - mu)^2) + (eta_est - eta)^2 + (sigma2_est - mean(diag(sigma)))^2)
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
    rho_est = rho_est,
    sigma2_est = sigma2_est,
    loglikelihood = loglik_est,
    convergence = converged,
    iterations = iter
  )
}
