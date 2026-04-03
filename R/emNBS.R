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
emNBS <- function(y, eta_min = 1e-4, eta_max = 20, alpha_eta = 0.25, ridge = 1e-8) {
  if (!is.matrix(y)) stop("'y' must be a numeric matrix.")
  Sy <- cov(y)
  if (!isSymmetric(Sy)) stop("Covariance matrix of 'y' must be symmetric.")
  if (!matrixcalc::is.positive.definite(Sy)) stop("Covariance matrix of 'y' must be positive definite.")

  clamp_eta <- function(x) pmax(eta_min, pmin(x, eta_max))

  n <- nrow(y)
  p <- ncol(y)
  a <- (p - 1) / 2

  mu_est <- colMeans(y)
  eta_est <- 1
  sigma_est <- Sy / (1 + 1 / (2 * eta_est))
  sigma_est <- (sigma_est + t(sigma_est)) / 2 + ridge * diag(p)

  loglik_est <- sum(apply(
    y, 1,
    function(yi) dNBS(yi, mu_est, sigma_est, clamp_eta(eta_est), log_density = TRUE)
  ))

  iter <- 0
  maxiter <- 2000
  tol_loglik <- 1e-4
  tol_param <- 1e-4
  converged <- FALSE

  while (iter < maxiter && !converged) {
    mu <- mu_est
    sigma <- sigma_est
    eta <- clamp_eta(eta_est)

    chol_sigma <- tryCatch(
      chol(sigma),
      error = function(e) stop("Sigma is not positive definite during iteration.")
    )
    sigma_inv <- chol2inv(chol_sigma)
    log_det_sigma <- 2 * sum(log(diag(chol_sigma)))

    diffs <- sweep(y, 2, mu, FUN = "-")
    delta <- rowSums((diffs %*% sigma_inv) * diffs)

    w <- sqrt(eta * (eta + delta))
    w <- pmax(w, 1e-12)

    K_w_ma   <- besselK(w, -a, expon.scaled = TRUE)
    K_w_ap   <- besselK(w, a - p, expon.scaled = TRUE)
    K_eta_p  <- besselK(eta, 0.5, expon.scaled = TRUE)
    K_eta_m  <- besselK(eta, -0.5, expon.scaled = TRUE)

    logd1 <- log(K_w_ma) -
      (p / 2) * log(2 * pi) - 0.5 * log_det_sigma -
      log(K_eta_p) +
      a * log(eta / w)

    logd2 <- log(K_w_ap) -
      (p / 2) * log(2 * pi) - 0.5 * log_det_sigma -
      log(K_eta_m) -
      (a - p) * log(eta / w)

    log_ratio <- pmax(pmin(logd2 - logd1, 700), -700)
    py <- 1 / (1 + exp(log_ratio))

    K_w_1ma   <- besselK(w, 1 - a, expon.scaled = TRUE)
    K_w_m1ma  <- besselK(w, -1 - a, expon.scaled = TRUE)
    K_w_1ap   <- besselK(w, 1 + a - p, expon.scaled = TRUE)
    K_w_m1ap  <- besselK(w, -1 + a - p, expon.scaled = TRUE)

    term1_v1 <- py * (K_w_1ma / K_w_ma)
    term2_v1 <- (1 - py) * (K_w_1ap / K_w_ap)
    v1 <- sqrt(1 + delta / eta) * (term1_v1 + term2_v1)

    term1_v2 <- py * (K_w_m1ma / K_w_ma)
    term2_v2 <- (1 - py) * (K_w_m1ap / K_w_ap)
    v2 <- (term1_v2 + term2_v2) / sqrt(1 + delta / eta)

    if (any(!is.finite(v1)) || any(!is.finite(v2))) {
      stop("Non-finite values in E-step.")
    }

    mu_est <- colSums(v2 * y) / sum(v2)

    diffs_new <- sweep(y, 2, mu_est, FUN = "-")
    sigma_est <- matrix(0, p, p)
    for (i in 1:n) {
      sigma_est <- sigma_est + v2[i] * tcrossprod(diffs_new[i, ])
    }
    sigma_est <- sigma_est / n
    sigma_est <- (sigma_est + t(sigma_est)) / 2 + ridge * diag(p)

    denom <- mean(v1) + mean(v2) - 2
    if (!is.finite(denom) || denom <= 0) {
      stop("Invalid value encountered in eta estimation.")
    }

    eta_raw <- 1 / denom
    eta_raw <- clamp_eta(eta_raw)

    # actualización amortiguada
    eta_est <- clamp_eta((1 - alpha_eta) * eta + alpha_eta * eta_raw)

    loglik_new <- sum(apply(
      y, 1,
      function(yi) dNBS(yi, mu_est, sigma_est, eta_est, log_density = TRUE)
    ))

    param_change <- sqrt(
      sum((mu_est - mu)^2) +
        sum((sigma_est - sigma)^2) +
        (eta_est - eta)^2
    )

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
