#' EM Algorithm for Multivariate Normal-Birnbaum-Saunders
#' with Equicorrelation Structure
#'
#' Estimates the parameters of the multivariate Normal-Birnbaum-Saunders (NBS)
#' distribution under the null hypothesis that the covariance matrix has an
#' equicorrelation structure, i.e., Sigma = sigma^2 R with
#' R = (1 - rho) I + rho J.
#'
#' @param y A numeric matrix of size n x p.
#' @param eta_min Lower bound for eta.
#' @param eta_max Upper bound for eta.
#' @param alpha_eta Damping factor for eta update in (0,1].
#' @param ridge Small positive ridge added to Sigma for numerical stability.
#' @param maxiter Maximum number of EM iterations.
#' @param tol_loglik Tolerance for log-likelihood convergence.
#' @param tol_param Tolerance for parameter convergence.
#'
#' @return A list with components:
#' \describe{
#'   \item{mu_est}{Estimated mean vector.}
#'   \item{sigma_est}{Estimated covariance matrix Sigma = sigma^2 R.}
#'   \item{eta_est}{Estimated shape parameter eta.}
#'   \item{rho_est}{Estimated equicorrelation parameter rho.}
#'   \item{sigma2_est}{Estimated common variance component sigma^2.}
#'   \item{loglikelihood}{Final log-likelihood value.}
#'   \item{convergence}{Logical; TRUE if converged.}
#'   \item{iterations}{Number of EM iterations performed.}
#' }
#'
#' @importFrom matrixcalc is.positive.definite
#' @export
emNBS_equicorrelation <- function(y,
                                  eta_min = 1e-4,
                                  eta_max = 20,
                                  alpha_eta = 0.25,
                                  ridge = 1e-8,
                                  maxiter = 2000,
                                  tol_loglik = 1e-4,
                                  tol_param = 1e-4) {
  # -----------------------------
  # 1. Input validation
  # -----------------------------
  if (!is.matrix(y)) stop("'y' must be a numeric matrix.")

  Sy <- stats::cov(y)
  if (!isSymmetric(Sy)) {
    stop("The covariance matrix of 'y' is not symmetric.")
  }
  if (!matrixcalc::is.positive.definite(Sy)) {
    stop("The covariance matrix of 'y' is not positive definite.")
  }

  if (!is.numeric(eta_min) || !is.numeric(eta_max) || eta_min <= 0 || eta_min >= eta_max) {
    stop("Require 0 < eta_min < eta_max.")
  }
  if (!is.numeric(alpha_eta) || alpha_eta <= 0 || alpha_eta > 1) {
    stop("'alpha_eta' must be in (0,1].")
  }
  if (!is.numeric(ridge) || ridge < 0) {
    stop("'ridge' must be nonnegative.")
  }

  # -----------------------------
  # 2. Helpers
  # -----------------------------
  clamp_eta <- function(x) pmax(eta_min, pmin(x, eta_max))

  clamp_rho <- function(rho, p) {
    # rho in (-1/(p-1), 1)
    lower <- -1 / (p - 1) + 1e-6
    upper <- 0.99
    pmax(lower, pmin(rho, upper))
  }

  safe_log_besselK <- function(x, nu) {
    val <- besselK(x, nu, expon.scaled = TRUE)
    if (any(!is.finite(val)) || any(val <= 0)) {
      stop("Non-finite or non-positive BesselK value encountered.")
    }
    log(val)
  }

  # -----------------------------
  # 3. Constants and initialization
  # -----------------------------
  n <- nrow(y)
  p <- ncol(y)
  a <- (p - 1) / 2
  I <- diag(p)
  J <- matrix(1, p, p)

  mu_est <- colMeans(y)
  eta_est <- 1

  sigma2_est <- mean(diag(Sy))
  rho_est <- mean(Sy[lower.tri(Sy)]) / sigma2_est
  rho_est <- clamp_rho(rho_est, p)

  R_est <- (1 - rho_est) * I + rho_est * J
  sigma_est <- sigma2_est * R_est
  sigma_est <- (sigma_est + t(sigma_est)) / 2 + ridge * I

  loglik_est <- sum(apply(
    y, 1,
    function(yi) dNBS(yi, mu_est, sigma_est, clamp_eta(eta_est), log_density = TRUE)
  ))

  iter <- 0
  converged <- FALSE

  # -----------------------------
  # 4. EM loop
  # -----------------------------
  while (iter < maxiter && !converged) {
    mu <- mu_est
    sigma <- sigma_est
    eta <- clamp_eta(eta_est)
    rho <- rho_est
    sigma2 <- sigma2_est

    chol_sigma <- tryCatch(
      chol(sigma),
      error = function(e) stop("Sigma is not positive definite during iteration.")
    )
    sigma_inv <- chol2inv(chol_sigma)
    log_det_sigma <- 2 * sum(log(diag(chol_sigma)))

    # Mahalanobis distances
    diffs <- sweep(y, 2, mu, FUN = "-")
    delta <- rowSums((diffs %*% sigma_inv) * diffs)

    # Stable w
    w <- sqrt(eta * (eta + delta))
    w <- pmax(w, 1e-12)

    # -----------------------------
    # E-step
    # -----------------------------
    # Use scaled BesselK to improve stability
    logK_w_ma   <- safe_log_besselK(w, -a)
    logK_w_ap   <- safe_log_besselK(w, a - p)
    logK_eta_p  <- safe_log_besselK(eta,  0.5)
    logK_eta_m  <- safe_log_besselK(eta, -0.5)

    logd1 <- logK_w_ma -
      (p / 2) * log(2 * pi) -
      0.5 * log_det_sigma -
      logK_eta_p +
      a * log(eta / w)

    logd2 <- logK_w_ap -
      (p / 2) * log(2 * pi) -
      0.5 * log_det_sigma -
      logK_eta_m -
      (a - p) * log(eta / w)

    # Stable py = 1 / (1 + exp(logd2 - logd1))
    log_ratio <- logd2 - logd1
    log_ratio <- pmax(pmin(log_ratio, 700), -700)
    py <- 1 / (1 + exp(log_ratio))

    # Bessel terms for v1 and v2
    K_w_ma   <- besselK(w, -a, expon.scaled = TRUE)
    K_w_ap   <- besselK(w, a - p, expon.scaled = TRUE)
    K_w_1ma  <- besselK(w, 1 - a, expon.scaled = TRUE)
    K_w_m1ma <- besselK(w, -1 - a, expon.scaled = TRUE)
    K_w_1ap  <- besselK(w, 1 + a - p, expon.scaled = TRUE)
    K_w_m1ap <- besselK(w, -1 + a - p, expon.scaled = TRUE)

    if (any(!is.finite(K_w_ma))   || any(!is.finite(K_w_ap)) ||
        any(!is.finite(K_w_1ma))  || any(!is.finite(K_w_m1ma)) ||
        any(!is.finite(K_w_1ap))  || any(!is.finite(K_w_m1ap))) {
      stop("Non-finite BesselK values encountered in E-step.")
    }

    term1_v1 <- py * (K_w_1ma / K_w_ma)
    term2_v1 <- (1 - py) * (K_w_1ap / K_w_ap)
    v1 <- sqrt(1 + delta / eta) * (term1_v1 + term2_v1)

    term1_v2 <- py * (K_w_m1ma / K_w_ma)
    term2_v2 <- (1 - py) * (K_w_m1ap / K_w_ap)
    v2 <- (term1_v2 + term2_v2) / sqrt(1 + delta / eta)

    if (any(!is.finite(v1)) || any(!is.finite(v2))) {
      stop("Non-finite values encountered in v1 or v2.")
    }
    if (sum(v2) <= 0 || !is.finite(sum(v2))) {
      stop("Invalid denominator in mu update.")
    }

    # -----------------------------
    # M-step: update mu
    # -----------------------------
    mu_est <- colSums(v2 * y) / sum(v2)

    # -----------------------------
    # M-step: update rho and sigma2
    # -----------------------------
    diffs_new <- sweep(y, 2, mu_est, FUN = "-")
    norm2 <- rowSums(diffs_new^2)
    sum_vec <- rowSums(diffs_new)

    cumsum1 <- sum(v2 * norm2)
    cumsum2 <- sum(v2 * sum_vec^2)

    rho_est <- (cumsum2 - cumsum1) / ((p - 1) * cumsum1)
    rho_est <- clamp_rho(rho_est, p)

    denom_rho <- (1 - rho_est) * (1 - rho_est + rho_est * p)
    if (!is.finite(denom_rho) || denom_rho <= 0) {
      stop("Invalid value in sigma2 estimation.")
    }

    sigma2_est <- ((1 - rho_est + rho_est * p) * cumsum1 - rho_est * cumsum2) /
      (n * p * denom_rho)

    if (!is.finite(sigma2_est) || sigma2_est <= 0) {
      stop("Sigma2 must be positive and finite.")
    }

    R_est <- (1 - rho_est) * I + rho_est * J
    sigma_est <- sigma2_est * R_est
    sigma_est <- (sigma_est + t(sigma_est)) / 2 + ridge * I

    # -----------------------------
    # M-step: update eta
    # -----------------------------
    denom_eta <- mean(v1) + mean(v2) - 2
    if (!is.finite(denom_eta) || denom_eta <= 0) {
      stop("Error in eta estimation (invalid division).")
    }

    eta_raw <- 1 / denom_eta
    eta_raw <- clamp_eta(eta_raw)

    # Damped update
    eta_est <- (1 - alpha_eta) * eta + alpha_eta * eta_raw
    eta_est <- clamp_eta(eta_est)

    # -----------------------------
    # Log-likelihood
    # -----------------------------
    loglik_new <- sum(apply(
      y, 1,
      function(yi) dNBS(yi, mu_est, sigma_est, eta_est, log_density = TRUE)
    ))

    if (!is.finite(loglik_new)) {
      stop("Non-finite log-likelihood encountered.")
    }

    param_change <- sqrt(
      sum((mu_est - mu)^2) +
        (eta_est - eta)^2 +
        (rho_est - rho)^2 +
        (sigma2_est - sigma2)^2
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
    rho_est = rho_est,
    sigma2_est = sigma2_est,
    loglikelihood = loglik_est,
    convergence = converged,
    iterations = iter
  )
}
