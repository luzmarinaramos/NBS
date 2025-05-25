#' Density function of the Multivariate Normal Birnbaum-Saunders (NBS) distribution
#'
#' Computes the density (or log-density) of the Multivariate Normal Birnbaum-Saunders distribution
#' at a given point \code{y}.
#'
#' @param y Numeric vector. Point at which to evaluate the density.
#' @param mu Numeric vector. Mean vector of the distribution.
#' @param sigma Numeric matrix. Covariance matrix (must be symmetric positive definite).
#' @param eta Positive numeric scalar. Shape parameter of the Birnbaum-Saunders mixture.
#' @param log_density Logical; if TRUE, returns the log-density. Default is FALSE.
#'
#' @return Numeric scalar. Density (or log-density) evaluated at \code{y}.
#'
#' @examples
#' mu <- c(0, 0)
#' sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' eta <- 1
#' y <- c(1, 1)
#' dNBS(y, mu, sigma, eta)
#'
#' @export
dNBS <- function(y, mu, sigma, eta, log_density = FALSE) {
  # Check inputs
  if (!is.numeric(y)) stop("y must be numeric vector")
  if (!is.numeric(mu)) stop("mu must be numeric vector")
  if (!is.matrix(sigma)) stop("sigma must be a numeric matrix")
  if (nrow(sigma) != ncol(sigma)) stop("sigma must be a square matrix")
  if (length(mu) != length(y)) stop("mu and y must have the same length")
  if (length(mu) != nrow(sigma)) stop("Dimensions of mu and sigma do not match")
  if (!is.numeric(eta) || length(eta) != 1 || eta <= 0) stop("eta must be a positive numeric scalar")

  # Check if sigma is symmetric
  if (!all.equal(sigma, t(sigma), tolerance = 1e-8)) {
    stop("sigma must be symmetric")
  }

  # Check positive definiteness of sigma
  det_sigma <- determinant(sigma, logarithm = TRUE)
  if (det_sigma$sign <= 0) stop("sigma must be positive definite")

  p <- length(mu)
  a <- (p - 1) / 2

  diff <- y - mu

  # Use solve(sigma, diff) instead of solve(sigma) %*% diff
  delta <- as.numeric(t(diff) %*% solve(sigma, diff))
  w <- sqrt(eta * (eta + delta))

  # Compute G term carefully to avoid division by zero
  eps <- .Machine$double.eps
  w_safe <- ifelse(w == 0, eps, w)

  G <- (eta / w_safe)^a * besselK(w_safe, a) +
    (eta / w_safe)^(p - a) * besselK(w_safe, p - a)

  # Calculate log density
  log_det_sigma <- as.numeric(det_sigma$modulus)
  d <- -log(2) - p * log(2 * pi) / 2 - 0.5 *log_det_sigma +
    (-log(besselK(eta, 0.5)) + log(G))

  if (log_density) return(d) else return(exp(d))
}
