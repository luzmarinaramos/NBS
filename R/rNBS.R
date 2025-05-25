
#' Random sampling from the Multivariate Normal Birnbaum-Saunders (NBS) distribution
#'
#' Generates \code{n} random samples from the Multivariate Normal Birnbaum-Saunders distribution.
#'
#' @param n Integer. Number of random samples to generate.
#' @param mu Numeric vector. Mean vector of the multivariate normal component.
#' @param sigma Numeric matrix. Covariance matrix (must be symmetric positive definite).
#' @param eta Positive numeric scalar. Shape parameter of the Birnbaum-Saunders distribution.
#'
#' @return A numeric matrix with \code{n} rows and \code{length(mu)} columns, where each row is a random sample from the multivariate NBS distribution.
#'
#' @examples
#' mu <- c(0, 0)
#' sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' eta <- 1
#' rNBS(5, mu, sigma, eta)
#'
#' @importFrom MASS mvrnorm
#' @importFrom VGAM rbisa
#' @export
rNBS <- function(n, mu, sigma, eta) {
  # Input checks
  if (!is.numeric(n) || length(n) != 1 || n <= 0) stop("n must be a positive integer")
  if (!is.numeric(mu)) stop("mu must be a numeric vector")
  if (!is.matrix(sigma)) stop("sigma must be a numeric matrix")
  if (nrow(sigma) != ncol(sigma)) stop("sigma must be a square matrix")
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
  # Generate MVN(0, sigma) samples
  z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sigma)

  # Calculate shape parameter for Birnbaum-Saunders mixing variable
  alpha <- 1 / sqrt(eta)

  # Generate Birnbaum-Saunders samples for scaling
  v <- VGAM::rbisa(n, scale = 1, shape = alpha)

  # Scale rows of z by sqrt(v) and add mu
  y <- sweep(z, 1, sqrt(v), "*") + matrix(mu, n, p, byrow = TRUE)

  return(y)
}
