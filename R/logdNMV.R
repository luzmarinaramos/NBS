#' Log-density of the Multivariate Normal distribution
#'
#' Computes the log-density of a multivariate normal distribution
#' for a given observation vector.
#'
#' @param y A numeric vector of length p, representing the observation.
#' @param mu A numeric vector of length p, representing the mean vector.
#' @param sigma A positive-definite covariance matrix of dimension p x p.
#'
#' @return A numeric scalar, the log-density value.
#'
#' @examples
#' mu <- c(0, 0)
#' sigma <- diag(2)
#' y <- c(1, 1)
#' logdNMV(y, mu, sigma)
#'
#' @import MASS
#' @import matrixcalc
#' @export
logdNMV <- function(y, mu, sigma) {
  sigma_inv <- solve(sigma)
  p <- length(mu)
  delta <- t(matrix(y, ncol = 1) - matrix(mu, ncol = 1)) %*%
    sigma_inv %*%
    (matrix(y, ncol = 1) - matrix(mu, ncol = 1))
  l <- (-p / 2) * log(2 * pi) - 0.5 * log(det(sigma)) - 0.5 * delta
  return(l)
}
