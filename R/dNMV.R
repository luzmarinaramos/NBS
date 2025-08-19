#' Density of the Multivariate Normal distribution
#'
#' Computes the probability density of a multivariate normal distribution
#' for a given observation vector.
#'
#' @param y A numeric vector of length p, representing the observation.
#' @param mu A numeric vector of length p, representing the mean vector.
#' @param sigma A positive-definite covariance matrix of dimension p x p.
#'
#' @return A numeric scalar, the density value.
#'
#' @examples
#' mu <- c(0, 0)
#' sigma <- diag(2)
#' y <- c(1, 1)
#' dNMV(y, mu, sigma)
#'
#' @export
dNMV <- function(y, mu, sigma) {
  exp(logdNMV(y, mu, sigma))
}
