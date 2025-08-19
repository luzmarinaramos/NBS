#' Log-likelihood of the Multivariate Normal distribution
#'
#' Computes the log-likelihood of a sample assumed to come from
#' a multivariate normal distribution with parameters \eqn{mu} and \eqn{sigma}.
#'
#' @param y A numeric matrix of dimension n x p, where each row represents an observation.
#' @param mu A numeric vector of length p, representing the mean vector.
#' @param sigma A positive-definite covariance matrix of dimension p x p.
#'
#' @return A numeric scalar, the log-likelihood value.
#'
#' @examples
#' mu <- c(0, 0)
#' sigma <- diag(2)
#' y <- matrix(rnorm(20), ncol = 2)
#' loglikeNMV(y, mu, sigma)
#'
#' @export
loglikeNMV <- function(y, mu, sigma) {
  n <- nrow(y)
  logvero <- 0
  for (i in 1:n) {
    yi <- unlist(y[i, ])
    logvero <- logvero + logdNMV(yi, mu, sigma)
  }
  return(logvero)
}
