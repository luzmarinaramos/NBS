#' Maximum Likelihood Estimator of Covariance with Fixed Mean
#'
#' This function computes the maximum likelihood estimator (MLE) of
#' the covariance matrix \eqn{\Sigma} for a multivariate normal distribution,
#' assuming that the mean vector \eqn{\mu} is fixed at a given value \eqn{\mu_0}.
#'
#' @param y A numeric matrix of size \eqn{n \times p}, where each row
#'   corresponds to an observation and each column to a variable.
#' @param mu_0 A numeric vector of length \eqn{p}, the fixed mean vector.
#'
#' @return A list with two elements:
#' \item{mu_est}{The fixed mean vector \eqn{\mu_0}.}
#' \item{sigma_est}{A \eqn{p \times p} covariance matrix, the estimated covariance matrix.}
#'
#' @details
#' The covariance is estimated as the average of the outer products of the data
#' centered at the fixed mean \eqn{\mu_0}.
#'
#' @examples
#' set.seed(123)
#' y <- MASS::mvrnorm(50, mu = c(1,2), Sigma = matrix(c(1,0.3,0.3,1),2,2))
#' mu0 <- c(1,2)
#' result <- emNMV_fixedmu(y, mu0)
#' result$mu_est
#' result$sigma_est
#'
#' @export
emNMV_fixedmu <- function(y, mu_0) {
  n <- nrow(y)
  p <- ncol(y)
  cumsum <- matrix(0, nrow = p, ncol = p)

  for (i in 1:n) {
    diff <- matrix(y[i, ], ncol = 1) - matrix(mu_0, ncol = 1)
    cumsum <- cumsum + diff %*% t(diff)
  }

  mu_est <- mu_0
  sigma_est <- cumsum / n

  return(list(mu_est = mu_est, sigma_est = sigma_est))
}
