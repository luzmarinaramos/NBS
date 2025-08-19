#' Maximum Likelihood Estimators for Multivariate Normal
#'
#' This function computes the maximum likelihood estimators (MLEs) of
#' the mean vector \eqn{\mu} and covariance matrix \eqn{\Sigma} for a
#' multivariate normal distribution, based on a data matrix \eqn{Y}.
#'
#' @param y A numeric matrix of size \eqn{n \times p}, where each row
#'   corresponds to an observation and each column to a variable.
#'
#' @return A list with two elements:
#' \item{mu_est}{A numeric vector of length \eqn{p}, the estimated mean vector.}
#' \item{sigma_est}{A \eqn{p \times p} covariance matrix, the estimated covariance matrix.}
#'
#' @details
#' The mean is estimated as the sample mean of the data, and the covariance
#' is computed as the average of the outer products of the centered data.
#'
#' @examples
#' set.seed(123)
#' y <- MASS::mvrnorm(50, mu = c(1,2), Sigma = matrix(c(1,0.3,0.3,1),2,2))
#' result <- emNMV(y)
#' result$mu_est
#' result$sigma_est
#'
#' @export
emNMV <- function(y) {
  n <- dim(y)[1]
  mu_est <- colMeans(y)
  sigma_est <- matrix(0, ncol = dim(y)[2], nrow = dim(y)[2])

  for (i in 1:n) {
    diff <- matrix(y[i, ], ncol = 1) - matrix(mu_est, ncol = 1)
    sigma_est <- sigma_est + diff %*% t(diff)
  }

  sigma_est <- sigma_est / n
  return(list(mu_est = mu_est, sigma_est = sigma_est))
}
