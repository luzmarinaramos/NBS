#' Maximum Likelihood Estimator of Mean with Fixed Covariance
#'
#' This function computes the maximum likelihood estimator (MLE) of
#' the mean vector \eqn{\mu} for a multivariate normal distribution,
#' assuming that the covariance matrix \eqn{\Sigma} is fixed at a given value \eqn{\Sigma_0}.
#'
#' @param y A numeric matrix of size \eqn{n \times p}, where each row
#'   corresponds to an observation and each column to a variable.
#' @param sigma_0 A numeric positive-definite matrix of size \eqn{p \times p},
#'   the fixed covariance matrix.
#'
#' @return A list with two elements:
#' \item{mu_est}{The estimated mean vector \eqn{\mu}.}
#' \item{sigma_est}{The fixed covariance matrix \eqn{\Sigma_0}.}
#'
#' @details
#' The mean vector is estimated as the sample mean. The covariance matrix is kept
#' fixed at the provided value \eqn{\Sigma_0}.
#'
#' @examples
#' set.seed(123)
#' y <- MASS::mvrnorm(50, mu = c(1,2), Sigma = matrix(c(1,0.3,0.3,1),2,2))
#' Sigma0 <- diag(2)
#' result <- emNMV_fixedsigma(y, Sigma0)
#' result$mu_est
#' result$sigma_est
#'
#' @export
emNMV_fixedsigma <- function(y, sigma_0) {
  mu_est <- colMeans(y)
  sigma_est <- sigma_0

  return(list(mu_est = mu_est, sigma_est = sigma_est))
}
