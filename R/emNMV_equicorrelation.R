#' Estimation under Equicorrelation Structure
#'
#' This function computes the maximum likelihood estimators (MLE)
#' of the mean vector \eqn{\mu} and the covariance matrix \eqn{\Sigma}
#' of a multivariate normal distribution, under the constraint that
#' \eqn{\Sigma} has an equicorrelation structure.
#'
#' @param y A numeric matrix of size \eqn{n \times p}, where each row
#'   corresponds to an observation and each column to a variable.
#'
#' @return A list with two elements:
#' \item{mu_est}{The estimated mean vector \eqn{\mu}.}
#' \item{sigma_est}{The estimated covariance matrix \eqn{\Sigma} with equicorrelation structure.}
#'
#' @details
#' The covariance matrix is assumed to have the form:
#' \deqn{\Sigma = \sigma^2 \left[ (1-\rho) I_p + \rho J_p \right],}
#' where \eqn{I_p} is the identity matrix, \eqn{J_p} is the all-ones matrix,
#' \eqn{\rho} is the correlation parameter, and \eqn{\sigma^2} is a variance parameter.
#'
#' The estimators are obtained by maximizing the likelihood under this restriction.
#'
#' @examples
#' set.seed(123)
#' Sigma <- 2 * ((1 - 0.5) * diag(3) + 0.5 * matrix(1, 3, 3))
#' y <- MASS::mvrnorm(50, mu = c(1,2,3), Sigma = Sigma)
#' result <- emNMV_equicorrelation(y)
#' result$mu_est
#' result$sigma_est
#'
#' @export
emNMV_equicorrelation <- function(y) {
  n <- nrow(y)
  p <- ncol(y)
  I <- diag(p)
  J <- matrix(1, nrow = p, ncol = p)

  mu_est <- colMeans(y)

  # Compute quadratic forms
  cumsum1 <- 0
  cumsum2 <- 0
  for (i in 1:n) {
    diff <- matrix(y[i, ] - mu_est, nrow = 1)
    cumsum1 <- cumsum1 + diff %*% t(diff)
    cumsum2 <- cumsum2 + diff %*% J %*% t(diff)
  }

  cumsum1 <- as.numeric(cumsum1)
  cumsum2 <- as.numeric(cumsum2)

  # Estimate rho
  rho_est <- (cumsum2 - cumsum1) / ((p - 1) * cumsum1)
  R_est <- (1 - rho_est) * I + rho_est * J
  R_est_inv <- solve(R_est)

  # Estimate sigma^2
  cumsum3 <- 0
  for (i in 1:n) {
    diff <- matrix(y[i, ] - mu_est, nrow = 1)
    cumsum3 <- cumsum3 + diff %*% R_est_inv %*% t(diff)
  }
  cumsum3 <- as.numeric(cumsum3)
  sigma2_est <- cumsum3 / (n * p)

  sigma_est <- sigma2_est * R_est

  return(list(mu_est = mu_est, sigma_est = sigma_est))
}
