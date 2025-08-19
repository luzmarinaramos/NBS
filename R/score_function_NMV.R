#' Score Function for the Multivariate Normal Distribution
#'
#' Computes the score vector for the parameters \eqn{(\mu, \Sigma)} of a
#' multivariate normal distribution given a sample.
#'
#' @param y A numeric matrix of size \eqn{n \times p}, where each row is an
#'   observation from a \eqn{p}-dimensional multivariate normal distribution.
#' @param mu A numeric vector of length \eqn{p}, the mean vector.
#' @param sigma A \eqn{p \times p} positive definite covariance matrix.
#'
#' @return A numeric vector containing the score contributions:
#'   - First \eqn{p} elements: score with respect to \eqn{\mu}.
#'   - Next \eqn{p(p+1)/2} elements: score with respect to the unique elements
#'     of \eqn{\Sigma}, vectorized with \code{vech()}.
#'
#' @details
#' The score is defined as the gradient of the log-likelihood with respect to
#' the parameters. For the multivariate normal, the score components are:
#' \deqn{U_\mu = \sum_{i=1}^n \Sigma^{-1}(y_i - \mu),}
#' \deqn{U_\Sigma = \tfrac{1}{2}\sum_{i=1}^n \mathrm{vech}\big[
#' \Sigma^{-1}(y_i-\mu)(y_i-\mu)^T \Sigma^{-1} - \Sigma^{-1}\big].}
#'
#' @examples
#' set.seed(123)
#' y <- matrix(rnorm(20), ncol = 2)
#' mu <- c(0, 0)
#' sigma <- diag(2)
#' score_function_NMV(y, mu, sigma)
#'
#' @references
#' Magnus, J.R., & Neudecker, H. (1999). \emph{Matrix Differential Calculus with
#' Applications in Statistics and Econometrics}. Wiley.
#'
#' @export
score_function_NMV <- function(y, mu, sigma) {
  n <- nrow(y)
  p <- length(mu)
  q <- p * (p + 1) / 2

  # More stable inversion
  sigma_inv <- chol2inv(chol(sigma))

  # Centered data
  y_centered <- sweep(y, 2, mu)

  U_mu <- colSums(y_centered %*% sigma_inv)

  U_phi <- rep(0, q)
  for (i in 1:n) {
    yi <- matrix(y_centered[i, ], ncol = 1)
    U_phi <- U_phi + 0.5 * vech(sigma_inv %*% yi %*% t(yi) %*% sigma_inv - sigma_inv)
  }

  U <- c(U_mu, U_phi)
  return(U)
}
