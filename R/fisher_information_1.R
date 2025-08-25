#' Fisher Information Matrix for the NBS Distribution
#'
#' Computes the Fisher Information matrix for the Normal-Birnbaum-Saunders (NBS) distribution
#' given a sample \eqn{y}, mean vector \eqn{\mu}, covariance matrix \eqn{\Sigma},
#' and shape parameter \eqn{\eta}.
#'
#' @param y A numeric matrix of dimension \eqn{n \times p}, where each row corresponds to an observation
#'          and each column to a variable.
#' @param mu A numeric vector of length \eqn{p}, representing the mean vector of the distribution.
#' @param sigma A \eqn{p \times p} positive-definite symmetric covariance matrix.
#' @param eta A positive numeric scalar, the shape parameter of the NBS distribution.
#'
#' @details
#' The function evaluates the Fisher Information matrix by computing the score contributions
#' for each observation and aggregating them. It uses Bessel functions of the second kind
#' (`besselK`) and their derivatives to handle the NBS density structure.
#'
#' Internally, auxiliary functions \eqn{G(\eta,\delta)}, \eqn{\partial G / \partial \delta},
#' and \eqn{\partial G / \partial \eta} are used, where \eqn{\delta = (y_i - \mu)^T \Sigma^{-1} (y_i - \mu)}.
#'
#' @return A Fisher Information matrix of dimension \eqn{(p + q + 1) \times (p + q + 1)},
#' where \eqn{q = p(p+1)/2}.
#'
#' @references
#' - Birnbaum, Z. W., & Saunders, S. C. (1969). A new family of life distributions.
#'   *Journal of Applied Probability*, 6(2), 319–327.
#' - Fang, K. T., Kotz, S., & Ng, K. W. (1990). *Symmetric Multivariate and Related Distributions*.
#'   Chapman and Hall.
#' - Gupta, A. K., & Varga, T. (1993). *Elliptically Contoured Models in Statistics*.
#'   Springer.
#'
#' @examples
#' set.seed(123)
#' y <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' mu <- rep(0, 5)
#' sigma <- diag(5)
#' eta <- 2
#' fisher_information_1(y, mu, sigma, eta)
#'
#' @seealso \code{\link{scoreNBS}}, \code{\link{emNBS}}
#'
#' @export
fisher_information_1 <- function(y, mu, sigma, eta) {
  # Input validations
  if (!is.matrix(y)) stop("Error: 'y' must be a matrix.")
  if (length(mu) != ncol(y)) stop("Error: Incompatible dimensions between 'y' and 'mu'.")
  if (!isSymmetric(sigma)) stop("Error: 'sigma' must be a symmetric matrix.")
  if (!is.positive.definite(sigma)) stop("Error: 'sigma' must be positive definite.")
  if (eta <= 0) stop("Error: 'eta' must be positive.")

  n <- nrow(y)
  p <- length(mu)
  q <- p * (p + 1) / 2
  a <- (p - 1) / 2
  sigma_inv <- solve(sigma)
  D_p <- duplication.matrix(p)

  # Pre-calculate constant values of the Bessel function for eta
  bessel_eta <- besselK(eta, 0.5)
  dbessel_eta <- -0.5 * (besselK(eta, 1.5) + besselK(eta, -0.5))

  # Derivative of the besselK function
  dbesselK <- function(u, d) {
    -0.5 * (besselK(u, d + 1) + besselK(u, d - 1))
  }

  # Auxiliary functions G, G_delta, G_eta
  G <- function(eta, delta) {
    w <- sqrt(eta * (eta + delta))
    (eta / w)^a * besselK(w, a) + (eta / w)^(p - a) * besselK(w, p - a)
  }

  G_delta <- function(eta, delta) {
    w <- sqrt(eta * (eta + delta))
    w_delta <- eta / (2 * w)
    t1 <- (eta / w)^a * (dbesselK(w, a) - a * besselK(w, a) / w)
    t2 <- (eta / w)^(p - a) * (dbesselK(w, p - a) - (p - a) * besselK(w, p - a) / w)
    (t1 + t2) * w_delta
  }

  G_eta <- function(eta, delta) {
    w <- sqrt(eta * (eta + delta))
    w_eta <- (2 * eta + delta) / (2 * w)
    t1 <- (eta / w)^a * (w_eta * dbesselK(w, a) + a * (1 / eta - w_eta / w) * besselK(w, a))
    t2 <- (eta / w)^(p - a) * (w_eta * dbesselK(w, p - a) + (p - a) * (1 / eta - w_eta / w) * besselK(w, p - a))
    t1 + t2
  }

  # Initialize fisher information matrix
  var_U <- matrix(0,nrow=(p+q+1),ncol=(p+q+1))
  # Loop to calculate components for each observation
  for (i in 1:n) {
    yi <- y[i, ]
    delta <- as.numeric(t(yi - mu) %*% sigma_inv %*% (yi - mu))
    G_val <- G(eta, delta)
    W_delta <- G_delta(eta, delta) / G_val
    W_eta <- G_eta(eta, delta) / G_val - dbessel_eta / bessel_eta
    U_mu <-  - 2 * W_delta * (sigma_inv %*% (yi - mu))
    M <- 2 * W_delta * sigma_inv %*% (yi - mu) %*% t(yi - mu) %*% sigma_inv + sigma_inv
    U_phi <-  - 0.5 * t(D_p) %*% as.vector(M)
    U_eta <-   W_eta
    U <- c(U_mu, U_phi, U_eta)
    var_U <- var_U + U%*%t(U)
  }
  I <- var_U/n

  # Return the concatenated score vector
  return(I)
}
