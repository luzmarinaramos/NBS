#' Score function for the Normal Birnbaum-Saunders (NBS) distribution
#'
#' Computes the score vector (gradient of the log-likelihood) for a sample \code{y}
#' drawn from a Normal Birnbaum-Saunders (NBS) distribution with parameters \code{mu}, \code{sigma}, and \code{eta}.
#'
#' @param y A numeric matrix of observed data, where each row is an observation and each column corresponds to a variable.
#' @param mu A numeric vector representing the location parameter of the NBS distribution.
#' @param sigma A positive definite symmetric matrix representing the scale (covariance) parameter.
#' @param eta A positive numeric scalar representing the shape parameter.
#'
#' @return A numeric vector representing the concatenated score vector with respect to \code{mu}, \code{phi = vech(sigma)} (the half-vectorization of \code{sigma}), and \code{eta}.
#'
#' @examples
#' \dontrun{
#' library(matrixcalc)
#' set.seed(1)
#' y <- matrix(rnorm(30), nrow = 10, ncol = 3)
#' mu <- rep(0, 3)
#' sigma <- diag(3)
#' eta <- 1
#' scoreNBS(y, mu, sigma, eta)
#' }
#'
#' @importFrom matrixcalc duplication.matrix is.positive.definite
#' @export

scoreNBS <- function(y, mu, sigma, eta) {
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

  # Initialize components of the score vector
  U_mu <- numeric(p)
  U_phi <- numeric(q)
  U_eta <- 0

  # Loop to calculate components for each observation
  for (i in 1:n) {
    yi <- y[i, ]
    delta <- as.numeric(t(yi - mu) %*% sigma_inv %*% (yi - mu))
    G_val <- G(eta, delta)
    W_delta <- G_delta(eta, delta) / G_val
    W_eta <- G_eta(eta, delta) / G_val - dbessel_eta / bessel_eta

    U_mu <- U_mu - 2 * W_delta * (sigma_inv %*% (yi - mu))
    M <- 2 * W_delta * sigma_inv %*% (yi - mu) %*% t(yi - mu) %*% sigma_inv + sigma_inv
    U_phi <- U_phi - 0.5 * t(D_p) %*% as.vector(M)
    U_eta <- U_eta + W_eta
  }

  # Return the concatenated score vector
  c(U_mu, U_phi, U_eta)
}
