#' Fisher Information Matrix for the Multivariate Normal
#'
#' Computes the Fisher information matrix of a multivariate normal distribution
#' with mean vector `mu` and covariance matrix `sigma`.
#'
#' @param mu A numeric vector of length p (mean vector).
#' @param sigma A symmetric, positive-definite covariance matrix of dimension p x p.
#'
#' @return A list with:
#' \item{I_mu_mu}{Block corresponding to the mean parameters.}
#' \item{I_phi_phi}{Block corresponding to the covariance parameters (vech(sigma)).}
#' \item{Information}{The full Fisher information matrix.}
#'
#' @importFrom matrixcalc duplication.matrix
#' @export
fisher_information_NMV <- function(mu, sigma) {
  p <- length(mu)
  q <- p * (p + 1) / 2

  # checks
  if (!isSymmetric(sigma)) stop("Sigma must be symmetric")
  if (any(eigen(sigma, only.values = TRUE)$values <= 0))
    stop("Sigma must be positive definite")

  D_p <- duplication.matrix(p)
  sigma_inv <- solve(sigma)

  I_mu_mu <- sigma_inv
  I_phi_phi <- (1 / 2) * t(D_p) %*% kronecker(sigma_inv, sigma_inv) %*% D_p

  zero_12 <- matrix(0, nrow = p, ncol = q)
  zero_21 <- matrix(0, nrow = q, ncol = p)

  top_row <- cbind(I_mu_mu, zero_12)
  bottom_row <- cbind(zero_21, I_phi_phi)
  Information <- rbind(top_row, bottom_row)

  return(list(I_mu_mu = I_mu_mu,
              I_phi_phi = I_phi_phi,
              Information = Information))
}
