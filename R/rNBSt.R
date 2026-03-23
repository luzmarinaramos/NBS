#' Random sampling from the NBSt distribution
#'
#' Generates \code{n} random samples from the Multivariate Normal Birnbaum--Saunders--t
#' (NBSt) distribution using the hierarchical representation
#' \deqn{
#' Y_i \mid U_i, V_i \sim \mathcal{N}_p\!\left(\mu,\ \frac{U_i}{V_i}\Sigma\right),\quad
#' U_i \sim \mathrm{BS}(\alpha, 1),\quad
#' V_i \sim \mathrm{Gamma}(\nu/2,\nu/2),
#' }
#' where the Gamma distribution is parameterized with \emph{shape} \eqn{\nu/2} and
#' \emph{rate} \eqn{\nu/2}.
#'
#' @param n Integer. Number of random samples to generate.
#' @param mu Numeric vector of length \eqn{p}. Location vector.
#' @param Sigma Numeric \eqn{p \times p} matrix. Scale matrix (must be symmetric
#'   positive definite).
#' @param alpha Positive numeric scalar. Birnbaum--Saunders shape parameter
#'   (with scale fixed at 1).
#' @param nu Positive numeric scalar. Degrees of freedom of the t component.
#'
#' @return A numeric matrix with \code{n} rows and \code{length(mu)} columns. Each row
#'   is a draw from the NBSt distribution.
#'
#' @examples
#' set.seed(123)
#' mu <- c(0, 0)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' alpha <- 0.8
#' nu <- 5
#' rNBSt(5, mu, Sigma, alpha, nu)
#'
#' @importFrom MASS mvrnorm
#' @importFrom VGAM rbisa
#' @export
rNBSt <- function(n, mu, Sigma, alpha, nu) {

  p <- length(mu)

  U <- VGAM::rbisa(n, scale = 1, shape = alpha)

  V <- stats::rgamma(n, shape = nu / 2, rate = nu / 2)

  Z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

  s <- sqrt(U / V)
  Y <- sweep(Z, 1, s, "*") + matrix(mu, n, p, byrow = TRUE)

  Y
}
