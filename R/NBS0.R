#' Density, CDF, and CDF derivative of the standardized univariate Normal-Birnbaum-Saunders (NBS) distribution
#'
#' These functions compute the density, cumulative distribution function (CDF),
#' and derivative of the CDF with respect to the shape parameter \eqn{\eta},
#' for the standardized univariate Normal-Birnbaum-Saunders distribution
#' (i.e., \eqn{\mu = 0}, \eqn{\sigma = 1}).
#'
#' @param y Numeric. Evaluation point.
#' @param eta Positive numeric. Shape parameter of the standardized NBS distribution.
#'
#' @details
#' - `dNBS0(y, eta)` returns the probability density function:
#'   \deqn{f(y;\eta) = \frac{1}{2\sqrt{2\pi}K_{0.5}(\eta)} \left[ K_0(w) + \frac{\eta}{w} K_1(w)\right],}
#'   where \eqn{w = \sqrt{\eta(\eta+y^2)}} and \eqn{K_\nu(\cdot)} is the modified Bessel function of the third kind.
#'
#' - `pNBS0(y, eta)` returns the cumulative distribution function, computed numerically as:
#'   \deqn{F(y;\eta) = \frac{\sqrt{\eta}}{2} \int_0^\infty \left( v^{-1/2} + v^{-3/2} \right)
#'   \phi\left( \sqrt{\eta}(v^{1/2} - v^{-1/2}) \right) \Phi\left( \frac{y}{\sqrt{v}} \right) dv.}
#'
#' - `dpNBS0(y, eta)` returns the derivative of the CDF with respect to \eqn{\eta}, given by
#'   \deqn{\frac{\partial}{\partial \eta} F(y;\eta) = \frac{1}{2\eta}F(y;\eta) +
#'   \frac{1}{4} \int_0^\infty (v^{-1/2}+v^{-3/2})(v^{1/2}-v^{-1/2})
#'   \phi(\sqrt{\eta}(v^{1/2}-v^{-1/2})) \Phi(y/\sqrt{v}) \, dv.}
#'
#' @return
#' A numeric value:
#' - For `dNBS0()`, the density at point `y`.
#' - For `pNBS0()`, the CDF at point `y`.
#' - For `dpNBS0()`, the derivative of the CDF with respect to `eta` at point `y`.
#'
#' @examples
#' dNBS0(0, 1)
#' pNBS0(0, 1)
#' dpNBS0(0, 1)
#'
#' @export
dNBS0 <- function(y, eta) {
  delta <- y^2
  w <- sqrt(eta * (eta + delta))
  d <- (1 / (2 * sqrt(2 * pi) * besselK(eta, 0.5))) *
    (besselK(w, 0) + (eta / w) * besselK(w, 1))
  return(d)
}

#' @rdname dNBS0
#' @export
pNBS0 <- function(y, eta) {
  integrand <- function(v) {
    term1 <- v^(-1/2) + v^(-3/2)
    phi_val <- dnorm(sqrt(eta) * (sqrt(v) - 1/sqrt(v)))
    Phi_val <- pnorm(y / sqrt(v))
    term1 * phi_val * Phi_val
  }
  result <- (sqrt(eta) / 2) * integrate(
    integrand, lower = 0, upper = Inf,
    rel.tol = 1e-8, abs.tol = 1e-10
  )$value
  return(result)
}

#' @rdname dNBS0
#' @export
dpNBS0 <- function(y, eta) {
  integrand <- function(v) {
    term1 <- v^(-1/2) + v^(-3/2)
    term2 <- v^(1/2) - v^(-1/2)
    phi_val <- dnorm(sqrt(eta) * (sqrt(v) - 1/sqrt(v)))
    Phi_val <- pnorm(y / sqrt(v))
    term1 * term2^2 * phi_val * Phi_val
  }
  result1 <- (1 / (2 * eta)) * pNBS0(y, eta)
  result2 <- (sqrt(eta) / 4) * integrate(
    integrand, lower = 0, upper = Inf,
    rel.tol = 1e-8, abs.tol = 1e-10
  )$value
  result <- result1 - result2
  return(result)
}
