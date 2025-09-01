#' Estimation of the Agreement Probability in the Bivariate Normal Distribution
#'
#' This function estimates the agreement probability \eqn{\psi_c} for two measurements
#' that follow a bivariate Normal joint distribution.
#' The estimation is based on an observed sample \eqn{y} and a threshold value \eqn{c}.
#'
#' The function uses the EM algorithm to obtain the maximum likelihood estimators of the parameters,
#' and then computes the point estimate of \eqn{\psi_c}, its asymptotic standard error,
#' and an asymptotic confidence interval based on the normal approximation.
#'
#' @param y A matrix or data.frame of dimension \eqn{n \times 2}, containing a sample from the bivariate NBS distribution.
#' Each row corresponds to one observation of the two measurements.
#' @param c A positive threshold value \eqn{c > 0} that defines the agreement band.
#'
#' @details
#' The agreement probability \eqn{\psi_c} is defined as
#' \deqn{ \psi_c = P(|Y_1 - Y_2| \leq c) }
#' where \eqn{(Y_1, Y_2)} follows a bivariate normal distribution.
#'
#' The procedure is as follows:
#' \enumerate{
#'   \item Estimate the parameters of the bivariate normal distribution using the EM algorithm.
#'   \item Compute the mean difference and the variance of the difference.
#'   \item Obtain the point estimate of \eqn{\psi_c}.
#'   \item Derive the asymptotic standard error from the Fisher information matrix.
#'   \item Construct the 95\% asymptotic confidence interval.
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{estimate}}{Point estimate of \eqn{\psi_c}.}
#'   \item{\code{std_error}}{Asymptotic standard error.}
#'   \item{\code{lower_limit}}{Lower bound of the 95\% confidence interval.}
#'   \item{\code{upper_limit}}{Upper bound of the 95\% confidence interval.}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate a sample of size 200 from the bivariate normal distribution (illustrative example)
#' set.seed(123)
#' y <- matrix(rnorm(400), ncol = 2)
#'
#' # Compute the agreement probability with threshold c = 1
#' prob_agreement(y, c = 1)
#' }
#'
#' @seealso \code{\link{emNBS}}, \code{\link{fisher_information}}, \code{\link{pNBS0}}, \code{\link{dNBS0}}
#'
#' @export
prob_agreement_NMV <- function(y, c) {
  n <- nrow(y)

  # Estimación de parámetros por EM
  theta_est <- emNMV(y)
  mu_est <- theta_est$mu_est
  sigma_est <- theta_est$sigma_est

  # Diferencia de medias
  mu_D <- mu_est[1] - mu_est[2]

  # Varianza y desviación estándar de la diferencia
  sigma2_D <- sigma_est[1,1] + sigma_est[2,2] - 2 * sigma_est[1,2]   # ERROR: antes usabas sigma
  sigma_D  <- sqrt(sigma2_D)

  # Estimador de phi_c
  psi_est <- pnorm((c - mu_D) / sigma_D) -
    pnorm((-c - mu_D) / sigma_D)

  # Información de Fisher
  FI_est <- fisher_information_NMV(mu_est, sigma_est)

  # Derivadas parciales para el Jacobiano
  psi_mu <- ( dnorm((c - mu_D) / sigma_D) -
                dnorm((-c - mu_D) / sigma_D) ) * c(1, -1)

  # Revisa esta parte: parece que intentas derivar respecto a sigma,
  # pero la fórmula original no estaba clara. Te la dejo estructurada:
  psi_sigma <- ( -( (c - mu_D) / sigma_D^2 ) * dnorm((c - mu_D) / sigma_D) +
                   ( (-c - mu_D) / sigma_D^2 ) * dnorm((-c - mu_D) / sigma_D) ) *
    (1 / (2 * sigma_D)) * c(1, -2, 1)

  # Jacobiano completo
  jacobiano <- c(psi_mu, psi_sigma)

  # Varianza asintótica
  V <- t(jacobiano) %*% solve(FI_est) %*% jacobiano   # ERROR: antes pusiste jacobiano^T
  sd <- sqrt(V / n)

  # Intervalo de confianza
  li <- psi_est - qnorm(0.975) * sd
  ls <- psi_est + qnorm(0.975) * sd

  result <- list(
    estimate = psi_est,
    error_estandar = sd,
    limite_inferior = li,
    limite_superior = ls
  )

  return(result)
}
