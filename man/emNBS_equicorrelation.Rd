% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emNBS_equicorrelation.R
\name{emNBS_equicorrelation}
\alias{emNBS_equicorrelation}
\title{EM Algorithm for Multivariate Normal-Birnbaum-Saunders with Equicorrelation Structure}
\usage{
emNBS_equicorrelation(y)
}
\arguments{
\item{y}{A numeric matrix of size \eqn{n \times p}, where each row is an observation from
the multivariate NBS distribution.}
}
\value{
A list with the following components:
\describe{
  \item{\code{mu_est}}{Estimated mean vector \eqn{\mu}.}
  \item{\code{sigma_est}}{Estimated covariance matrix \eqn{\Sigma = \sigma^2 R}.}
  \item{\code{eta_est}}{Estimated shape parameter \eqn{\eta}.}
  \item{\code{rho_est}}{Estimated equicorrelation parameter \eqn{\rho}.}
  \item{\code{sigma2_est}}{Estimated common variance component \eqn{\sigma^2}.}
  \item{\code{loglikelihood}}{Log-likelihood value at convergence.}
  \item{\code{convergence}}{Logical indicating whether the algorithm converged.}
  \item{\code{iterations}}{Number of iterations until convergence.}
}
}
\description{
Estimates the parameters of the multivariate Normal-Birnbaum-Saunders (NBS) distribution
under the null hypothesis that the covariance matrix has an equicorrelation structure,
i.e., \eqn{\Sigma = \sigma^2 R} with \eqn{R = (1 - \rho)I + \rho J}.
}
\details{
This function uses the Expectation-Maximization (EM) algorithm to obtain maximum likelihood
estimates of the mean vector \eqn{\mu}, the scale parameter \eqn{\eta}, the common variance
\eqn{\sigma^2}, and the equicorrelation parameter \eqn{\rho}.


This function assumes the covariance structure is of the form \eqn{\Sigma = \sigma^2 R}, where
\eqn{R = (1 - \rho)I + \rho J} is the equicorrelation matrix. This is useful for hypothesis testing
involving equicorrelation versus a general covariance structure. The function implements safeguards
to ensure numerical stability and positive-definiteness of the estimated matrices.
}
\examples{
\dontrun{
set.seed(123)
n <- 100
p <- 3
mu <- rep(0, p)
eta <- 1
sigma2 <- 2
rho <- 0.3
R <- (1 - rho) * diag(p) + rho * matrix(1, p, p)
Sigma <- sigma2 * R
y <- rNBS(n, mu, Sigma, eta)
fit <- emNBS_equicorrelation(y)
}

}
\seealso{
\code{\link{emNBS}} for unrestricted NBS estimation.
}
