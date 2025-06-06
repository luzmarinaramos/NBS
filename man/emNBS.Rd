% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emNBS.R
\name{emNBS}
\alias{emNBS}
\title{EM algorithm to estimate parameters of the NBS distribution}
\usage{
emNBS(y)
}
\arguments{
\item{y}{Numeric matrix of dimension \eqn{n \times p}, the observed sample.}
}
\value{
A list with elements:
\describe{
  \item{\code{mu_est}}{Estimated mean vector \eqn{\mu}.}
  \item{\code{sigma_est}}{Estimated covariance matrix \eqn{\Sigma}.}
  \item{\code{eta_est}}{Estimated shape parameter \eqn{\eta}.}
  \item{\code{loglikelihood}}{Log-likelihood of the final estimates.}
  \item{\code{convergence}}{Logical, TRUE if converged before max iterations.}
  \item{\code{iterations}}{Number of iterations performed.}
}
}
\description{
Estimates the parameters \code{mu}, \code{sigma}, and \code{eta} of the
Normal-Birnbaum-Saunders (NBS) distribution from a sample \code{y} using the EM algorithm.
}
\details{
The function uses the Expectation-Maximization (EM) algorithm to maximize the likelihood of the NBS model.
The Cholesky decomposition is used to compute matrix inverses and determinants efficiently.
Vectorized calculations improve computational speed for large datasets.

Convergence is declared when both the change in log-likelihood and the change in parameters
are below specified tolerances or the maximum number of iterations is reached.
}
\examples{
\dontrun{
set.seed(1)
y <- matrix(rnorm(100 * 3), 100, 3)  # Example data
fit <- emNBS(y)
print(fit$mu_est)
}
}
