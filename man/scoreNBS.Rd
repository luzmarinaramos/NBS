% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreNBS.R
\name{scoreNBS}
\alias{scoreNBS}
\title{Score function for the Normal Birnbaum-Saunders (NBS) distribution}
\usage{
scoreNBS(y, mu, sigma, eta)
}
\arguments{
\item{y}{A numeric matrix of observed data, where each row is an observation and each column corresponds to a variable.}

\item{mu}{A numeric vector representing the location parameter of the NBS distribution.}

\item{sigma}{A positive definite symmetric matrix representing the scale (covariance) parameter.}

\item{eta}{A positive numeric scalar representing the shape parameter.}
}
\value{
A numeric vector representing the concatenated score vector with respect to \code{mu}, \code{phi = vech(sigma)} (the half-vectorization of \code{sigma}), and \code{eta}.
}
\description{
Computes the score vector (gradient of the log-likelihood) for a sample \code{y}
drawn from a Normal Birnbaum-Saunders (NBS) distribution with parameters \code{mu}, \code{sigma}, and \code{eta}.
}
\examples{
\dontrun{
library(matrixcalc)
set.seed(1)
y <- matrix(rnorm(30), nrow = 10, ncol = 3)
mu <- rep(0, 3)
sigma <- diag(3)
eta <- 1
scoreNBS(y, mu, sigma, eta)
}

}
