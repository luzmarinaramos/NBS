% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dNBS.R
\name{dNBS}
\alias{dNBS}
\title{Density function of the Multivariate Normal Birnbaum-Saunders (NBS) distribution}
\usage{
dNBS(y, mu, sigma, eta, log_density = FALSE)
}
\arguments{
\item{y}{Numeric vector. Point at which to evaluate the density.}

\item{mu}{Numeric vector. Mean vector of the distribution.}

\item{sigma}{Numeric matrix. Covariance matrix (must be symmetric positive definite).}

\item{eta}{Positive numeric scalar. Shape parameter of the Birnbaum-Saunders mixture.}

\item{log_density}{Logical; if TRUE, returns the log-density. Default is FALSE.}
}
\value{
Numeric scalar. Density (or log-density) evaluated at \code{y}.
}
\description{
Computes the density (or log-density) of the Multivariate Normal Birnbaum-Saunders distribution
at a given point \code{y}.
}
\examples{
mu <- c(0, 0)
sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
eta <- 1
y <- c(1, 1)
dNBS(y, mu, sigma, eta)

}
