dNBS <- function(y, mu, sigma, eta, log_density = FALSE) {
  p <- length(mu)
  a <- (p - 1) / 2
  delta <- t(y-mu)%*%solve(sigma)%*%(y-mu)
  w <- sqrt(eta * (eta + delta))

  G <- (eta / w)^a * besselK(w, a) +
    (eta / w)^(p - a) * besselK(w, p - a)

  d <- -log(2) - p * log(2 * pi) / 2 -0.5 * log(det(sigma)) -
    log(besselK(eta, 0.5)) +log(G)

  if (log_density) {return(d)} else {return(exp(d))}
}


