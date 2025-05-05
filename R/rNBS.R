rNBS <- function(n, mu, sigma, eta) {
  p <- length(mu)
  z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sigma)
  alpha <- 1 / sqrt(eta)
  v <- VGAM::rbisa(n, scale = 1, shape = alpha)
  y <- sweep(z, 1, sqrt(v), "*") + matrix(mu, n, p, byrow = TRUE)
  return(y)
}
