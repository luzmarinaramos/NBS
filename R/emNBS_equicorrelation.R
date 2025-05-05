emNBS_equicorrelation <- function(y){
  n <- nrow(y)
  p <- ncol(y)
  a <- (p - 1) / 2
  I <- diag(p)
  J <- matrix(1,p,p)

  mu_est <- colMeans(y)
  eta_est <- 1
  S <- cov(y)
  sigma2_est <- mean(diag(S))
  rho_est <- mean(S[lower.tri(S)])/sigma2_est
  R_est <- (1-rho_est)*I+rho_est*J
  sigma_est <- sigma2_est*R_est

  loglik_est <- sum(apply(y, 1, function(yi) dNBS(yi, mu_est, sigma_est, eta_est, log_density = TRUE)))

  iter <- 0
  dif <- 1
  maxiter <- 2000

  while (dif > 1e-4 && iter < maxiter) {
    mu <- mu_est
    sigma <- sigma_est
    eta <- eta_est
    sigma_inv <- solve(sigma)

    v1 <- numeric(n)
    v2 <- numeric(n)

    loglike <- loglik_est

    for (i in 1:n) {
      yi <- y[i, ]
      delta <- t(yi - mu) %*% sigma_inv %*% (yi - mu)
      w <- sqrt(eta * (eta + delta))

      logd1 <- log(besselK(w, -a)) - p * log(2 * pi) / 2 - 0.5 * log(det(sigma)) -
        log(besselK(eta, 0.5)) + a * log(eta / w)
      logd2 <- log(besselK(w, a - p)) - p * log(2 * pi) / 2 - 0.5 * log(det(sigma)) -
        log(besselK(eta, -0.5)) - (a - p) * log(eta / w)

      max_log <- max(logd1, logd2)
      denom <- exp(logd1 - max_log) + exp(logd2 - max_log)
      py <- exp(logd1 - max_log) / denom

      v1[i] <- (1 + delta / eta)^(1/2) * (py * besselK(w, 1 - a) / besselK(w, -a) +
                                            (1 - py) * besselK(w, 1 + a - p) / besselK(w, a - p))

      v2[i] <- (1 + delta / eta)^(-1/2) * (py * besselK(w, -1 - a) / besselK(w, -a) +
                                             (1 - py) * besselK(w, -1 + a - p) / besselK(w, a - p))
    }

    mu_est <- colSums(v2 * y) / sum(v2)

    cumsum1 <- 0
    cumsum2 <- 0
    for (i in 1:n) {
      diff <- y[i, ] - mu_est
      cumsum1 <- cumsum1 + v2[i] * sum(diff^2)
      cumsum2 <- cumsum2 + v2[i] * (sum(diff))^2
    }
    rho_est <- (cumsum2-cumsum1)/((p-1)*cumsum1)
    sigma2_est <- ((1-rho_est+rho_est*p)*cumsum1-rho_est*cumsum2)/(n*p*(1-rho_est)*(1-rho_est+rho_est*p))
    R_est <- (1-rho_est)*I+rho_est*J
    sigma_est <- sigma2_est*R_est

    eta_est <- 1 / (mean(v1) + mean(v2) - 2)

    loglik_est <- sum(apply(y, 1, function(yi) dNBS(yi, mu_est, sigma_est, eta_est, log_density = TRUE)))
    dif <- loglik_est - loglike
    iter <- iter + 1
  }

  convergence <- iter < maxiter

  return(list(
    mu_est = mu_est,
    sigma_est = sigma_est,
    eta_est = eta_est,
    rho_est = rho_est,
    sigma2_est = sigma2_est,
    loglikelihood = loglik_est,
    convergence = convergence,
    iterations = iter
  ))
}
