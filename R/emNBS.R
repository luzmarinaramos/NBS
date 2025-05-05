emNBS <- function(y){
  n <- nrow(y)
  p <- ncol(y)
  a <- (p - 1) / 2

  mu_est <- colMeans(y)
  eta_est <- 1
  sigma_est <- cov(y) / (1 + 1 / (2 * eta_est))

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

    sigma_est <- matrix(0, p, p)
    for (i in 1:n) {
      diff <- y[i, ] - mu_est
      sigma_est <- sigma_est + v2[i] * tcrossprod(diff)
    }
    sigma_est <- sigma_est / n

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
    loglikelihood = loglik_est,
    convergence = convergence,
    iterations = iter
  ))
}
