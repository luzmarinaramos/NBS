emNBSt_nu_fixed <- function(y,nu){
  n <- nrow(y); p <- ncol(y)

  mu_est <- colMeans(y)
  Sigma_est <- cov(y)
  alpha_est <- 2
  nu_est <- nu

  tol <- 1e-4
  maxiter <- 200
  error <- Inf
  iter <- 0


  while(error > tol && iter < maxiter){
    iter <- iter + 1
    mu <- mu_est; Sigma <- Sigma_est; alpha <- alpha_est; nu <- nu_est

    invSigma <- solve(Sigma)

    s1 <- rep(0,p); s2 <- 0; s4 <- 0; s5 <- 0
    w1 <- numeric(n)

    for(i in 1:n){
      yi <- as.numeric(y[i,])
      di <- as.numeric(t(yi-mu) %*% invSigma %*% (yi-mu))
      w <- pesos(di, p, alpha, nu)   # asumo vector length 5
      w1[i] <- w[1]
      s1 <- s1 + w1[i]*yi
      s2 <- s2 + w1[i]
      s4 <- s4 + (w[2] + w[3] - 2)
      s5 <- s5 + (w[4] - w[5])
    }

    mu_est <- s1/s2

    # Sigma centrada en mu_est (nuevo)
    S <- matrix(0, p, p)
    for(i in 1:n){
      zi <- as.numeric(y[i,]) - mu_est
      S <- S + w1[i] * tcrossprod(zi)
    }
    Sigma_est <- S / n

    alpha_est <- sqrt(max(s4/n, 1e-12))
    nu_est <- nu

    error <- sqrt(sum((mu_est-mu)^2) + sum((Sigma_est-Sigma)^2) +
                    (alpha_est-alpha)^2 + (nu_est-nu)^2)
  }

  list(mu_est=mu_est, Sigma_est=Sigma_est, alpha_est=alpha_est,
       nu_est=nu_est, iterations=iter, converged = (error <= tol))
}
