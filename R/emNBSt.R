emNBSt <- function(y){
  n <- nrow(y); p <- ncol(y)

  mu_est <- colMeans(y)
  Sigma_est <- cov(y)
  alpha_est <- 2
  nu_est <- 3

  tol <- 1e-4
  maxiter <- 200
  error <- Inf
  iter <- 0

  update_nu <- function(nu, s5, n, lower=1e-6, upper=1e6, max_expand=60){
    f <- function(x) log(x/2) + 1 - digamma(x/2) + (s5/n)
    a <- max(lower, nu/2); b <- min(upper, nu*2)
    fa <- f(a); fb <- f(b)
    k <- 0
    while(k < max_expand && is.finite(fa) && is.finite(fb) && sign(fa) == sign(fb)){
      a <- max(lower, a/2); b <- min(upper, b*2)
      fa <- f(a); fb <- f(b); k <- k + 1
    }
    if(!is.finite(fa) || !is.finite(fb) || sign(fa) == sign(fb)) return(nu)
    uniroot(f, interval=c(a,b), tol=1e-10)$root
  }

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
    nu_est <- update_nu(nu, s5, n)

    error <- sqrt(sum((mu_est-mu)^2) + sum((Sigma_est-Sigma)^2) +
                    (alpha_est-alpha)^2 + (nu_est-nu)^2)
  }

  list(mu_est=mu_est, Sigma_est=Sigma_est, alpha_est=alpha_est,
       nu_est=nu_est, iterations=iter, converged = (error <= tol))
}
