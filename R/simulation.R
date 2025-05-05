set.seed(NULL)
#nn <- c(60, 200, 800)
nn <- c(60)
#pp <- c(5, 15, 45)
pp <- c(45)
B <- 100  # número de simulaciones
significancia <- 0.05

results <- array(0, dim = c(length(nn), length(pp), 4),
                 dimnames = list(
                   n = as.character(nn),
                   p = as.character(pp),
                   test = c("LR", "Wald", "Score", "Gradiente")
                 ))

for (i in seq_along(nn)) {
  for (j in seq_along(pp)) {
    n <- nn[i]
    p <- pp[j]

    mu <- rep(1, p)
    I <- diag(p)
    J <- matrix(1, p, p)
    rho <- 0.5
    S2 <- 1
    sigma <- S2^2 * ((1 - rho) * I + rho * J)
    eta <- 0.5

    rechazo <- matrix(0, nrow = B, ncol = 4)

    for (b in 1:B) {
      y <- rNBS(n, mu, sigma, eta)
      test1 <- testNBS_equicorrelation(y)
      t1 <- test1$p_value[1] < significancia
      t2 <- test1$p_value[2] < significancia
      t3 <- test1$p_value[3] < significancia
      t4 <- test1$p_value[4] < significancia
      rechazo[b, ] <- c(t1, t2, t3, t4)
    }

    results[i, j, ] <- colMeans(rechazo)
  }
}

