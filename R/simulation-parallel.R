run_simulation_parallel <- function(n, B,mu,sigma,eta,significancia,rho) {
library(parallel)

n_cores <- detectCores() - 1  # deja un core libre

# Función que ejecuta 1 réplica completa (generar muestra, test, captura errores)
run_one_sim <- function(seed) {
  set.seed(seed)
  repeat {
    y <- rNBS(n, mu, sigma, eta)
    res <- tryCatch({
      test_res <- testNBS_equicorrelation(y)
      test_res$p_value < significancia
    }, error = function(e) NULL)
    if (!is.null(res)) return(as.numeric(res))  # vector lógico convertido a numérico
  }
}

# Semillas para reproducibilidad
seeds <- sample.int(1e7, B)

# Paralelización
cl <- makeCluster(n_cores)
clusterExport(cl, c("n", "mu", "sigma", "eta", "significancia", "rNBS", "testNBS_equicorrelation"))
results <- parLapply(cl, seeds, run_one_sim)
stopCluster(cl)

# Convertir lista a matriz
results_mat <- do.call(rbind, results)

# Proporción de rechazos para cada test
colMeans(results_mat)

}
