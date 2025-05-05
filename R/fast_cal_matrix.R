# matrix functions
commutation_matrix <- function(p) {
  K_p <- matrix(0, nrow = p^2, ncol = p^2)
  for(i in 1:p) {
    for(j in 1:p) {
      K_p[(i-1)*p + j, (j-1)*p + i] <- 1
    }
  }
  return(K_p)
}
# vec <- function(A) as.vector(A)

