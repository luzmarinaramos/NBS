make_pd <- function(sigma, epsilon = 1e-8) {
  try_chol <- try(chol(sigma), silent = TRUE)

  if (inherits(try_chol, "try-error")) {
    sigma <- sigma + epsilon * diag(ncol(sigma))
  }
  return(sigma)
}
