#' Generate the commutation matrix of size \eqn{p^2 \times p^2}
#'
#' Constructs the commutation matrix \eqn{K_p}, which is a
#' \eqn{p^2 \times p^2} permutation matrix satisfying
#' \eqn{\mathrm{vec}(A^\top) = K_p \mathrm{vec}(A)}
#' for any \eqn{p \times p} matrix \eqn{A}.
#'
#' @param p Positive integer specifying the dimension of the square matrix.
#'
#' @return A \eqn{p^2 \times p^2} commutation matrix.
#'
#' @examples
#' K2 <- commutation_matrix(2)
#' print(K2)
#'
#' @export
commutation_matrix <- function(p) {
  if (!(is.numeric(p) && length(p) == 1 && p == floor(p) && p > 0)) {
    stop("p must be a positive integer")
  }

  K_p <- matrix(0, nrow = p^2, ncol = p^2)
  for (i in 1:p) {
    for (j in 1:p) {
      K_p[(i - 1) * p + j, (j - 1) * p + i] <- 1
    }
  }
  return(K_p)
}
duplication_matrix_safe <- function(n) {
  if (n == 1) {
    return(matrix(1, nrow = 1, ncol = 1))
  } else {
    return(matrixcalc::duplication.matrix(n))
  }
}

