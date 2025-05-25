#' Generate the commutation matrix of size p^2 x p^2
#'
#' Constructs the commutation matrix \(K_p\) which is a \(p^2 \times p^2\) permutation matrix
#' satisfying \(\text{vec}(A^T) = K_p \text{vec}(A)\) for any \(p \times p\) matrix \(A\).
#'
#' @param p Positive integer specifying the dimension of the square matrix.
#'
#' @return A \(p^2 \times p^2\) commutation matrix.
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

