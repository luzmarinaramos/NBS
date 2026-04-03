#' Restriction matrix for equicorrelation structure
#'
#' Constructs the restriction matrix \eqn{H} associated with the equicorrelation
#' structure of a covariance matrix \eqn{\Sigma}. Under equicorrelation,
#' all diagonal elements of \eqn{\Sigma} are equal and all off-diagonal elements
#' are equal. These conditions can be expressed as linear constraints on
#' \eqn{\phi = \mathrm{vech}(\Sigma)}, namely \eqn{H \phi = 0}.
#'
#' The resulting matrix \eqn{H} has dimension \eqn{(q - 2) \times q}, where
#' \eqn{q = p(p+1)/2} is the number of unique elements in \eqn{\Sigma}.
#'
#' @param p Integer. Dimension of the covariance matrix \eqn{\Sigma} (number of variables).
#'
#' @return A numeric matrix of dimension \eqn{(q - 2) \times q}, where each row
#' represents a linear restriction enforcing equicorrelation. The columns correspond
#' to the elements of \eqn{\phi = \mathrm{vech}(\Sigma)}.
#'
#' @details
#' The restriction matrix is constructed as follows:
#' \itemize{
#'   \item Diagonal constraints: \eqn{\sigma_{ii} - \sigma_{11} = 0}, for \eqn{i = 2, \dots, p}.
#'   \item Off-diagonal constraints: \eqn{\sigma_{ij} - \sigma_{21} = 0}, for all
#'   \eqn{i > j}, excluding the reference pair \eqn{(i,j) = (2,1)}.
#' }
#'
#' The function assumes that the operator \code{vech()} stacks the lower triangular
#' part of \eqn{\Sigma} column-wise.
#'
#' @examples
#' # Example for p = 3
#' H <- restriction_matrix_equicorrelation(3)
#' H
#'
#' @export
restriction_matrix_equicorrelation <- function(p) {
  q <- p * (p + 1) / 2

  # index matrix for vech(Sigma)
  idx <- matrix(NA_integer_, p, p)
  k <- 1
  for (j in 1:p) {
    for (i in j:p) {
      idx[i, j] <- k
      idx[j, i] <- k
      k <- k + 1
    }
  }

  rows <- list()

  # 1) Equal diagonal elements: sigma_ii - sigma_11 = 0, i=2,...,p
  base_diag <- idx[1, 1]
  for (i in 2:p) {
    r <- rep(0, q)
    r[idx[i, i]] <- 1
    r[base_diag] <- -1
    rows[[length(rows) + 1]] <- r
  }

  # 2) Equal off-diagonal elements: sigma_ij - sigma_21 = 0
  base_off <- idx[2, 1]
  for (j in 1:(p - 1)) {
    for (i in (j + 1):p) {
      if (!(i == 2 && j == 1)) {
        r <- rep(0, q)
        r[idx[i, j]] <- 1
        r[base_off] <- -1
        rows[[length(rows) + 1]] <- r
      }
    }
  }

  H <- do.call(rbind, rows)
  rownames(H) <- paste0("r", seq_len(nrow(H)))
  colnames(H) <- paste0("phi", seq_len(q))
  H
}
