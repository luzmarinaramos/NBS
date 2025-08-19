#' Test for Equicorrelation Structure in the Multivariate Normal Distribution
#'
#' Performs hypothesis testing for the covariance matrix of a multivariate
#' normal distribution under the null hypothesis of equicorrelation structure.
#'
#' @param y A numeric matrix of size \eqn{n \times p}, where each row is an
#'   observation from a \eqn{p}-dimensional multivariate normal distribution.
#' @param nivel_significancia Numeric value in \eqn{(0,1)}, significance level
#'   used for the test (default is \code{0.05}).
#'
#' @details
#' The null hypothesis is that the covariance matrix has an
#' equicorrelation structure, i.e.
#' \deqn{\Sigma = \sigma^2 \left[ (1-\rho)I_p + \rho J_p \right],}
#' where \eqn{I_p} is the identity matrix and \eqn{J_p} is a \eqn{p \times p}
#' matrix of ones.
#'
#' The function computes four classical test statistics:
#' \itemize{
#'   \item \strong{Likelihood Ratio (LR)} test
#'   \item \strong{Wald} test
#'   \item \strong{Score (Rao)} test
#'   \item \strong{Gradient (Terrell)} test
#' }
#'
#' Each statistic is compared with the \eqn{\chi^2} distribution with
#' \eqn{df = p(p+1)/2 - 2} degrees of freedom.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{Test}{Name of the test (LR, Wald, Score, Gradient).}
#'   \item{Statistic}{Value of the test statistic.}
#'   \item{p_value}{Associated p-value under the \eqn{\chi^2(df)} distribution.}
#' }
#'
#' @examples
#' set.seed(123)
#' y <- matrix(rnorm(50), ncol = 5)
#' testNMV_equicorrelation(y)
#'
#' @references
#' Silvey, S.D. (1975). \emph{Statistical Inference}. Chapman and Hall.
#'
#' Terrell, G.R. (2002). The gradient statistic. \emph{Computing Science and
#' Statistics}, 34, 206-215.
#'
#' @export
testNMV_equicorrelation <- function(y, nivel_significancia = 0.05) {
  n <- nrow(y)
  p <- ncol(y)
  q <- p * (p + 1) / 2
  df <- q - 2  # degrees of freedom under H0

  # 1. Unrestricted fit (alternative hypothesis)
  fit <- emNMV(y)
  mu_est <- fit[[1]]
  sigma_est <- fit[[2]]
  phi_est <- vech(sigma_est)

  # 2. Restricted fit (null hypothesis: equicorrelation)
  fit_H0 <- emNMV_equicorrelation(y)
  mu_est_H0 <- fit_H0[[1]]
  sigma_est_H0 <- fit_H0[[2]]
  phi_est_H0 <- vech(sigma_est_H0)

  # 3. Fisher information under H1
  FI <- fisher_information_NMV(mu_est, sigma_est)
  FI <- (FI + t(FI)) / 2
  if (!is.positive.definite(FI)) stop("Fisher information matrix under H1 is not positive definite.")
  I_phi_phi <- FI[(p + 1):(p + q), (p + 1):(p + q)]

  # 4. Fisher information under H0
  FI_H0 <- fisher_information_NMV(mu_est_H0, sigma_est_H0)
  FI_H0 <- (FI_H0 + t(FI_H0)) / 2
  if (!is.positive.definite(FI_H0)) stop("Fisher information matrix under H0 is not positive definite.")
  I_phi_phi_H0 <- FI_H0[(p + 1):(p + q), (p + 1):(p + q)]

  # 5. Score vector under H0
  U_H0 <- score_function_NMV(y, mu_est_H0, sigma_est_H0)
  U_phi_H0 <- U_H0[(p + 1):(p + q)]

  # 6. Test statistics
  LR <- 2 * (loglikeNMV(y, mu_est, sigma_est) - loglikeNMV(y, mu_est_H0, sigma_est_H0))
  W  <- n * crossprod(phi_est - phi_est_H0, I_phi_phi %*% (phi_est - phi_est_H0))

  chol_I_phi_phi_H0 <- tryCatch(
    chol(I_phi_phi_H0),
    error = function(e) {
      stop("Error in Cholesky decomposition: I_phi_phi_H0 is not positive definite.")
    })
  inv_I_phi_phi_H0 <- chol2inv(chol_I_phi_phi_H0)
  S <- (1 / n) * crossprod(U_phi_H0, inv_I_phi_phi_H0 %*% U_phi_H0)

  G <- as.numeric(t(U_phi_H0) %*% (phi_est - phi_est_H0))

  # 7. p-values
  p_LR <- 1 - pchisq(LR, df)
  p_W  <- 1 - pchisq(W, df)
  p_S  <- 1 - pchisq(S, df)
  p_G  <- 1 - pchisq(G, df)

  # 8. Output
  return(data.frame(
    Test = c("Likelihood Ratio", "Wald", "Score", "Gradient"),
    Statistic = c(LR, W, S, G),
    p_value = c(p_LR, p_W, p_S, p_G)
  ))
}
