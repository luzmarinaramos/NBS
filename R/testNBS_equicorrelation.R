#' @title Equicorrelation Test for the Normal-Birnbaum-Saunders Distribution
#' @description Performs four classical asymptotic hypothesis tests (Likelihood Ratio, Wald, Score, and Gradient)
#' to test whether the covariance matrix of a multivariate Normal-Birnbaum-Saunders distribution is equicorrelated.
#'
#' @param y A numeric matrix of dimension n x p containing the sample.
#' @param nivel_significancia Significance level (default is 0.05).
#'
#' @return A list with the following elements:
#' \item{order}{A character vector with the names of the tests.}
#' \item{statistics}{A numeric vector with the test statistics.}
#' \item{p_value}{A numeric vector with the corresponding p-values.}
#'
#' @details The function assumes that \code{emNBS} and \code{emNBS_equicorrelation} return objects
#' that include a logical element \code{converged}. The tests are performed only if both models converge
#' and the corresponding Fisher information matrices are positive definite. Symmetrization is applied
#' to the Fisher matrix if necessary. The inverse is computed via Cholesky decomposition for numerical stability.
#'
#' @importFrom matrixcalc is.positive.definite
#' @importFrom stats pchisq
#' @importFrom base chol chol2inv crossprod
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3)
#' result <- testNBS_equicorrelation(y)
#' print(result)
#' }
#'
#' @export
testNBS_equicorrelation <- function(y, nivel_significancia = 0.05) {
  n <- nrow(y)
  p <- ncol(y)
  q <- p * (p + 1) / 2
  df <- q - 2  # degrees of freedom

  # 1. Unrestricted fit (alternative hypothesis)
  fit <- emNBS(y)
  if (!fit$convergence) stop("EM algorithm (unrestricted model) did not converge.")

  mu_est <- fit$mu_est
  sigma_est <- fit$sigma_est
  phi_est <- vech(sigma_est)
  eta_est <- fit$eta_est

  # 2. Restricted fit (null hypothesis: equicorrelation)
  fit_H0 <- emNBS_equicorrelation(y)
  if (!fit_H0$convergence) stop("EM algorithm (equicorrelation model) did not converge.")

  mu_est_H0 <- fit_H0$mu_est
  sigma_est_H0 <- fit_H0$sigma_est
  phi_est_H0 <- vech(sigma_est_H0)
  eta_est_H0 <- fit_H0$eta_est

  # 3. Fisher information matrix under H1
  FI <- fisher_information(mu_est, sigma_est, eta_est)
  FI <- (FI + t(FI)) / 2  # symmetrize
  if (!is.positive.definite(FI)) stop("Fisher information matrix under H1 is not positive definite.")

  I_phi_phi <- FI[(p + 1):(p + q), (p + 1):(p + q)]
  I_phi_eta <- FI[(p + 1):(p + q), (p + q + 1)]
  I_eta_eta <- FI[(p + q + 1), (p + q + 1)]

  FF_phi_phi <- I_phi_phi - I_phi_eta %*% t(I_phi_eta) / I_eta_eta

  # 4. Fisher information matrix under H0
  FI_H0 <- fisher_information(mu_est_H0, sigma_est_H0, eta_est_H0)
  FI_H0 <- (FI_H0 + t(FI_H0)) / 2  # symmetrize
  if (!is.positive.definite(FI_H0)) stop("Fisher information matrix under H0 is not positive definite.")

  I_phi_phi_H0 <- FI_H0[(p + 1):(p + q), (p + 1):(p + q)]
  I_phi_eta_H0 <- FI_H0[(p + 1):(p + q), (p + q + 1)]
  I_eta_eta_H0 <- FI_H0[(p + q + 1), (p + q + 1)]

  FF_phi_phi_H0 <- I_phi_phi_H0 - I_phi_eta_H0 %*% t(I_phi_eta_H0) / I_eta_eta_H0

  # 5. Score vector under H0
  U_H0 <- scoreNBS(y, mu_est_H0, sigma_est_H0, eta_est_H0)
  U_phi_H0 <- U_H0[(p + 1):(p + q)]

  # 6. Test statistics
  LR <- 2 * (fit$loglikelihood - fit_H0$loglikelihood)
  W <- n * crossprod(phi_est - phi_est_H0, FF_phi_phi %*% (phi_est - phi_est_H0))

  chol_FF_H0 <- chol(FF_phi_phi_H0)
  inv_FF_H0 <- chol2inv(chol_FF_H0)
  S <- (1 / n) * crossprod(U_phi_H0, inv_FF_H0 %*% U_phi_H0)

  G <- as.numeric(t(U_phi_H0) %*% (phi_est - phi_est_H0))

  # 7. p-values
  p_LR <- 1 - pchisq(LR, df)
  p_W  <- 1 - pchisq(W, df)
  p_S  <- 1 - pchisq(S, df)
  p_G  <- 1 - pchisq(G, df)

  # 8. Output
  return(list(
    order = c("Likelihood Ratio", "Wald", "Score", "Gradient"),
    statistics = c(LR, W, S, G),
    p_value = c(p_LR, p_W, p_S, p_G)
  ))
}
