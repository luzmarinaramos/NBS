#' @title Fisher Information Matrix for Multivariate Normal Birnbaum-Saunders Distribution
#' @description Computes the expected Fisher Information Matrix for the multivariate normal Birnbaum-Saunders distribution.
#' Some expectations are approximated using Monte Carlo integration.
#' @param mu Mean vector of length \code{p}.
#' @param sigma Covariance matrix of dimension \code{p} x \code{p}, must be symmetric and positive definite.
#' @param eta Positive scalar shape parameter.
#' @return A numeric matrix of dimension \code{p + p*(p+1)/2 + 1} by \code{p + p*(p+1)/2 + 1}, with attribute \code{"valid"} indicating if the matrix is symmetric and positive definite.
#' @examples
#' set.seed(1)
#' p <- 2
#' mu <- rep(0,p)
#' sigma <- diag(p)
#' eta <- 1
#' fisher_information(mu, sigma, eta)
#'
#' @importFrom matrixcalc duplication.matrix
#' @importFrom matrixcalc is.positive.definite
#' @export
fisher_information <- function(mu, sigma, eta) {
  # Checks
  if (!isSymmetric(sigma, tol = 1e-10)) stop("sigma must be symmetric.")
  if (!is.positive.definite(sigma)) stop("sigma must be positive definite.")
  if (!is.numeric(eta) || eta <= 0) stop("eta must be a positive scalar.")

  # Preparations
  p <- length(mu)
  q <- p * (p + 1) / 2
  a <- (p - 1) / 2
  sigma_inv <- solve(sigma)
  vec_sigma_inv <- as.vector(sigma_inv)
  D_p <- duplication.matrix(p)
  I_p2 <- diag(p^2)
  K_p <- commutation_matrix(p)
  N_p <- 0.5 * (I_p2 + K_p)

  # Derivatives of Bessel functions
  dbesselK <- function(u,d){ -0.5 * (besselK(u,d+1) + besselK(u,d-1)) }
  ddbesselK <- function(u,d){ 0.25 * (besselK(u,d+2) + 2*besselK(u,d) + besselK(u,d-2)) }

  # Functions used in expectations
  G <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    (eta/w)^a * besselK(w,a) + (eta/w)^(p-a) * besselK(w,p-a)
  }
  G_delta <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    w_delta <- eta/(2*w)
    t1 <- (eta/w)^a * (dbesselK(w,a) - a*besselK(w,a)/w)
    t2 <- (eta/w)^(p-a) * (dbesselK(w,p-a) - (p-a)*besselK(w,p-a)/w)
    (t1 + t2) * w_delta
  }
  G_eta <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    w_eta <- (2*eta+delta)/(2*w)
    t1 <- (eta/w)^a * (w_eta*dbesselK(w,a) + a*(1/eta - w_eta/w)*besselK(w,a))
    t2 <- (eta/w)^(p-a) * (w_eta*dbesselK(w,p-a) + (p-a)*(1/eta - w_eta/w)*besselK(w,p-a))
    t1 + t2
  }
  W_delta_eta <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    w_eta <- (2*eta+delta)/(2*w)
    w_delta <- eta/(2*w)
    gama1 <- function(x){(eta/w)^x * (x/eta - x*w_eta/w)}
    gama2 <- function(x){w_eta*ddbesselK(w,x) - (x*w_eta/w)*(dbesselK(w,x) - besselK(w,x)/w)}
    g <- G(eta,delta)/besselK(eta,0.5)
    g_delta <- G_delta(eta,delta)/besselK(eta,0.5)
    g_eta <- G_eta(eta,delta)/besselK(eta,0.5) - dbesselK(eta,0.5)*G(eta,delta)/besselK(eta,0.5)^2
    G_delta_eta <- ((eta/w)^a * gama2(a) + gama1(a)*(dbesselK(w,a) - a*besselK(w,a)/w) +
                      (eta/w)^(p-a) * gama2(p-a) + gama1(p-a)*(dbesselK(w,p-a) - (p-a)*besselK(w,p-a)/w)) * w_delta +
      G_delta(eta,delta)*(w - eta*w_eta)/(2*w^2*w_delta)
    g_delta_eta <- G_delta_eta/besselK(eta,0.5) - G_delta(eta,delta)*dbesselK(eta,0.5)/besselK(eta,0.5)^2
    (g_delta_eta - g_delta * g_eta / g) / g
  }
  W_eta_eta <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    w_eta <- (2*eta+delta)/(2*w)
    w_eta_eta <- (1 - w_eta^2)/w
    gama1 <- function(x){(eta/w)^x * (x/eta - x*w_eta/w)}
    gama3 <- function(x){x * (2*(w_eta/w)^2 - 1/eta^2 - 1/w^2)}
    gama5 <- function(x){w_eta_eta + x*w_eta*(1/eta - w_eta/w)}
    gama4 <- function(x){w_eta^2*ddbesselK(w,x) + gama5(x)*dbesselK(w,x) + gama3(x)*besselK(w,x)}
    g <- G(eta,delta)/besselK(eta,0.5)
    g_eta <- G_eta(eta,delta)/besselK(eta,0.5) - dbesselK(eta,0.5)*G(eta,delta)/besselK(eta,0.5)^2
    G_eta_eta <- (eta/w)^a * gama4(a) + gama1(a)*(w_eta*dbesselK(w,a) + a*(1/eta - w_eta/w)*besselK(w,a)) +
      (eta/w)^(p-a) * gama4(p-a) + gama1(p-a)*(w_eta*dbesselK(w,p-a) + (p-a)*(1/eta - w_eta/w)*besselK(w,p-a))
    g_eta_eta <- (1/besselK(eta,0.5)) * (G_eta_eta - 2*G_eta(eta,delta)*dbesselK(eta,0.5)/besselK(eta,0.5) -
                                           G(eta,delta)*(ddbesselK(eta,0.5)/besselK(eta,0.5) - 2*dbesselK(eta,0.5)^2/besselK(eta,0.5)^2))
    (g_eta_eta - g_eta^2 / g) / g
  }

  # Monte Carlo expectations
  N <- 10000
  z <- rNBS(N, mu, sigma, eta)  # z: N x p matrix

  stats <- apply(z, 1, function(zi) {
    zi <- as.vector(zi)
    delta <- as.numeric(t(zi - mu) %*% sigma_inv %*% (zi - mu))
    W_delta <- G_delta(eta, delta) / G(eta, delta)
    c(W_delta^2 * delta,
      W_delta^2 * delta^2,
      W_delta_eta(eta, delta) * delta,
      W_eta_eta(eta, delta))
  })

  means <- rowMeans(stats)
  d_g <- means[1]; f_g <- means[2]
  E_W_delta_eta_U <- means[3]; E_W_eta_eta <- means[4]
  c_g <- f_g / (p * (p + 2))

  # Fisher information matrix
  I_mu_mu   <- (4 * d_g / p) * sigma_inv
  ceros_p_q <- matrix(0, p, q)
  ceros_p_1 <- matrix(0, p, 1)
  ceros_q_p <- matrix(0, q, p)
  ceros_1_p <- matrix(0, 1, p)

  I_phi_eta <- (1 / p) * E_W_delta_eta_U * t(D_p) %*% vec_sigma_inv
  I_eta_phi <- t(I_phi_eta)
  I_eta_eta <- -E_W_eta_eta

  kron_part  <- 8 * c_g * kronecker(sigma_inv, sigma_inv) %*% N_p
  outer_part <- (4 * c_g - 1) * (vec_sigma_inv %*% t(vec_sigma_inv))
  I_phi_phi  <- 0.25 * t(D_p) %*% (kron_part + outer_part) %*% D_p

  top_row    <- cbind(I_mu_mu, ceros_p_q, ceros_p_1)
  mid_row    <- cbind(ceros_q_p, I_phi_phi, I_phi_eta)
  bottom_row <- cbind(ceros_1_p, I_eta_phi, c(I_eta_eta))
  Information <- rbind(top_row, mid_row, bottom_row)

  return(Information)
}
