library(matrixcalc)
fisher_information <- function(mu,sigma,eta){
  p <- length(mu)
  q <- p*(p + 1) / 2
  a <- (p-1)/2
  sigma_inv <- solve(sigma)
  vec_sigma_inv <- as.vector(sigma_inv)
  D_p <- duplication.matrix(p)
  I_p2 <- diag(p^2)
  K_p <- commutation_matrix(p)
  N_p <- 0.5 * (I_p2 + K_p)

  # bessel function derivatives
  dbesselK <- function(u,d){-0.5*(besselK(u,d+1)+besselK(u,d-1))}
  ddbesselK <- function(u,d){0.25*(besselK(u,d+2)+2*besselK(u,d)+besselK(u,d-2))}
  # some functions
  G <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    out <- (eta/w)^(a)*besselK(w,a)+(eta/w)^(p-a)*besselK(w,p-a)
    return(out)
  }
  G_delta <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    w_delta <- eta/(2*w)
    t1 <- (eta/w)^(a)*(dbesselK(w,a)-a*besselK(w,a)/w)
    t2 <- (eta/w)^(p-a)*(dbesselK(w,p-a)-(p-a)*besselK(w,p-a)/w)
    out <- (t1+t2)*w_delta
    return(out)
  }
  G_eta <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    w_eta <- (2*eta+delta)/(2*w)
    t1 <- (eta/w)^(a)*(w_eta*dbesselK(w,a)+a*(1/eta-w_eta/w)*besselK(w,a))
    t2 <- (eta/w)^(p-a)*(w_eta*dbesselK(w,p-a)+(p-a)*(1/eta-w_eta/w)*besselK(w,p-a))
    out <- t1+t2
    return(out)
  }
  W_delta_eta <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    w_eta <- (2*eta+delta)/(2*w)
    w_delta <- eta/(2*w)
    gama1 <- function(x){(eta/w)^(x)*(x/eta-x*w_eta/w)}
    gama2 <- function(x){w_eta*ddbesselK(w,x)-(x*w_eta/w)*(dbesselK(w,x)-besselK(w,x)/w)}
    g <- G(eta,delta)/besselK(eta,0.5)
    g_delta <- G_delta(eta,delta)/besselK(eta,0.5)
    g_eta <- G_eta(eta,delta)/besselK(eta,0.5)-dbesselK(eta,0.5)*G(eta,delta)/besselK(eta,0.5)^2
    G_delta_eta <- ((eta/w)^(a)*gama2(a)+gama1(a)*(dbesselK(w,a)-a*besselK(w,a)/w)+(eta/w)^(p-a)*gama2(p-a)+gama1(p-a)*(dbesselK(w,p-a)-(p-a)*besselK(w,p-a)/w))*w_delta+G_delta(eta,delta)*(w-eta*w_eta)/(2*w^2*w_delta)
    g_delta_eta <- G_delta_eta/besselK(eta,0.5)-G_delta(eta,delta)*dbesselK(eta,0.5)/besselK(eta,0.5)^2
    out <- (g_delta_eta-g_delta*g_eta/g)/g
    return(out)
  }
  W_eta_eta <- function(eta,delta){
    w <- sqrt(eta*(eta+delta))
    w_eta <- (2*eta+delta)/(2*w)
    w_eta_eta <- (1-w_eta^2)/w
    gama1 <- function(x){(eta/w)^(x)*(x/eta-x*w_eta/w)}
    gama3 <- function(x){x*(2*(w_eta/w)^2-1/eta^2-1/w^2)}
    gama5 <- function(x){w_eta_eta+x*w_eta*(1/eta-w_eta/w)}
    gama4 <- function(x){w_eta^2*ddbesselK(w,x)+gama5(x)*dbesselK(w,x)+gama3(x)*besselK(w,x)}
    g <- G(eta,delta)/besselK(eta,0.5)
    g_eta <- G_eta(eta,delta)/besselK(eta,0.5)-dbesselK(eta,0.5)*G(eta,delta)/besselK(eta,0.5)^2
    G_eta_eta <- (eta/w)^(a)*gama4(a)+gama1(a)*(w_eta*dbesselK(w,a)+a*(1/eta-w_eta/w)*besselK(w,a))+ (eta/w)^(p-a)*gama4(p-a)+gama1(p-a)*(w_eta*dbesselK(w,p-a)+(p-a)*(1/eta-w_eta/w)*besselK(w,p-a))
    g_eta_eta <- (1/besselK(eta,0.5))*(G_eta_eta-2*G_eta(eta,delta)*dbesselK(eta,0.5)/besselK(eta,0.5)-G(eta,delta)*(ddbesselK(eta,0.5)/besselK(eta,0.5)-2*dbesselK(eta,0.5)^2/besselK(eta,0.5)^2))
    out <- (g_eta_eta-g_eta^2/g)/g
    return(out)
  }


  # some expectantias
  N <- 10000
  z <- rNBS(N,mu,sigma,eta)

  cumsum1 <- 0
  cumsum2 <- 0
  cumsum3 <- 0
  cumsum4 <- 0
  for(i in 1:N){
    zi <- as.vector(z[i,])
    delta <- as.numeric(t(zi-mu)%*%sigma_inv%*%(zi-mu))
  W_delta <- G_delta(eta,delta)/ G(eta,delta)
  cumsum1 <- cumsum1 + W_delta^2*delta
  cumsum2 <- cumsum2 + W_delta^2*delta^2
  cumsum3 <- cumsum3 + W_delta_eta(eta,delta)*delta
  cumsum4 <- cumsum4 + W_eta_eta(eta,delta)
  }
  d_g <- cumsum1 /N
  f_g <- cumsum2/N
  E_W_delta_eta_U <- cumsum3/N
  E_W_eta_eta <- cumsum4/N
  c_g <- f_g/(p*(p+2))

    # fisher information

  I_mu_mu   <- (4*d_g/p)* sigma_inv
  ceros_p_q <- matrix(0,nrow=p,ncol=q)
  ceros_p_1 <- matrix(0,nrow=p,ncol=1)
  ceros_q_p <- matrix(0,nrow=q,ncol=p)

  I_phi_eta <- (1 / p) * E_W_delta_eta_U * t(D_p) %*% vec_sigma_inv
  ceros_1_p <- matrix(0,nrow=1,ncol=p)
  I_eta_phi <- t(I_phi_eta)
  I_eta_eta <- - E_W_eta_eta

  kron_part <- 8 * c_g * kronecker(sigma_inv, sigma_inv) %*% N_p
  outer_part <- (4 * c_g - 1) * (vec_sigma_inv %*% t(vec_sigma_inv))
  I_phi_phi <- 0.25 * t(D_p) %*% (kron_part + outer_part) %*% D_p
  top_row <- cbind(I_mu_mu,ceros_p_q,ceros_p_1)
  mid_row <- cbind(ceros_q_p,I_phi_phi,I_phi_eta)
  bottom_row <- cbind(ceros_1_p,I_eta_phi,I_eta_eta)
  Information <- rbind(top_row,mid_row,bottom_row)
  return(Information)
}
