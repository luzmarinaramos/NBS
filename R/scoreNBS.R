library(matrixcalc)
scoreNBS <- function(y,mu,sigma,eta){
  n <- nrow(y)
  p <- length(mu)
  q <- p*(p + 1) / 2
  a <- (p-1)/2
  sigma_inv <- solve(sigma)
  D_p <- duplication.matrix(p)

  # bessel function derivatives
  dbesselK <- function(u,d){-0.5*(besselK(u,d+1)+besselK(u,d-1))}
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


  # contructing components of the score function
  U_mu <- rep(0,p)
  U_phi <- rep(0,q)
  U_eta <- 0
  for(i in 1:n){
    yi <- as.vector(y[i,])
    delta <- as.numeric(t(yi-mu)%*%sigma_inv%*%(yi-mu))
    W_delta <- as.numeric(G_delta(eta,delta)/ G(eta,delta))
    W_eta <- G_eta(eta,delta)/G(eta,delta)-dbesselK(eta,0.5)/besselK(eta,0.5)
    U_mu <- U_mu - (2* W_delta)*(sigma_inv%*%(yi-mu))
    M <- (2*W_delta)* sigma_inv%*%(yi-mu)%*%t(yi-mu)%*%sigma_inv + sigma_inv
    U_phi <- U_phi - 0.5*t(D_p)%*%vec(M)
    U_eta <- U_eta + W_eta
  }
  U <- c(U_mu,U_phi,U_eta)
  out <- list(score=U,U_mu=U_mu,U_phi=U_phi,U_eta=U_eta)
  return(out)
}
