testNBS_equicorrelation <- function(y, nivel_significancia = 0.05) {

    n <- nrow(y)
    p <- ncol(y)
    q <- p*(p+1)/2
    a <- (p-1)/2
    gl <- q-2
  # 1. Calcular estadístico de prueba
    fit <- emNBS(y)
    mu_est <- fit$mu_est
    sigma_est <- fit$sigma_est
    phi_est <- vech(sigma_est)
    eta_est <- fit$eta_est
    theta_est <- c(mu_est,phi_est,eta_est)
    fit_H0 <- emNBS_equicorrelation(y)
    mu_est_H0 <- fit_H0$mu_est
    sigma_est_H0 <- fit_H0$sigma_est
    phi_est_H0 <- vech(sigma_est_H0)
    eta_est_H0 <- fit_H0$eta_est
    theta_est_H0 <- c(mu_est_H0,phi_est_H0,eta_est_H0)

    FI <- fisher_information(mu_est,sigma_est,eta_est)
    if (!is.symmetric.matrix(FI)) {I <- (FI+t(FI))/2}
    if (!is.positive.definite(FI)) {stop("FI matrix is not positive definite.")}
    I_phi_phi <- FI[(p+1):(p+q),(p+1):(p+q)]
    I_phi_eta <- FI[(p+1):(p+q),(p+q+1)]
    I_eta_eta <- FI[(p+q+1),(p+q+1)]
    FI_H0 <- fisher_information(mu_est_H0,sigma_est_H0,eta_est_H0)
    if (!is.symmetric.matrix(FI_H0)) {I <- (FI_H0+t(FI_H0))/2}
    if (!is.positive.definite(FI_H0)) {stop("FI matrix is not positive definite.")}
    I_phi_phi_H0 <- FI_H0[(p+1):(p+q),(p+1):(p+q)]
    I_phi_eta_H0 <- FI_H0[(p+1):(p+q),(p+q+1)]
    I_eta_eta_H0 <- FI_H0[(p+q+1),(p+q+1)]
    FF_phi_phi <- I_phi_phi-I_phi_eta%*%t(I_phi_eta)/I_eta_eta
    FF_phi_phi_H0 <- I_phi_phi_H0-I_phi_eta_H0%*%t(I_phi_eta_H0)/I_eta_eta_H0
    Score_H0 <- scoreNBS(y,mu_est_H0,sigma_est_H0,eta_est_H0)
    U_H0 <- Score_H0$score
    U_phi_H0 <- Score_H0$U_phi

    LR <- 2*(fit$loglikelihood-fit_H0$loglikelihood)
    LR <- as.numeric(LR)
    W <- n*t(phi_est-phi_est_H0)%*%FF_phi_phi%*%(phi_est-phi_est_H0)
    W <- as.numeric(W)
    S <- (1/n)*t(U_phi_H0)%*%solve(FF_phi_phi_H0)%*%U_phi_H0
    S <- as.numeric(S)
    G <- t(U_phi_H0)%*%(phi_est - phi_est_H0)
    G <- as.numeric(G)
  # 2. Calcular valor-p (p-value)
  valor_p_LR <- 1-pchisq(LR,gl)
  valor_p_W <- 1-pchisq(W,gl)
  valor_p_S <- 1-pchisq(S,gl)
  valor_p_G <- 1-pchisq(G,gl)

  #orden <- c('Likelihood Ratio','Wald','Score','Gradient')
  orden <- c('Likelihood Ratio','Wald','Score','Gradient')
  #estadistico <- c(LR,W,S,G)
  estadistico <- c(LR,W,S,G)
  #p_value <- c(valor_p_LR,valor_p_W,valor_p_S,valor_p_G)
  p_value <- c(valor_p_LR,valor_p_W,valor_p_S,valor_p_G)

  # 4. Devolver resultados
  return(list(order=orden,statistics=estadistico,p_value=p_value))
}
