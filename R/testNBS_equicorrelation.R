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
    I <- FI$Information
    I_phi_phi <- FI$I_phi_phi
    I_phi_eta <- FI$I_phi_eta
    I_eta_eta <- FI$I_eta_eta
    FI_H0 <- fisher_information(mu_est_H0,sigma_est_H0,eta_est_H0)
    I_H0 <- FI_H0$Information
    I_phi_phi_H0 <- FI_H0$I_phi_phi
    I_phi_eta_H0 <- FI_H0$I_phi_eta
    I_eta_eta_H0 <- FI_H0$I_eta_eta
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
    S1 <- (1/n)*t(U_H0)%*%solve(I_H0)%*%U_H0
    S1 <- as.numeric(S1)
    G <- t(U_phi_H0)%*%(phi_est - phi_est_H0)
    G <- as.numeric(G)
    G1 <- t(U_H0)%*%(theta_est - theta_est_H0)
    G1 <- as.numeric(G1)
  # 2. Calcular valor-p (p-value)
  valor_p_LR <- 1-pchisq(LR,gl)
  valor_p_W <- 1-pchisq(W,gl)
  valor_p_S <- 1-pchisq(S,gl)
  valor_p_S1 <- 1-pchisq(S1,gl)
  valor_p_G <- 1-pchisq(G,gl)
  valor_p_G1 <- 1-pchisq(G1,gl)

  #orden <- c('Likelihood Ratio','Wald','Score','Gradient')
  orden <- c('Likelihood Ratio','Wald','Score','Score corregido','Gradient','Gradient corregido')
  #estadistico <- c(LR,W,S,G)
  estadistico <- c(LR,W,S,S1,G,G1)
  #p_value <- c(valor_p_LR,valor_p_W,valor_p_S,valor_p_G)
  p_value <- c(valor_p_LR,valor_p_W,valor_p_S,valor_p_S1,valor_p_G,valor_p_G1)

  # 4. Devolver resultados
  return(list(order=orden,statistics=estadistico,p_value=p_value))
}
