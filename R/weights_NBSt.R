#### calculo de los pesos #####

pesos <- function(deltai, p, alpha, nu){

  d_BS <- function(u){
    ((u + 1) / (2 * alpha * u^(3/2))) * dnorm((sqrt(u) - 1/sqrt(u)) / alpha)
  }

  gi <- function(u){d_BS(u) * u^(nu/2) * (nu*u + deltai)^(-(nu + p)/2)}

  integrando1 <- function(u){((nu+p) / (nu*u + deltai)) * gi(u)  }
  integrando2 <- function(u){ u * gi(u)   }
  integrando3 <- function(u){ u^(-1) * gi(u)   }
  integrando4 <- function(u){ log(nu+deltai/u) * gi(u)   }
  integrando5 <- function(u){(((nu+p)*u) / (nu*u + deltai)) * gi(u)  }


  numerator1  <- integrate(integrando1, lower = 0, upper = Inf)$value
  numerator2  <- integrate(integrando2, lower = 0, upper = Inf)$value
  numerator3  <- integrate(integrando3, lower = 0, upper = Inf)$value
  numerator4  <- integrate(integrando4, lower = 0, upper = Inf)$value
  numerator5  <- integrate(integrando5, lower = 0, upper = Inf)$value
  denominator <- integrate(gi,         lower = 0, upper = Inf)$value

  w1 <-  numerator1 / denominator
  w2 <-  numerator2 / denominator
  w3 <-  numerator3 / denominator
  w4 <-  digamma((nu+p)/2)+log(2) - numerator4 / denominator
  w5 <-  numerator5 / denominator
  return(c(w1,w2,w3,w4,w5))
}


