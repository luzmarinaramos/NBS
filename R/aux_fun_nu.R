update_nu <- function(nu, s5, n,
                      lower = 1e-6, upper = 1e6,
                      max_expand = 60, tol = 1e-10) {

  f <- function(x) log(x/2) + 1 - digamma(x/2) + (s5/n)

  # intervalo inicial alrededor de nu
  a <- max(lower, nu/2)
  b <- min(upper, nu*2)

  fa <- f(a); fb <- f(b)

  # expandir hasta lograr cambio de signo
  k <- 0
  while (k < max_expand && is.finite(fa) && is.finite(fb) && sign(fa) == sign(fb)) {
    a <- max(lower, a/2)
    b <- min(upper, b*2)
    fa <- f(a); fb <- f(b)
    k <- k + 1
  }

  # si no hay bracket, no hay raíz (o no la encontramos): devolver nu actual
  if (!is.finite(fa) || !is.finite(fb) || sign(fa) == sign(fb)) return(nu)

  uniroot(f, interval = c(a, b), tol = tol)$root
}
