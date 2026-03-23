#' Density of the NBSt distribution (Normal–Birnbaum–Saunders–t mixture)
#'
#' Computes the density (or log-density) of the NBSt distribution defined through
#' the elliptical representation
#' \deqn{ f_Y(y)=C_{p,\nu}\,|\Sigma|^{-1/2}\,g(\delta), }
#' where \eqn{\delta=(y-\mu)^\top \Sigma^{-1}(y-\mu)} and the generator \eqn{g(\delta)}
#' is given by the one-dimensional integral
#' \deqn{
#' g(\delta)=
#' \frac{e^{1/\alpha^2}}{2\alpha\sqrt{2\pi}}
#' \int_0^\infty
#' \left(1+\frac{\delta}{\nu u}\right)^{-(p+\nu)/2}
#' \left(u^{-(p+1)/2}+u^{-(p+3)/2}\right)
#' \exp\!\left\{-\frac{1}{2\alpha^2}\left(u+u^{-1}\right)\right\}du.
#' }
#'
#' The integral is evaluated numerically using the stable change of variables
#' \eqn{u=\exp(s)} with a default truncation \eqn{s\in[-L,L]} (the integrand decays
#' super-exponentially due to the \eqn{\exp\{-\cosh(s)/\alpha^2\}} term).
#'
#' @param y A numeric vector of length \eqn{p}, or a numeric matrix/data.frame with
#'   \eqn{p} columns (each row is an observation).
#' @param mu Numeric vector of length \eqn{p} giving the location parameter.
#' @param Sigma Positive-definite \eqn{p\times p} scale matrix.
#' @param nu Degrees of freedom, a positive scalar.
#' @param alpha Birnbaum--Saunders shape parameter \eqn{\alpha>0} (with scale fixed at 1).
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @param L Positive scalar; truncation limit for the transformed integral in
#'   \eqn{s=\log(u)}. The integration is performed over \eqn{[-L,L]}.
#' @param rel.tol Relative tolerance passed to \code{\link[stats]{integrate}}.
#' @param abs.tol Absolute tolerance passed to \code{\link[stats]{integrate}}.
#' @param subdivisions Maximum number of subintervals used by
#'   \code{\link[stats]{integrate}}.
#'
#' @details
#' The normalizing constant \eqn{C_{p,\nu}} is
#' \deqn{
#' C_{p,\nu}=(2\pi)^{-p/2}\left(\frac{\nu}{2}\right)^{-p/2}
#' \frac{\Gamma\left(\frac{p+\nu}{2}\right)}{\Gamma\left(\frac{\nu}{2}\right)}.
#' }
#'
#' Internally, \eqn{\delta} is computed using a Cholesky factorization of \eqn{\Sigma}
#' for numerical stability. The integrand is evaluated on the log-scale where possible
#' to reduce overflow/underflow.
#'
#' If numerical integration fails for a given observation (rare; typically due to
#' extreme parameter values), the function returns \code{NA} for that observation.
#' Increasing \code{L} and/or \code{subdivisions} usually resolves it.
#'
#' @return A numeric scalar if \code{y} is a vector, or a numeric vector of length
#'   \code{nrow(y)} if \code{y} is a matrix/data.frame. Values are densities unless
#'   \code{log=TRUE}.
#'
#' @examples
#' set.seed(1)
#' p <- 3
#' y  <- rnorm(p)
#' mu <- rep(0, p)
#' Sigma <- diag(p)
#' nu <- 5
#' alpha <- 0.8
#'
#' dNBSt(y, mu, Sigma, nu, alpha)
#' dNBSt(y, mu, Sigma, nu, alpha, log = TRUE)
#'
#' # Multiple observations (rows)
#' Y <- matrix(rnorm(30), ncol = 3)
#' dNBSt(Y, mu, Sigma, nu, alpha)
#'
#' @seealso \code{\link[stats]{integrate}} for numerical integration.
#'
#' @export
dNBSt <- function(y, mu, Sigma, nu, alpha, log = FALSE,
                  L = 30, rel.tol = 1e-8, abs.tol = 0, subdivisions = 2000) {

  # --- stable helpers ---
  log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
  log1p_exp <- function(logx) ifelse(logx > 35, logx, log1p(exp(logx))) # log(1+exp(logx))

  # --- input handling ---
  if (is.data.frame(y)) y <- as.matrix(y)
  y_mat <- if (is.matrix(y)) y else matrix(as.numeric(y), nrow = 1)

  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)

  p <- ncol(y_mat)
  if (length(mu) != p) stop("mu must have length p (same dimension as y).")
  if (!all(dim(Sigma) == c(p, p))) stop("Sigma must be a p x p matrix.")
  if (!is.numeric(nu) || length(nu) != 1 || nu <= 0) stop("nu must be > 0.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) stop("alpha must be > 0.")

  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) stop("Sigma must be symmetric positive-definite (chol failed).")
  logdetSigma <- 2 * sum(log(diag(R)))

  # C_{p,nu} from your eq. (C_pnu)
  logC_pnu <- -(p/2) * log(2*pi) - (p/2) * log(nu/2) +
    lgamma((p + nu)/2) - lgamma(nu/2)

  # prefactor of g(delta) from your eq. (g_delta_integral_explicit)
  log_pref_g <- (1/alpha^2) - log(2*alpha) - 0.5*log(2*pi)

  # compute delta(y) = (y-mu)' Sigma^{-1} (y-mu)
  delta_one <- function(y_vec) {
    z <- as.numeric(y_vec) - mu
    v <- backsolve(R, z, transpose = TRUE)
    w <- backsolve(R, v)
    sum(z * w)
  }

  # g(delta) via u = exp(s): du = exp(s) ds
  g_of_delta <- function(delta) {
    # integrand in s
    integrand <- function(s) {
      # stable log(1 + delta/(nu*u)), u=exp(s) => delta/(nu*u) = (delta/nu)*exp(-s)
      if (delta == 0) {
        log_term <- 0
      } else {
        logx <- log(delta/nu) - s
        log_term <- log1p_exp(logx)
      }

      # (u^{-(p+1)/2} + u^{-(p+3)/2}) = u^{-(p+3)/2} (u+1)
      # log = -((p+3)/2)s + log(1+exp(s))
      log_upart <- -((p + 3)/2) * s + log1pexp(s)

      # exp{-(u + u^{-1})/(2 alpha^2)} with u=exp(s) => exp{ -cosh(s)/alpha^2 }
      # since (exp(s)+exp(-s))/2 = cosh(s)
      log_exp_part <- -cosh(s) / (alpha^2)

      # (1 + delta/(nu u))^{-(p+nu)/2}
      log_pow <- -((p + nu)/2) * log_term

      # Jacobian du = exp(s) ds => +s
      log_int <- log_pow + log_upart + log_exp_part + s

      out <- exp(log_int)
      out[!is.finite(out)] <- 0
      out
    }

    val <- tryCatch(
      integrate(integrand, lower = -L, upper = L,
                rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions),
      error = function(e) NULL
    )

    if (is.null(val) || !is.finite(val$value) || val$value < 0) return(NA_real_)
    exp(log_pref_g) * val$value
  }

  dens_one <- function(y_vec) {
    dlt <- delta_one(y_vec)
    g   <- g_of_delta(dlt)
    if (!is.finite(g) || g <= 0) return(NA_real_)

    logf <- logC_pnu - 0.5*logdetSigma + log(g)
    if (log) logf else exp(logf)
  }

  out <- apply(y_mat, 1, dens_one)
  out
}
