#' Generate nonlinear example data
#' 
#' Mainly intended to demonstrate the use of basis expansion models for sparse additive modeling; intended for use with [expand_spline()].
#' 
#' @param n      Sample size (numeric; default = 100).
#' @param p      Number of features (numeric; default = 16).
#' @param seed   Change to get different random data sets (numeric; default = 1).
#' 
#' @examples
#' Data <- gen_nonlinear_data()

gen_nonlinear_data <- function(n=100, p=16, seed=1) {
  if (!(is.numeric(p) && p >= 6)) stop('p must be at least 6', call.=FALSE)
  X <- matrix(runif(n*p), nrow=n, ncol=p)
  w <- floor(log10(p)) + 1
  colnames(X) <- sprintf(paste0('V%0', w, 'd'), 1:p)
  f <- list(
    function(x){2*(exp(-10*x)-exp(-10))/(1-exp(-10)) - 1},
    function(x){-2*(exp(-10*x)-exp(-10))/(1-exp(-10)) + 1},
    function(x){2*x-1},
    function(x){-2*x+1},
    function(x){8*(x-0.5)^2 - 1},
    function(x){-8*(x-0.5)^2 + 1})
  eta <- matrix(NA, nrow=n, ncol=6)
  for (j in 1:6) eta[,j] <- f[[j]](X[,j])
  mu <- 5 + apply(eta, 1, sum)
  y <- rnorm(n, mean=mu, sd=0.25)
  list(X=X, y=y, mu=mu)
}
