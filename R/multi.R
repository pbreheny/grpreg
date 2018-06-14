multiX <- function(X, m) {
  p <- ncol(X)
  n <- nrow(X)
  A <- matrix(0, m*n, m*p)
  for (i in 1:m) {
    A[m*(1:n)-i+1, m*(1:p)-i+1] <- X
  }
  cbind(matrix(as.numeric(diag(m)),m*n,m,byrow=TRUE)[,2:m],A)
}
multiG <- function(g, ncolY) {
  structure(c(rep(0, ncolY-1), rep(g, each=ncolY)),
            levels=attr(g, 'levels'),
            m=attr(g, 'm'))
}
