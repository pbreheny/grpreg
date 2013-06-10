multiX <- function(X, m) {
  p <- ncol(X)
  n <- nrow(X)
  A <- matrix(0, m*n, m*p)
  for (i in 1:m) {
    A[m*(1:n)-2,m*(1:p)-2] <- X
    A[m*(1:n)-1,m*(1:p)-1] <- X
    A[m*(1:n),m*(1:p)] <- X
  }
  X <- cbind(matrix(as.numeric(diag(m)),m*n,m,byrow=TRUE)[,2:m],A)
}
multiY <- function(y) {
  y <- as.numeric(t(y))
}
