orthogonalize <- function(X, group)
{
  n <- nrow(X)
  J <- max(group)
  Uinv <- vector("list", J)
  XX <- matrix(NA, nrow=n, ncol=ncol(X))
  for (j in 1:J) {
    ind <- which(group==j)
    U <- chol(crossprod(X[,ind])/n)
    Uinv[[j]] <- backsolve(U,diag(length(ind)))
    XX[,ind] <- X[,ind] %*% Uinv[[j]]
  }
  attr(XX,"Uinv") <- Uinv
  XX
}
unorthogonalize <- function(b, XX, group)
{
  beta <- b
  J <- max(group)
  for (j in 1:J) {
    ind <- which(group==j)
    beta[ind+1,] <- attr(XX,"Uinv")[[j]] %*% b[ind+1,]
  }
  beta
}