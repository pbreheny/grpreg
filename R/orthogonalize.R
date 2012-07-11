orthogonalize <- function(X, group)
{
  n <- nrow(X)
  p <- ncol(X)
  J <- group[p]
  Uinv <- vector("list", J)
  XX <- matrix(NA, nrow=n, ncol=p)
  for (j in 1:J) {
    ind <- which(group==j)
    U <- chol(crossprod(X[,ind])/n)
    Uinv[[j]] <- backsolve(U,diag(length(ind)))
    XX[,ind] <- X[,ind] %*% Uinv[[j]]
  }
  attr(XX,"Uinv") <- Uinv
  XX
}
