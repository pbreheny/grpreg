orthogonalize <- function(X, group)
{
  n <- nrow(X)
  J <- max(group)
  T <- vector("list", J)
  XX <- matrix(0, nrow=nrow(X), ncol=ncol(X))
  ##Uinv <- vector("list", J)
  XX[,which(group==0)] <- X[,which(group==0)]
  for (j in 1:J) {
    ind <- which(group==j)
    if (length(ind)==0) next
    SVD <- svd(crossprod(X[, ind, drop=FALSE])/n)
    r <- which(SVD$d > 1e-10)
    T[[j]] <- sweep(SVD$u[,r,drop=FALSE], 2, SVD$d[r]^(-1/2), "*")
    XX[,ind[r]] <- X[,ind]%*%T[[j]]
    ##result <- try(U <- chol(crossprod(X[, ind, drop=FALSE])/n))
    ##if (class(result)=="try-error") stop(paste("Group ", j," is not full-rank"))
    ##Uinv[[j]] <- backsolve(U,diag(length(ind)))
    ##XX[,ind] <- X[,ind] %*% Uinv[[j]]
  }
  nz <- !apply(XX==0,2,all)
  XX <- XX[, nz, drop=FALSE]
  attr(XX, "T") <- T
  attr(XX, "group") <- group[nz]
  ##attr(XX,"Uinv") <- Uinv
  XX
}
unorthogonalize <- function(b, XX, group.full, group)
{
  beta <- matrix(NA, nrow=length(group.full)+1, ncol=ncol(b))
  beta[c(1,1+which(group==0)),] <- b[c(1,1+which(group==0)),]
  J <- max(group)
  for (j in 1:J) {
    ind.full <- which(group.full==j)
    ind <- which(group==j)
    if (length(ind)==0) next
    beta[ind.full+1,] <- attr(XX,"T")[[j]] %*% b[ind+1,]
  }
  beta
}
