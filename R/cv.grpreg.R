cv.grpreg <- function(X, y, ..., nfolds=10, seed, trace=FALSE)
{
  if (!missing(seed)) set.seed(seed)
  fit <- grpreg(X=X, y=y, ...)
  E <- matrix(NA, nrow=length(y), ncol=length(fit$lambda))
  
  n <- length(y)
  if (fit$family=="gaussian") {
    cv.ind <- ceiling(sample(1:n)/n*nfolds)
  } else if (fit$family=="binomial") {
    ind1 <- which(y==1)
    ind0 <- which(y==0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
    cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
    cv.ind <- numeric(n)
    cv.ind[y==1] <- cv.ind1
    cv.ind[y==0] <- cv.ind0
  }
  
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #",i,sep="","\n")
    X1 <- X[cv.ind!=i,]
    y1 <- y[cv.ind!=i]
    X2 <- X[cv.ind==i,]
    y2 <- y[cv.ind==i]

    fit.i <- grpreg(X1, y1, warn=FALSE, ...)
    yhat <- predict(fit.i, X2, type="response")
    E[cv.ind==i, 1:ncol(yhat)] <- loss.grpreg(y2, yhat, fit$family)
  }

  ## Remove saturated values
  ind <- which(!apply(is.na(E),2,any))
  E <- E[,ind]
  lambda <- fit$lambda[ind]
  
  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)
  structure(list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min]), class="cv.grpreg")
}
