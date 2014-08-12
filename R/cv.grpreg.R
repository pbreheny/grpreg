cv.grpreg <- function(X, y, group=1:ncol(X), ..., nfolds=10, seed, trace=FALSE) {
  if (!missing(seed)) set.seed(seed)
  fit <- grpreg(X=X, y=y, group=group, ...)
  multi <- FALSE
  if (is.matrix(y) && ncol(y) > 1) {
    multi <- TRUE
    m <- ncol(y)
    p <- ncol(X)
    n <- nrow(X)
    y <- as.numeric(t(y))
    A <- matrix(0, m*n, m*p)
    for (i in 1:m) {
      A[m*(1:n)-2,m*(1:p)-2] <- X
      A[m*(1:n)-1,m*(1:p)-1] <- X
      A[m*(1:n),m*(1:p)] <- X
    }
    xnames <- colnames(X)
    X <- cbind(matrix(as.numeric(diag(m)),m*n,m,byrow=TRUE),A)
  }
  g <- if (multi) c(rep(0,m), fit$group) else fit$group
  E <- matrix(NA, nrow=length(y), ncol=length(fit$lambda))
  if (fit$family=="binomial") PE <- E
  
  n <- length(y)
  if (fit$family=="binomial") {
    ind1 <- which(y==1)
    ind0 <- which(y==0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
    cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
    cv.ind <- numeric(n)
    cv.ind[y==1] <- cv.ind1
    cv.ind[y==0] <- cv.ind0    
  } else {
    cv.ind <- ceiling(sample(1:n)/n*nfolds)    
  }

  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #",i,sep="","\n")
    X1 <- X[cv.ind!=i, , drop=FALSE]
    y1 <- y[cv.ind!=i]
    X2 <- X[cv.ind==i, , drop=FALSE]
    y2 <- y[cv.ind==i]

    fit.i <- grpreg(X1, y1, g, lambda=fit$lambda, warn=FALSE, ...)
    yhat <- predict(fit.i, X2, type="response")
    E[cv.ind==i, 1:length(fit.i$lambda)] <- loss.grpreg(y2, yhat, fit$family)
    if (fit$family=="binomial") PE[cv.ind==i, 1:length(fit.i$lambda)] <- (yhat < 0.5) == y2
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[,ind,drop=FALSE]
  lambda <- fit$lambda[ind]
  
  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)
  null.dev <- calcNullDev(X, y, group=g, family=fit$family)
  
  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=null.dev)
  if (fit$family=="binomial") val$pe <- apply(PE[,ind], 2, mean)
  structure(val, class="cv.grpreg")
}
