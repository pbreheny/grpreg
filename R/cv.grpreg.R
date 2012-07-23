cv.grpreg <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gMCP", "gLasso"), family=c("gaussian","binomial"), nlambda=100, lambda, lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, alpha=1, nfolds=10, seed, trace=FALSE, group.multiplier=rep(1,J),...)
{
  ## Check for errors
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (alpha > 1 | alpha < 0) stop("alpha must be in [0,1]")
  if (length(group)!=ncol(X)) stop("group does not match X")
  if (is.null(colnames(X))) colnames(X) <- paste("V",1:ncol(X),sep="")
  J <- max(group)
  K <- as.numeric(table(group))
  if (!(identical(as.integer(sort(unique(group))),as.integer(1:J)) | identical(as.integer(sort(unique(group))),as.integer(0:J)))) stop("Groups must be consecutively numbered 1,2,3,...")
  if (length(group.multiplier)!=J) stop("Length of group.multiplier must equal number of penalized groups")
  if (!missing(seed)) set.seed(seed)
  
  ## Set up XX, yy, lambda
  XX <- standardize(X)
  center <- attr(XX, "center")
  scale <- attr(XX, "scale")
  if (strtrim(penalty,2)=="gr") XX <- orthogonalize(XX, group)
  yy <- if (family=="gaussian") y - mean(y) else y
  if (missing(lambda)) lambda <- setupLambda(XX, yy, group, family, penalty, alpha, lambda.min, nlambda, group.multiplier)
  rm(XX)
  
  error <- array(NA,dim=c(nfolds,length(lambda)))
  
  n <- length(y)
  if (family=="gaussian") {
    cv.ind <- ceiling(sample(1:n)/n*nfolds)
  } else if (family=="binomial") {
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

    fit.i <- grpreg(X1, y1, group=group, penalty=penalty, family=family, alpha=alpha, lambda=lambda, warn=FALSE, ...)
    yhat <- predict(fit.i, X2, type="response")
    error[i, 1:ncol(yhat)] <- loss.grpreg(y2, yhat, family)
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(error),2,all))
  E <- error[,ind]
  lambda <- lambda[ind]
  
  val <- list(E=E, cve=apply(E,2,sum)/n, lambda=lambda)
  val$min <- which.min(val$cve)
  val$lambda.min <- lambda[val$min]
  class(val) <- "cv.grpreg"
  val
}
