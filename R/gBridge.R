gBridge <- function(X, y, group=1:ncol(X), family=c("gaussian","binomial"), nlambda=100, lambda, lambda.min={if (nrow(X) > ncol(X)) .001 else .05}, lambda.max, alpha=1, eps=.005, delta=1e-7, max.iter=1000, gamma=0.5, group.multiplier=rep(1,J), warn=TRUE)
{
  ## Check for errors
  family <- match.arg(family)
  if (alpha > 1 | alpha < 0) stop("alpha must be in [0,1]")
  if (length(group)!=ncol(X)) stop("group does not match X")
  if (delta <= 0) stop("Delta must be a positive number")
  J <- max(group)
  K <- as.numeric(table(group))
  if (!(identical(as.integer(sort(unique(group))),as.integer(1:J)) | identical(as.integer(sort(unique(group))),as.integer(0:J)))) stop("Groups must be consecutively numbered 1,2,3,...")
  if (length(group.multiplier)!=J) stop("Length of group.multiplier must equal number of penalized groups")

  ## Set up X, y, lambda
  XX <- standardize(X)
  center <- attr(XX, "center")
  scale <- attr(XX, "scale")
  yy <- if (family=="gaussian") y - mean(y) else y
  if (missing(lambda)) {
    lambda <- setupLambda.gBridge(XX, yy, group, family, alpha, lambda.min, lambda.max, nlambda, gamma, group.multiplier)
  } else {
    nlambda <- length(lambda)
  }
  
  ## Fit
  n <- length(yy)
  p <- ncol(XX)
  if (family=="gaussian") {
    fit <- .C("gpPathFit_gaussian", double(p*nlambda), integer(nlambda), double(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(group), as.integer(n), as.integer(p), "gBridge", as.integer(J), as.integer(K), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.double(delta), as.integer(max.iter), as.double(gamma), as.double(0), as.integer(p), as.integer(group.multiplier), as.integer(TRUE))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1 ## Intercept
    loss <- fit[[4]]
  }
  if (family=="binomial") {
    fit <- .C("gpPathFit_binomial", double(nlambda), double(p*nlambda), integer(nlambda), double(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(group), as.integer(n), as.integer(p), "gBridge", as.integer(J), as.integer(K), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.double(delta), as.integer(max.iter), as.double(gamma), as.double(0), as.double(group.multiplier), as.integer(p), as.integer(warn), as.integer(TRUE))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    iter <- fit[[3]]
    df <- fit[[4]]
    loss <- fit[[5]]
  }

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(b[p,])
  b <- b[,ind,drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda")
  beta <- unstandardize(b, center, scale)
  
  
  ## Names
  if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
  else varnames <- colnames(X)
  varnames <- c("(Intercept)",varnames)
  dimnames(beta) <- list(varnames, round(lambda,digits=4))
  
  structure(list(beta = beta,
                 family = family,
                 group = group,
                 lambda = lambda,
                 alpha = alpha,
                 loss = loss,
                 n = length(y),
                 penalty = "gBridge",
                 df = df,
                 iter = iter),
            class = "grpreg")
}
