grpreg <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gMCP", "gBridge", "gLasso"), family=c("gaussian","binomial"), nlambda=100, lambda, lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, alpha=1, eps=.005, max.iter=1000, dfmax=p, gamma = 3, group.multiplier={if (strtrim(penalty,2)=="gr") sqrt(table(group[group!=0])) else rep(1,J)}, warn=TRUE, ...)
{
  ## Check for errors
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (penalty=="gLasso") penalty <- "grLasso"
  if (penalty=="gBridge") stop("gBridge has been divorced from the grpreg function; use the gBridge() function instead")  
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")
  if (length(group)!=ncol(X)) stop("group does not match X")
  J <- max(group)
  if (!(identical(as.integer(sort(unique(group))),as.integer(1:max(group))) | identical(as.integer(sort(unique(group))),as.integer(0:max(group))))) stop("Groups must be consecutively numbered 1,2,3,...")
  if (length(group.multiplier)!=max(group)) stop("Length of group.multiplier must equal number of penalized groups")
  
  ## Set up X, y, lambda
  XX <- standardize(X)
  center <- attr(XX, "center")
  scale <- attr(XX, "scale")
  nz <- which(scale > 1e-6)
  nzg <- setdiff(unique(group),unique(group[nz]))
  if (length(nzg)) {
    J  <- J - length(nzg)
    group.multiplier <- group.multiplier[-nzg]    
  }
  XX <- XX[ ,nz, drop=FALSE]
  group.orig <- group
  group <- group[nz]
  K <- as.numeric(table(group))
  if (strtrim(penalty,2)=="gr") XX <- orthogonalize(XX, group)
  yy <- if (family=="gaussian") y - mean(y) else y
  if (missing(lambda)) {
    lambda <- setupLambda(XX, yy, group, family, penalty, alpha, lambda.min, nlambda, group.multiplier)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  n <- length(yy)
  p <- ncol(XX)
  K0 <- if (min(group)==0) K[1] else 0
  K1 <- if (min(group)==0) cumsum(K) else c(0, cumsum(K))
  tau <- 1/3
  if (family=="gaussian") {
    if (strtrim(penalty,2)=="gr") fit <- .C("grPathFit_gaussian", double(p*nlambda), integer(nlambda), double(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(n), as.integer(p), penalty, as.integer(J), as.integer(K1), as.integer(K0), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.integer(max.iter), as.double(gamma), as.double(group.multiplier), as.integer(dfmax), as.integer(user.lambda))
    else fit <- .C("gpPathFit_gaussian", double(p*nlambda), integer(nlambda), double(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(n), as.integer(p), penalty, as.integer(J), as.integer(K1), as.integer(K0), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.double(0), as.integer(max.iter), as.double(gamma), as.double(tau), as.integer(dfmax), as.double(group.multiplier), as.integer(user.lambda))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1 ## Intercept
    loss <- fit[[4]]
  }
  if (family=="binomial") {
    if (strtrim(penalty,2)=="gr") fit <- .C("grPathFit_binomial", double(nlambda), double(p*nlambda), integer(nlambda), double(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(n), as.integer(p), penalty, as.integer(J), as.integer(K1), as.integer(K0), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.integer(max.iter), as.double(gamma), as.double(group.multiplier), as.integer(dfmax), as.integer(warn), as.integer(user.lambda))
    else fit <- .C("gpPathFit_binomial", double(nlambda), double(p*nlambda), integer(nlambda), double(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(n), as.integer(p), penalty, as.integer(J), as.integer(K1), as.integer(K0), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.double(0), as.integer(max.iter), as.double(gamma), as.double(tau), as.double(group.multiplier), as.integer(dfmax), as.integer(warn), as.integer(user.lambda))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    iter <- fit[[3]]
    df <- fit[[4]]
    loss <- fit[[5]]
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(b[p,])
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda")

  ## Unstandardize
  if (strtrim(penalty,2)=="gr") b <- unorthogonalize(b, XX, group)
  b <- unstandardize(b, center[nz], scale[nz])
  beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
  beta[1,] <- b[1,]
  beta[nz+1,] <- b[-1,]
  
  ## Names
  if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
  else varnames <- colnames(X)
  varnames <- c("(Intercept)",varnames)
  dimnames(beta) <- list(varnames, round(lambda,digits=4))
  
  structure(list(beta=beta,
                 family=family,
                 group=group.orig,
                 lambda=lambda,
                 alpha=alpha,
                 loss = loss,
                 n = length(y),
                 penalty=penalty,
                 df=df,
                 iter=iter),
            class = "grpreg")
}
