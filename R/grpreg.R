##grpreg <- function(X, y, group=1:ncol(X), family=c("gaussian","binomial"), penalty=c("geMCP","geLasso","gMCP","gBridge","gLasso"), nlambda=100, lambda, lambda.min=ifelse(n>p,.001,.05), lambda.max, alpha=1, tau=1/3, eps=.005, max.iter=1000, dfmax=p, delta=1e-8, gamma, group.multiplier=rep(1,J), verbose=FALSE, warn.conv=TRUE)
grpreg <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gMCP", "gBridge", "gLasso"), family=c("gaussian","binomial"), nlambda=100, lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, alpha=1, eps=.005, max.iter=1000, dfmax=p, gamma=3, group.multiplier=rep(1,J), warn=TRUE, ...)
{
  ## Check for errors
  family <- match.arg(family)
  pen <- match.arg(penalty)
  if (penalty=="gLasso") penalty <- "grLasso"
  if (penalty=="gBridge") gBridge(...)
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")
  if (length(group)!=ncol(X)) stop("group does not match X")
  J <- max(group)
  K <- as.numeric(table(group))
  if (!(identical(as.integer(sort(unique(group))),as.integer(1:J)) | identical(as.integer(sort(unique(group))),as.integer(0:J)))) stop("Groups must be consecutively numbered 1,2,3,...")
  if (length(group.multiplier)!=J) stop("Length of group.multiplier must equal number of penalized groups")

  ## Set up X, y, lambda
  XX <- standardize(X)
  center <- attr(XX, "center")
  scale <- attr(XX, "scale")
  if (strtrim("penalty",2)=="gr") XX <- orthogonalize(XX)
  yy <- if (family=="gaussian") y - mean(y) else y
  lambda <- setupLambda(XX, yy, group, family, penalty, alpha, lambda.min, nlambda, gamma, group.multiplier)

  ## Fit
  n <- length(yy)
  p <- ncol(XX)
  tau <- 1/3
  if (family=="gaussian") {
    fit <- .C("gpPathFit_gaussian", double(p*nlambda), integer(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(group), as.integer(n), as.integer(p), penalty, as.integer(J), as.integer(K), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.integer(max.iter), as.double(gamma), as.double(tau), as.integer(dfmax), as.integer(group.multiplier))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    ## loss <- fit[[2]]
    iter <- fit[[2]]
    df <- fit[[3]]
  }
  if (family=="binomial") {
    fit <- .C("gpPathFit_binomial", double(nlambda), double(p*nlambda), integer(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(group), as.integer(n), as.integer(p), penalty, as.integer(J), as.integer(K), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.integer(max.iter), as.double(gamma), as.double(tau), as.double(group.multiplier), as.integer(dfmax), as.integer(warn))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    ## loss <- fit[[3]]
    iter <- fit[[3]]
    df <- fit[[4]]
  }
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(b[p,])
  b <- b[,ind,drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  
  ## Unstandardize
  if (strtrim("penalty",2)=="gr") beta <- unorthogonalize(b, XX, group)
  beta <- unstandardize(b, center, scale)
  
  ## print(mean(y) + XX %*% b[-1,])
  ## print(cbind(1,X) %*% beta)
  
  ## Names
  if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
  else varnames <- colnames(X)
  varnames <- c("(Intercept)",varnames)
  dimnames(beta) <- list(varnames, round(lambda,digits=4))
  
  structure(list(beta=beta,
                 family=family,
                 group=group,
                 lambda=lambda,
                 alpha=alpha,
                 loss = calcL(cbind(1,X),y,beta,family),
                 n = length(y),
                 penalty=penalty,
                 df=df,
                 iter=iter),
            class = "grpreg")
}
