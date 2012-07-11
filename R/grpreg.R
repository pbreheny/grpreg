##grpreg <- function(X, y, group=1:ncol(X), family=c("gaussian","binomial"), penalty=c("geMCP","geLasso","gMCP","gBridge","gLasso"), nlambda=100, lambda, lambda.min=ifelse(n>p,.001,.05), lambda.max, alpha=1, tau=1/3, eps=.005, max.iter=1000, dfmax=p, delta=1e-8, gamma, group.multiplier=rep(1,J), verbose=FALSE, warn.conv=TRUE)
grpreg <- function(X, y, group=1:ncol(X), inner=c("L2","L1","MCP"), outer=c("L1","MCP","SCAD"), penalty, family=c("gaussian","binomial"), nlambda=100, lambda, lambda.min={if (n>p) .001 else .05}, alpha=1, eps=.005, max.iter=1000, dfmax=p, gamma=3, group.multiplier=rep(1,J), verbose=FALSE, warn.conv=TRUE, ...)
{
  ## Check for errors
  family <- match.arg(family)
  inner <- match.arg(inner)
  outer <- match.arg(outer)
  if (!missing(penalty)) {
    message <- switch(penalty,
                      gMCP="Use inner='MCP', outer='MCP' instead.",
                      gBridge="Use the gBridge() function instead.",
                      gLasso="Use inner='L2', outer='L1' instead.")
    stop(paste("Use of the argument 'penalty' is obsolete.  ", message, sep=""))
  }
  if (alpha > 1 | alpha < 0) stop("alpha must be in [0,1]")
  if (length(group)!=ncol(X)) stop("group does not match X")
  if (is.null(colnames(X))) colnames(X) <- paste("V",1:ncol(X),sep="")
  J <- max(group)
  K <- as.numeric(table(group))
  if (!(identical(as.integer(sort(unique(group))),as.integer(1:J)) | identical(as.integer(sort(unique(group))),as.integer(0:J)))) stop("Groups must be consecutively numbered 1,2,3,...")
  if (length(group.multiplier)!=J) stop("Length of group.multiplier must equal number of penalized groups")

  ## Standardize X
  XX <- if (inner=="L2") orthogonalize(X, group) else standardize(X)

  ## Setup lambda
  if (missing(lambda)) {
    lambda1 <- setupLambda(XX, y, group, family, penalty, lambda.min, nlambda, gamma, group.multiplier)
    lambda <- lambda1/alpha
  }
  l <- length(lambda)
  
  ## Fit
  path <- .C("gpPathFit", double(p*l), integer(l), double(l), as.double(XX), as.double(y), as.integer(group), family, as.integer(n), as.integer(p), as.integer(J), as.integer(K), penalty, as.double(lambda*alpha), as.integer(l), as.double(eps), as.integer(max.iter), as.integer(verbose), as.double(delta), as.double(gamma), as.double(tau), as.double((1-alpha)*lambda), as.double(group.multiplier), as.integer(dfmax), as.integer(warn.conv))
  
  ## Eliminate saturated lambda values, if any
  b <- matrix(path[[1]], nrow=p, byrow=T)
  ind <- !is.na(b[p,])
  b <- b[, ind, drop=FALSE]
  iter <- path[[2]][ind]
  lambda <- lambda[ind]
  df <- path[[3]][ind]
  iter <- path[[2]][ind]
  
  ## Unstandardize
  beta <- unstandardize(b,meanx,normx)
  
  ## Names
  if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
  else varnames <- colnames(X)
  varnames <- c("(Intercept)",varnames)
  dimnames(beta) <- list(varnames,round(lambda,digits=4))
  
  val <- list(beta=beta,
              family=family,
              group=group,
              lambda=lambda,
              alpha=alpha,
              loss = calcL(cbind(1,X),y,beta,family),
              n = length(y),
              penalty=penalty,
              df=df,
              iter=iter)
  class(val) <- "grpreg"
  return(val)
}
