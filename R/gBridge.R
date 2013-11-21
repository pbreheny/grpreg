gBridge <- function(X, y, group=1:ncol(X), family=c("gaussian","binomial"), nlambda=100, lambda, lambda.min={if (nrow(X) > ncol(X)) .001 else .05}, lambda.max, alpha=1, eps=.005, delta=1e-7, max.iter=1000, gamma=0.5, group.multiplier=rep(1,J), warn=TRUE) {
  ## Check for errors
  family <- match.arg(family)
  if (alpha > 1 | alpha < 0) stop("alpha must be in [0,1]")
  if (length(group)!=ncol(X)) stop("group does not match X")
  if (delta <= 0) stop("Delta must be a positive number")
  J <- max(group)
  if (!(identical(as.integer(sort(unique(group))),as.integer(1:J)) | identical(as.integer(sort(unique(group))),as.integer(0:J)))) stop("Groups must be consecutively numbered 1,2,3,...")
  if (length(group.multiplier)!=J) stop("Length of group.multiplier must equal number of penalized groups")

  ## Set up XX, yy, lambda
  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  multi <- FALSE
  if (is.matrix(y) && ncol(y) > 1) {
    multi <- TRUE
    m <- ncol(y)
    response.names <- if (is.null(colnames(y))) paste("Y",1:m,sep="") else colnames(y)
    y <- multiY(y)
    X <- multiX(X, m)
    group <- c(rep(0, m-1), rep(group, rep(m,length(group))))
    group.multiplier <- rep(1,J)
  }
  XX <- standardize(X)
  center <- attr(XX, "center")
  scale <- attr(XX, "scale")
  nz <- which(scale > 1e-6)
  zg <- setdiff(unique(group), unique(group[nz]))
  if (length(zg)) {
    J  <- J - length(zg)
    group.multiplier <- group.multiplier[-zg]
  }
  XX <- XX[ ,nz, drop=FALSE]
  group.orig <- group
  group <- group[nz]
  K <- as.numeric(table(group))  
  yy <- if (family=="gaussian") y - mean(y) else y
  if (missing(lambda)) {
    lambda <- setupLambda.gBridge(XX, yy, group, family, alpha, lambda.min, lambda.max, nlambda, gamma, group.multiplier)
  } else {
    nlambda <- length(lambda)
  }
  
  ## Fit
  n <- length(yy)
  p <- ncol(XX)
  K0 <- if (min(group)==0) K[1] else 0
  K1 <- if (min(group)==0) cumsum(K) else c(0, cumsum(K))
  if (family=="gaussian") {
    fit <- .C("gpPathFit_gaussian", double(p*nlambda), integer(nlambda), double(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(n), as.integer(p), "gBridge", as.integer(J), as.integer(K1), as.integer(K0), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.double(delta), as.integer(max.iter), as.double(gamma), as.double(0), as.integer(p), as.integer(J), as.double(group.multiplier), as.integer(TRUE))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1 ## Intercept
    loss <- fit[[4]]
  }
  if (family=="binomial") {
    fit <- .C("gpPathFit_binomial", double(nlambda), double(p*nlambda), integer(nlambda), double(nlambda), double(nlambda), as.double(XX), as.double(yy), as.integer(n), as.integer(p), "gBridge", as.integer(J), as.integer(K1), as.integer(K0), as.double(lambda*alpha), as.double(lambda*(1-alpha)), as.integer(nlambda), as.double(eps), as.double(delta), as.integer(max.iter), as.double(gamma), as.double(0), as.double(group.multiplier), as.integer(p), as.integer(J), as.integer(warn), as.integer(TRUE))
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
  b <- unstandardize(b, center[nz], scale[nz])
  beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
  beta[1,] <- b[1,]
  beta[nz+1,] <- b[-1,]
    
  ## Names
  varnames <- c("(Intercept)", xnames)
  if (multi) {
    beta[2:m,] <- sweep(beta[2:m,], 2, beta[1,], FUN="+")
    beta <- array(beta, dim=c(m, nrow(beta)/m, ncol(beta)))
    group.orig <- group.orig[-(1:(m-1))]
    dimnames(beta) <- list(response.names, varnames, round(lambda,digits=4))
  } else {
    dimnames(beta) <- list(varnames, round(lambda,digits=4))
  }
  
  structure(list(beta = beta,
                 family = family,
                 group = group.orig,
                 lambda = lambda,
                 alpha = alpha,
                 loss = loss,
                 n = length(y),
                 penalty = "gBridge",
                 df = df,
                 iter = iter),
            class = "grpreg")
}
