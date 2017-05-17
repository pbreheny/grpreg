grpreg <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),
                   family=c("gaussian","binomial", "poisson"), nlambda=100, lambda,
                   lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, log.lambda = TRUE,
                   alpha=1, eps=1e-4, max.iter=10000, dfmax=p, gmax=length(unique(group)),
                   gamma=ifelse(penalty=="grSCAD", 4, 3), tau=1/3, group.multiplier,
                   warn=TRUE, returnX=FALSE, ...) {

  # Deprecation support / error checking
  if (!missing(penalty)) {
    if (penalty[1]=="gBridge") stop("gBridge has been divorced from the grpreg function; use the gBridge() function instead")
    if (penalty[1]=="gMCP") {
      writeLines(strwrap("penalty='gMCP' is deprecated and may not be supported in future versions.  Use penalty='cMCP' instead."))
      penalty <- "cMCP"
    }
    if (penalty[1]=="gLasso") {
      writeLines(strwrap("You have specified penalty='gLasso'; grpreg is assuming you mean group lasso (penalty='grLasso')"))
      penalty <- "grLasso"
    }
  }
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (gamma <= 1 & penalty %in% c("grMCP", "cMCP")) stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="grSCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")

  # Construct XG, yy
  bilevel <- strtrim(penalty,2) != "gr"
  yy <- newY(y, family)
  m <- attr(yy, "m")
  XG <- newXG(X, group, group.multiplier, m, bilevel)
  if (nrow(XG$X) != length(yy)) stop("X and y do not have the same number of observations")

  # Setup lambda
  if (missing(lambda)) {
    lambda <- setupLambda(XG$X, yy, XG$g, family, penalty, alpha, lambda.min, log.lambda, nlambda, XG$m)
    lam.max <- lambda[1]
    user.lambda <- FALSE
  } else {
    lam.max <- -1
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  n <- length(yy)
  p <- ncol(XG$X)
  K <- as.numeric(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  if (K0) {
    lambda[1] <- lambda[1] + 1e-5
    user.lambda <- TRUE
  }
  if (family=="gaussian") {
    if (bilevel) fit <- .Call("lcdfit_gaussian", XG$X, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), XG$m, as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    else fit <- .Call("gdfit_gaussian", XG$X, yy, penalty, K1, K0, lambda, lam.max, alpha, eps, as.integer(max.iter), gamma, XG$m, as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1 ## Intercept
    loss <- fit[[4]]
  } else {
    if (family=="binomial") {
      if (bilevel) fit <- .Call("lcdfit_binomial", XG$X, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
      else fit <- .Call("gdfit_binomial", XG$X, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    } else if (family=="poisson") {
      if (bilevel) fit <- .Call("lcdfit_poisson", XG$X, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
      else fit <- .Call("gdfit_poisson", XG$X, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    }
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    iter <- fit[[3]]
    df <- fit[[4]]
    loss <- fit[[5]]
  }

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda")

  ## Unstandardize
  if (strtrim(penalty,2)=="gr") b <- unorthogonalize(b, XG$X, XG$g)
  b <- unstandardize(b, XG$center, XG$scale)
  beta <- matrix(0, nrow=m*(ncol(X)+1), ncol=length(lambda))
  beta[1,] <- b[1,]
  if (XG$reorder) {
    beta[XG$nz+1,] <- b[-1,]
    beta[-1,] <- beta[1+XG$ord.inv,]
  } else {
    beta[XG$nz+1,] <- b[-1,]
  }

  ## Names
  varnames <- c("(Intercept)", XG$names)
  if (m > 1) {
    beta[2:m,] <- sweep(beta[2:m,,drop=FALSE], 2, beta[1,], FUN="+")
    beta <- array(beta, dim=c(m, nrow(beta)/m, ncol(beta)))
    group <- group[-(1:(m-1))]
    dimnames(beta) <- list(colnames(yy), varnames, round(lambda,digits=4))
  } else {
    dimnames(beta) <- list(varnames, round(lambda,digits=4))
  }

  val <- structure(list(beta = beta,
                        family = family,
                        group = group,
                        lambda = lambda,
                        alpha = alpha,
                        loss = loss,
                        n = n,
                        penalty = penalty,
                        df = df,
                        iter = iter,
                        group.multiplier = XG$m),
                   class = "grpreg")
  if (returnX) {
    val$XG = XG
    val$y = yy
  } else if (family=="poisson") {
    val$y <- y
  }
  val
}
