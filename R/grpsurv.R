grpsurv <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),
                   gamma=ifelse(penalty=="grSCAD", 4, 3), alpha=1, nlambda=100, lambda,
                   lambda.min={if (nrow(X) > ncol(X)) 0.001 else .05}, eps=.001, max.iter=10000,
                   dfmax=p, gmax=length(unique(group)), tau=1/3,
                   group.multiplier, warn=TRUE, returnX=FALSE, ...) {

  # Deprecation support / error checking
  if (penalty[1]=="gBridge") stop("gBridge has been divorced from the grpreg function; use the gBridge() function instead")
  if (penalty[1]=="gMCP") {
    writeLines(strwrap("penalty='gMCP' is deprecated and may not be supported in future versions.  Use penalty='cMCP' instead."))
    penalty <- "cMCP"
  }
  if (penalty[1]=="gLasso") {
    writeLines(strwrap("You have specified penalty='gLasso'; grpreg is assuming you mean group lasso (penalty='grLasso')"))
    penalty <- "grLasso"
  }
  penalty <- match.arg(penalty)
  if (gamma <= 1 & penalty %in% c("grMCP", "cMCP")) stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="grSCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")

  # Construct XG, Y
  bilevel <- strtrim(penalty,2) != "gr"
  Y <- newS(y)
  XG <- newXG(X[Y$ind,], group, group.multiplier, 1, bilevel)
  if (nrow(XG$X) != length(Y$fail)) stop("X and y do not have the same number of observations")

  # Set up lambda
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XG$X, Y$time, Y$fail, XG$g, penalty, alpha, lambda.min, nlambda, XG$m)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  n <- length(Y$time)
  p <- ncol(XG$X)
  K <- as.numeric(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  if (bilevel) {
    res <- .Call("lcdfit_cox", XG$X, Y$fail, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter),
                 XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
  } else {
    res <- .Call("gdfit_cox", XG$X, Y$fail, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter),
                 as.double(gamma), XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
  }
  b <- matrix(res[[1]], p, nlambda)
  iter <- res[[2]]
  df <- res[[3]]
  loss <- -1*res[[4]]
  Eta <- matrix(res[[5]], n, nlambda)

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for some values of lambda")

  ## Unstandardize
  if (!bilevel) b <- unorthogonalize(b, XG$X, XG$g, intercept=FALSE)
  b <- b/XG$scale
  beta <- matrix(0, nrow=ncol(X), ncol=length(lambda))
  if (XG$reorder) {
    beta[XG$nz,] <- b
    beta <- beta[XG$ord.inv,]
  } else {
    beta[XG$nz,] <- b
  }

  ## Names
  dimnames(beta) <- list(XG$names, round(lambda,digits=4))

  ## Output
  val <- structure(list(beta = beta,
                        group = group,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        loss = loss,
                        n = n,
                        df = df,
                        iter = iter,
                        group.multiplier = XG$m,
                        time = Y$time,
                        fail = Y$fail,
                        W = exp(Eta)),
                   class = c("grpsurv", "grpreg"))
  if (returnX) {
    val$XG <- XG
  }
  val
}
