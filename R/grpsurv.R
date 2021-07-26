grpsurv <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),
                   gamma=ifelse(penalty=="grSCAD", 4, 3), alpha=1, nlambda=100, lambda,
                   lambda.min={if (nrow(X) > ncol(X)) 0.001 else .05}, eps=.001, max.iter=10000,
                   dfmax=p, gmax=length(unique(group)), tau=1/3,
                   group.multiplier, warn=TRUE, returnX=FALSE, ...) {

  # Deprecation support / error checking
  if (penalty[1]=="gBridge") stop("gBridge has been divorced from the grpreg function; use the gBridge() function instead", call.=FALSE)
  if (penalty[1]=="gMCP") {
    writeLines(strwrap("penalty='gMCP' is deprecated and may not be supported in future versions.  Use penalty='cMCP' instead."))
    penalty <- "cMCP"
  }
  if (penalty[1]=="gLasso") {
    writeLines(strwrap("You have specified penalty='gLasso'; grpreg is assuming you mean group lasso (penalty='grLasso')"))
    penalty <- "grLasso"
  }
  penalty <- match.arg(penalty)
  if (gamma <= 1 & penalty %in% c("grMCP", "cMCP")) stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="grSCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0, 1]", call.=FALSE)

  # Check for expandedMatix
  if(inherits(X, "expandedMatrix")) {
    expanded <- TRUE
    group <- X$group
    knots <- X$knots
    boundary <- X$boundary
    degree <- X$degree
    originalx <- X$originalx
    type <- X$type
    X <- X$X
  } else {
    expanded <- FALSE
  }
  
  # Construct XG, Y
  bilevel <- strtrim(penalty, 2) != "gr"
  Y <- newS(y)
  XG <- newXG(X[Y$ind, , drop=FALSE], group, group.multiplier, 1, bilevel)
  if (nrow(XG$X) != length(Y$fail)) stop("X and y do not have the same number of observations", call.=FALSE)

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
  K <- as.integer(table(XG$g))
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
  loss <- -2*res[[4]]
  Eta <- matrix(res[[5]], n, nlambda)

  # Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  loss <- loss[ind]
  Eta <- Eta[, ind, drop=FALSE]
  if (iter[1] == max.iter) stop("Algorithm failed to converge for any values of lambda.  This indicates a combination of (a) an ill-conditioned feature matrix X and (b) insufficient penalization.  You must fix one or the other for your model to be identifiable.", call.=FALSE)
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda", call.=FALSE)

  # Unstandardize
  if (!bilevel) b <- unorthogonalize(b, XG$X, XG$g, intercept=FALSE)
  if (XG$reorder) b <- b[XG$ord.inv,]
  beta <- matrix(0, nrow=length(XG$scale), ncol=ncol(b))
  beta[XG$nz,] <- b / XG$scale[XG$nz]

  # Names
  dimnames(beta) <- list(XG$names, round(lambda, digits=4))
  colnames(Eta) <- round(lambda, digits=4)

  # Output
  val <- structure(list(beta = beta,
                        group = factor(group),
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
                        order = Y$ind,
                        linear.predictors = sweep(Eta, 2, colMeans(Eta), '-')),
                   class = c("grpsurv", "grpreg"))
  if (returnX) val$XG <- XG
  if (expanded) {
    val$meta <- list(knots = knots,
                     boundary = boundary,
                     degree = degree,
                     originalx = originalx,
                     type = type,
                     X = X)
    attr(val, "class") <- c("grpsurv", "grpreg", "expanded")
  }
  val
}
