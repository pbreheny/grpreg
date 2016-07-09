grpsurv <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gBridge", "gLasso", "gMCP"),
                   gamma=ifelse(penalty=="grSCAD", 4, 3), alpha=1, nlambda=100, lambda,
                   lambda.min={if (nrow(X) > ncol(X)) 0.001 else .05}, eps=.001, max.iter=1000,
                   dfmax=p, gmax=length(unique(group)), tau=1/3,
                   group.multiplier, warn=TRUE, returnX=FALSE, ...) {
  # Coersion
  penalty <- match.arg(penalty)
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (class(y) != "matrix") {
    tmp <- try(y <- as.matrix(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must be a matrix or able to be coerced to a matrix")
    if (ncol(y)!=2) stop("y must have two columns for survival data: time-on-study and a censoring indicator")
  }
  if (storage.mode(y)=="integer") storage.mode(y) <- "double"

  # Error checking
  if (gamma <= 1 & penalty %in% c("grMCP", "cMCP")) stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="grSCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to grpreg")
  if (penalty=="gBridge") stop("Survival modeling not yet available for group Bridge")
  if (length(group)!=ncol(X)) stop("group does not match X")

  # Reorder groups, if necessary
  grp <- reorderGroups(group, group.multiplier, strtrim(penalty,2)=="gr")
  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  if (grp$reorder) X <- X[,grp$ord]

  # Set up XX, yy, lambda
  ind <- order(y[,1])
  yy <- as.numeric(y[ind,1])
  Delta <- y[ind,2]
  std <- .Call("standardize", X)
  XX <- std[[1]][ind,,drop=FALSE]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)
  zg <- setdiff(unique(grp$g), unique(grp$g[nz]))
  if (length(zg)) grp$m <- grp$m[-zg]
  XX <- XX[ ,nz, drop=FALSE]
  grp$g <- grp$g[nz]
  if (strtrim(penalty,2)=="gr") {
    XX <- orthogonalize(XX, grp$g)
    grp$g <- attr(XX, "group")
  }
  K <- as.numeric(table(grp$g))
  if (nrow(XX) != length(yy)) stop("X and y do not have the same number of observations")
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XX, yy, Delta, grp$g, penalty, alpha, lambda.min, nlambda, grp$m)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  n <- length(yy)
  p <- ncol(XX)
  K0 <- as.integer(if (min(grp$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(grp$g)==0) cumsum(K) else c(0, cumsum(K)))
  if (strtrim(penalty,2)=="gr") {
    res <- .Call("gdfit_cox", XX, Delta, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter),
                 as.double(gamma), grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
  } else {
    res <- .Call("lcdfit_cox", XX, Delta, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter),
                 grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
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
  if (strtrim(penalty,2)=="gr") b <- unorthogonalize(b, XX, grp$g, intercept=FALSE)
  b <- b/scale[nz]
  beta <- matrix(0, nrow=ncol(X), ncol=length(lambda))
  if (grp$reorder) {
    beta[nz,] <- b
    beta <- beta[grp$ord.inv,]
  } else {
    beta[nz,] <- b
  }

  ## Names
  dimnames(beta) <- list(xnames, round(lambda,digits=4))

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
                        group.multiplier = grp$m),
                   class = c("grpsurv", "grpreg"))
  val$W <- exp(Eta)
  val$time <- yy
  val$fail <- Delta
  if (returnX) {
    val$X <- XX
    val$center <- center
    val$scale <- scale
    val$y <- yy
  }
  val
}
