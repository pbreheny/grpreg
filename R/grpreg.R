grpreg <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gBridge", "gLasso", "gMCP"),
                   family=c("gaussian","binomial", "poisson"), nlambda=100, lambda,
                   lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, alpha=1, eps=.001, max.iter=1000,
                   dfmax=p, gmax=length(unique(group)), gamma=ifelse(penalty=="grSCAD", 4, 3),
                   tau=1/3, group.multiplier, warn=TRUE, ...) {
  # Coersion
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (penalty=="gLasso") penalty <- "grLasso"
  if (penalty=="gMCP") penalty <- "cMCP"
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (is.matrix(y)) {
    for (j in 1:ncol(y)) {
      y[,j] <- coerceY(y[,j], family)
    }
  } else {
    y <- coerceY(y, family)
  }

  ## Error checking
  if (gamma <= 1 & penalty %in% c("grMCP", "cMCP")) stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="grSCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg")
  if (penalty=="gBridge") stop("gBridge has been divorced from the grpreg function; use the gBridge() function instead")
  if (length(group)!=ncol(X)) stop("group does not match X")

  # Reorder groups, if necessary
  grp <- reorderGroups(group, group.multiplier, strtrim(penalty,2)=="gr")
  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  if (grp$reorder) X <- X[,grp$ord]

  ## Set up XX, yy, lambda
  multi <- FALSE
  if (is.matrix(y) && ncol(y) > 1) {
    multi <- TRUE
    m <- ncol(y)
    response.names <- if (is.null(colnames(y))) paste("Y",1:m,sep="") else colnames(y)
    y <- multiY(y)
    X <- multiX(X, m)
    group <- grp$g <- c(rep(0, m-1), rep(grp$g, each=m))
    grp$m <- rep(grp$m, each=m)
  }
  std <- .Call("standardize", X)
  XX <- std[[1]]
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
  yy <- as.numeric(if (family=="gaussian") y - mean(y) else y)
  if (nrow(XX) != length(yy)) stop("X and y do not have the same number of observations")
  if (missing(lambda)) {
    lambda <- setupLambda(XX, yy, grp$g, family, penalty, alpha, lambda.min, nlambda, grp$m)
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
  if (K0) {
    lambda[1] <- lambda[1] + 1e-5
    user.lambda <- TRUE
  }
  if (family=="gaussian") {
    if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_gaussian", XX, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, grp$m, as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    else fit <- .Call("lcdfit_gaussian", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), grp$m, as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1 ## Intercept
    loss <- fit[[4]]
  }
  if (family=="binomial") {
    if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_binomial", XX, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    else fit <- .Call("lcdfit_binomial", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    iter <- fit[[3]]
    df <- fit[[4]]
    loss <- fit[[5]]
  }
  if (family=="poisson") {
    if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_poisson", XX, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    else fit <- .Call("lcdfit_poisson", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
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
  if (strtrim(penalty,2)=="gr") b <- unorthogonalize(b, XX, grp$g)
  b <- unstandardize(b, center[nz], scale[nz])
  beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
  beta[1,] <- b[1,]
  if (grp$reorder) {
    beta[nz+1,] <- b[-1,]
    beta[-1,] <- beta[1+grp$ord.inv,]
  } else {
    beta[nz+1,] <- b[-1,]
  }

  ## Names
  varnames <- c("(Intercept)", xnames)
  if (multi) {
    beta[2:m,] <- sweep(beta[2:m,,drop=FALSE], 2, beta[1,], FUN="+")
    beta <- array(beta, dim=c(m, nrow(beta)/m, ncol(beta)))
    group <- group[-(1:(m-1))]
    dimnames(beta) <- list(response.names, varnames, round(lambda,digits=4))
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
                        group.multiplier = grp$m),
                   class = "grpreg")
  if (family=="poisson") val$y <- y
  val
}
coerceY <- function(y, family) {
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }
  if (family=="binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data")
  if (family=="binomial" & !identical(sort(unique(y)), 0:1)) y <- as.numeric(y==max(y))
  y
}
