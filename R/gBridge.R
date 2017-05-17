gBridge <- function(X, y, group=1:ncol(X), family=c("gaussian","binomial","poisson"), nlambda=100, lambda,
                    lambda.min={if (nrow(X) > ncol(X)) .001 else .05}, lambda.max, alpha=1, eps=.001, delta=1e-7,
                    max.iter=10000, gamma=0.5, group.multiplier, warn=TRUE) {
  # Error checking
  family <- match.arg(family)
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg")
  if (length(group)!=ncol(X)) stop("group does not match X")
  if (delta <= 0) stop("Delta must be a positive number")

  # Construct XG, yy
  yy <- newY(y, family)
  m <- attr(yy, "m")
  XG <- newXG(X, group, group.multiplier, m, TRUE)
  if (nrow(XG$X) != length(yy)) stop("X and y do not have the same number of observations")

  ## Set up lambda
  if (missing(lambda)) {
    lambda <- setupLambda.gBridge(XG$X, yy, XG$g, family, alpha, lambda.min, lambda.max, nlambda, gamma, XG$m)
  } else {
    nlambda <- length(lambda)
  }

  ## Fit
  n <- length(yy)
  p <- ncol(XG$X)
  K <- as.numeric(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  if (family=="gaussian") {
    fit <- .Call("lcdfit_gaussian", XG$X, yy, "gBridge", K1, K0, lambda, alpha, eps, delta, gamma, 0, as.integer(max.iter), as.double(XG$m), as.integer(p), as.integer(max(XG$g)), as.integer(TRUE))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1 ## Intercept
    loss <- fit[[4]]
  } else {
    if (family=="binomial") {
      fit <- .Call("lcdfit_binomial", XG$X, yy, "gBridge", K1, K0, lambda, alpha, eps, delta, gamma, 0, as.integer(max.iter), as.double(XG$m), as.integer(p), as.integer(max(XG$g)), as.integer(warn), as.integer(TRUE))
    } else if (family=="poisson") {
      fit <- .Call("lcdfit_poisson", XG$X, yy, "gBridge", K1, K0, lambda, alpha, eps, delta, gamma, 0, as.integer(max.iter), as.double(XG$m), as.integer(p), as.integer(max(XG$g)), as.integer(warn), as.integer(TRUE))
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
                        penalty = "gBridge",
                        df = df,
                        iter = iter,
                        group.multiplier = XG$m),
                   class = "grpreg")
  if (family=="poisson") val$y <- y
  val
}
