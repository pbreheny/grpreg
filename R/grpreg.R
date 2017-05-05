grpreg <- function(X, y, group=1:ncol(X),
                   penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),
                   family=c("gaussian","binomial", "poisson"), nlambda=100, lambda,
                   lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, log.lambda = TRUE,
                   alpha=1, eps=.001, max.iter=10000,
                   dfmax=p, gmax=length(unique(group)), gamma=ifelse(penalty=="grSCAD", 4, 3),
                   tau=1/3, group.multiplier, warn=TRUE, return.time=TRUE, ...) {

  # Deprecation support / error checking
  if (penalty=="gBridge") stop("gBridge has been divorced from the grpreg function; use the gBridge() function instead")
  if (penalty=="gMCP") {
    writeLines(strwrap("penalty='gMCP' is deprecated and may not be supported in future versions.  Use penalty='cMCP' instead."))
    penalty <- "cMCP"
  }
  if (penalty=="gLasso") {
    writeLines(strwrap("You have specified penalty='gLasso'; grpreg is assuming you mean group lasso (penalty='grLasso')"))
    penalty <- "grLasso"
  }
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  screen <- match.arg(screen)
  if (gamma <= 1 & penalty %in% c("grMCP", "cMCP")) stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="grSCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0,1]")

  # Construct XG, yy
  yy <- newY(y, family)
  m <- attr(yy, "m")
  XG <- newXG(X, group, group.multiplier, m)

  if (nrow(XX) != length(yy)) stop("X and y do not have the same number of observations")
  if (missing(lambda)) {
    lambda <- setupLambda(XX, yy, grp$g, family, penalty, alpha, lambda.min, log.lambda, nlambda, grp$m)
    lam.max <- lambda[1]
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
    if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_gaussian", XX, yy, penalty, K1, K0, lambda, lam.max, alpha, eps, as.integer(max.iter), gamma, grp$m, as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    else fit <- .Call("lcdfit_gaussian", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), grp$m, as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter <- fit[[2]]
    df <- fit[[3]] + 1 ## Intercept
    loss <- fit[[4]]
  } else {
    if (family=="binomial") {
      if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_binomial", XX, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
      else fit <- .Call("lcdfit_binomial", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    } else if (family=="poisson") {
      if (strtrim(penalty,2)=="gr") fit <- .Call("gdfit_poisson", XX, yy, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
      else fit <- .Call("lcdfit_poisson", XX, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), grp$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
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
                        group.multiplier = grp$m),
                   class = "grpreg")
  if (screen == 'SSR' || screen == 'SEDPP' || screen == 'SSR-BEDPP') val$rejections <- rejections
  if (screen == 'SSR-BEDPP') val$safe_rejections <- safe_rejections
  if (family=="poisson") val$y <- y
  if (return.time) val$time <- as.numeric(time['elapsed'])

  val
}
