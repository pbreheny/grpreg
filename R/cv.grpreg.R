cv.grpreg <- function(X, y, group=1:ncol(X), ..., nfolds=10, seed, cv.ind, returnY=FALSE, trace=FALSE) {

  # Complete data fit
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$group <- group
  fit.args$returnX <- TRUE
  fit <- do.call("grpreg", fit.args)

  # Get standardized X, y
  X <- fit$XG$X
  y <- fit$y
  m <- attr(fit$y, "m")
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX) fit$XG <- NULL

  # Set up folds
  if (!missing(seed)) set.seed(seed)
  n <- length(y)
  if (missing(cv.ind)) {
    if (m > 1) {
      nn <- n/m
      cv.ind <- rep(ceiling(sample(1:nn)/nn*nfolds), each=m)
    } else if (fit$family=="binomial" & (min(table(y)) > nfolds)) {
      ind1 <- which(y==1)
      ind0 <- which(y==0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
      cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
      cv.ind <- numeric(n)
      cv.ind[y==1] <- cv.ind1
      cv.ind[y==0] <- cv.ind0
    } else {
      cv.ind <- ceiling(sample(1:n)/n*nfolds)
    }
  }

  # Do cross-validation
  E <- Y <- matrix(NA, nrow=length(y), ncol=length(fit$lambda))
  if (fit$family=="binomial") PE <- E
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- fit$XG$g
  cv.args$group.multiplier <- fit$XG$m
  cv.args$warn <- FALSE
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #",i,sep="","\n")
    res <- cvf(i, X, y, cv.ind, cv.args)
    Y[cv.ind==i, 1:res$nl] <- res$yhat
    E[cv.ind==i, 1:res$nl] <- res$loss
    if (fit$family=="binomial") PE[cv.ind==i, 1:res$nl] <- res$pe
  }

  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[,ind]
  lambda <- fit$lambda[ind]

  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)
  null.dev <- calcNullDev(X, y, group=fit$XG$g, family=fit$family)

  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=null.dev)
  if (fit$family=="binomial") val$pe <- apply(PE[,ind], 2, mean)
  if (returnY) {
    if (fit$family=="gaussian") val$Y <- Y + attr(y, "mean")
    else val$Y <- Y
  }
  structure(val, class="cv.grpreg")
}
cvf <- function(i, X, y, cv.ind, cv.args) {
  cv.args$X <- X[cv.ind!=i, , drop=FALSE]
  cv.args$y <- y[cv.ind!=i]
  fit.i <- do.call("grpreg", cv.args)

  X2 <- X[cv.ind==i, , drop=FALSE]
  y2 <- y[cv.ind==i]
  yhat <- predict(fit.i, X2, type="response")
  loss <- loss.grpreg(y2, yhat, fit.i$family)
  pe <- if (fit.i$family=="binomial") {(yhat < 0.5) == y2} else NULL
  list(loss=loss, pe=pe, nl=length(fit.i$lambda), yhat=yhat)
}
