cv.grpsurv <- function(X, y, group, ..., nfolds=10, seed, cv.ind, returnY=FALSE, trace=FALSE) {
  if (!missing(seed)) set.seed(seed)

  # Complete data fit
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$group <- group
  fit.args$returnX <- TRUE
  fit <- do.call("grpsurv", fit.args)

  # Get standardized X, y
  X <- fit$XG$X
  y <- cbind(fit$time, fit$fail)
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX) fit$X <- NULL

  # Set up folds
  n <- nrow(X)
  if (!missing(seed)) set.seed(seed)
  if (missing(cv.ind)) cv.ind <- ceiling(sample(1:n)/n*nfolds)
  Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))

  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- fit$XG$g
  cv.args$group.multiplier <- fit$XG$m
  cv.args$warn <- FALSE

  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #",i,sep="","\n")
    res <- cvf.surv(i, X, y, cv.ind, cv.args)
    Y[cv.ind==i, 1:res$nl] <- res$yhat
  }

  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(Y), 2, all))
  Y <- Y[,ind]
  lambda <- fit$lambda[ind]

  ## Return
  cve <- as.numeric(loss.grpsurv(y, Y))
  min <- which.min(cve)

  val <- list(cve=cve, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=cve[1])
  if (returnY) val$Y <- Y
  structure(val, class=c("cv.grpsurv", "cv.grpreg"))
}
cvf.surv <- function(i, XX, y, cv.ind, cv.args) {
  cv.args$X <- XX[cv.ind!=i, , drop=FALSE]
  cv.args$y <- y[cv.ind!=i,]
  fit.i <- do.call("grpsurv", cv.args)

  X2 <- XX[cv.ind==i, , drop=FALSE]
  y2 <- y[cv.ind==i,]
  nl <- length(fit.i$lambda)
  yhat <- predict(fit.i, X2)

  list(nl=length(fit.i$lambda), yhat=yhat)#, loss=loss)
}
