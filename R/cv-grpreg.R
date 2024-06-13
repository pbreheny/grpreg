#' Cross-validation for grpreg/grpsurv
#' 
#' Performs k-fold cross validation for penalized regression models with
#' grouped covariates over a grid of values for the regularization parameter
#' lambda.
#' 
#' The function calls [grpreg()] or [grpsurv()] `nfolds` times, each
#' time leaving out 1/`nfolds` of the data.  The cross-validation error is
#' based on the deviance;
#' [see here for more details](https://pbreheny.github.io/grpreg/articles/models.html).
#' 
#' For Gaussian and Poisson responses, the folds are chosen according to simple
#' random sampling.  For binomial responses, the numbers for each outcome class
#' are balanced across the folds; i.e., the number of outcomes in which
#' `y` is equal to 1 is the same for each fold, or possibly off by 1 if
#' the numbers do not divide evenly.  This approach is used for Cox regression
#' as well to balance the amount of censoring cross each fold.
#' 
#' For Cox models, `cv.grpsurv` uses the approach of calculating the full
#' Cox partial likelihood using the cross-validated set of linear predictors.
#' Other approaches to cross-validation for the Cox regression model have been
#' proposed in the literature; the strengths and weaknesses of the various
#' methods for penalized regression in the Cox model are the subject of current
#' research.  A simple approximation to the standard error is provided,
#' although an option to bootstrap the standard error (`se='bootstrap'`)
#' is also available.
#' 
#' As in [grpreg()], seemingly unrelated regressions/multitask learning can
#' be carried out by setting `y` to be a matrix, in which case groups are
#' set up automatically (see [grpreg()] for details), and
#' cross-validation is carried out with respect to rows of `y`.  As
#' mentioned in the details there, it is recommended to standardize the
#' responses prior to fitting.
#' 
#' @aliases cv.grpreg cv.grpsurv
#' 
#' @param X The design matrix, as in [grpreg()]/[grpsurv()].
#' @param y The response vector (or matrix), as in [grpreg()]/[grpsurv()].
#' @param group The grouping vector, as in [grpreg()]/[grpsurv()].
#' @param ... Additional arguments to [grpreg()]/[grpsurv()].
#' @param nfolds The number of cross-validation folds.  Default is 10.
#' @param seed You may set the seed of the random number generator in order to
#' obtain reproducible results.
#' @param fold Which fold each observation belongs to.  By default the
#' observations are randomly assigned.
#' @param returnY Should cv.grpreg()/cv.grpsurv() return the fitted
#' values from the cross-validation folds?  Default is FALSE; if TRUE, this
#' will return a matrix in which the element for row i, column j is the fitted
#' value for observation i from the fold in which observation i was excluded
#' from the fit, at the jth value of lambda.  NOTE: For `cv.grpsurv()`, the
#' rows of `Y` are ordered by time on study, and therefore will not
#' correspond to the original order of observations pased to `cv.grpsurv`.
#' @param trace If set to TRUE, cv.grpreg will inform the user of its progress
#' by announcing the beginning of each CV fold.  Default is FALSE.
#' @param se For `cv.grpsurv()`, the method by which the cross-valiation
#' standard error (CVSE) is calculated.  The 'quick' approach is based on a
#' rough approximation, but can be calculated more or less instantly.  The
#' 'bootstrap' approach is more accurate, but requires additional computing
#' time.
#' 
#' @returns An object with S3 class \code{"cv.grpreg"} containing:
#' \item{cve}{The error for each value of \code{lambda}, averaged across the
#' cross-validation folds.} \item{cvse}{The estimated standard error associated
#' with each value of for \code{cve}.} \item{lambda}{The sequence of
#' regularization parameter values along which the cross-validation error was
#' calculated.} \item{fit}{The fitted \code{grpreg} object for the whole data.}
#' \item{fold}{The fold assignments for cross-validation for each observation;
#' note that for \code{cv.grpsurv}, these are in terms of the ordered
#' observations, not the original observations.} \item{min}{The index of
#' \code{lambda} corresponding to \code{lambda.min}.} \item{lambda.min}{The
#' value of \code{lambda} with the minimum cross-validation error.}
#' \item{null.dev}{The deviance for the intercept-only model.}
#' \item{pe}{If `family="binomial"`, the cross-validation prediction error for
#' each value of `lambda`.}
#' 
#' @author Patrick Breheny
#' 
#' @seealso [grpreg()], [plot.cv.grpreg()], [summary.cv.grpreg()],
#' [predict.cv.grpreg()]
#' 
#' @examples
#' \dontshow{set.seed(1)}
#' data(Birthwt)
#' X <- Birthwt$X
#' y <- Birthwt$bwt
#' group <- Birthwt$group
#' 
#' cvfit <- cv.grpreg(X, y, group)
#' plot(cvfit)
#' summary(cvfit)
#' coef(cvfit) ## Beta at minimum CVE
#' 
#' cvfit <- cv.grpreg(X, y, group, penalty="gel")
#' plot(cvfit)
#' summary(cvfit)
#' 
#' @export cv.grpreg

cv.grpreg <- function(X, y, group=1:ncol(X), ..., nfolds=10, seed, fold, returnY=FALSE, trace=FALSE) {

  # Complete data fit
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  if (!inherits(X, "expandedMatrix")) fit.args$group <- group
  fit.args$returnX <- TRUE
  if ('penalty' %in% names(fit.args) && fit.args$penalty == 'gBridge') {
    fit <- do.call("gBridge", fit.args)
  } else {
    fit <- do.call("grpreg", fit.args)
  }

  # Get y, standardized X
  XG <- fit$XG
  X <- XG$X
  m <- attr(fit$y, "m")
  y <- if (fit$family=="gaussian") fit$y - attr(fit$y, "mean") else fit$y
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX) fit$XG <- NULL

  # Set up folds
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  n <- length(y)  
  if (missing(fold)) {
    if (m > 1) {
      nn <- n/m
      fold_ <- sample(1:nn %% (nfolds))
      fold_[fold_==0] <- nfolds
      fold <- rep(fold_, each=m)
    } else if (fit$family=="binomial") {
      ind1 <- which(y==1)
      ind0 <- which(y==0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      fold1 <- 1:n1 %% nfolds
      fold0 <- (n1 + 1:n0) %% nfolds
      fold1[fold1==0] <- nfolds
      fold0[fold0==0] <- nfolds
      fold <- integer(n)
      fold[y==1] <- sample(fold1)
      fold[y==0] <- sample(fold0)
    } else {
      fold <- sample(1:n %% nfolds)
      fold[fold==0] <- nfolds
    }
  } else {
    nfolds <- max(fold)
  }

  # Do cross-validation
  E <- Y <- matrix(NA, nrow=length(y), ncol=length(fit$lambda))
  if (fit$family=="binomial") PE <- E
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- XG$g
  cv.args$group.multiplier <- XG$m
  cv.args$warn <- FALSE
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #", i, sep="","\n")
    res <- cvf(i, X, y, fold, cv.args)
    Y[fold==i, 1:res$nl] <- res$yhat
    E[fold==i, 1:res$nl] <- res$deviance
    if (fit$family=="binomial") PE[fold==i, 1:res$nl] <- res$pe
  }

  # Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  # Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)
  null.dev <- calc_null_dev(X, y, group=XG$g, family=fit$family)

  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, fold=fold, min=min, lambda.min=lambda[min], null.dev=null.dev)
  if (fit$family=="binomial") val$pe <- apply(PE[, ind], 2, mean)
  if (returnY) {
    if (fit$family=="gaussian") val$Y <- Y + attr(y, "mean")
    else val$Y <- Y
  }
  structure(val, class="cv.grpreg")
}
cvf <- function(i, X, y, fold, cv.args) {
  cv.args$X <- X[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  if ('penalty' %in% names(cv.args) && cv.args$penalty == 'gBridge') {
    fit.i <- do.call("gBridge", cv.args)
  } else {
    fit.i <- do.call("grpreg", cv.args)
  }

  X2 <- X[fold==i, , drop=FALSE]
  y2 <- y[fold==i]
  yhat <- matrix(predict(fit.i, X2, type="response"), length(y2))
  deviance <- deviance_grpreg(y2, yhat, fit.i$family)
  pe <- if (fit.i$family=="binomial") {(yhat < 0.5) == y2} else NULL
  list(deviance=deviance, pe=pe, nl=length(fit.i$lambda), yhat=yhat)
}
