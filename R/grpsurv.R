#' Fit an group penalized survival model
#' 
#' Fit regularization paths for Cox models with grouped penalties over a grid
#' of values for the regularization parameter lambda.
#' 
#' The sequence of models indexed by the regularization parameter \code{lambda}
#' is fit using a coordinate descent algorithm.  In order to accomplish this,
#' the second derivative (Hessian) of the Cox partial log-likelihood is
#' diagonalized (see references for details).  The objective function is
#' defined to be \deqn{Q(\beta|X, y) = \frac{1}{n} L(\beta|X, y) + }{Q(\beta|X,
#' y) = (1/n)*L(\beta|X, y) + P(\beta, \lambda),}\deqn{
#' P_\lambda(\beta)}{Q(\beta|X, y) = (1/n)*L(\beta|X, y) + P(\beta, \lambda),}
#' where the loss function L is the negative partial log-likelihood (half the
#' deviance) from the Cox regression model.
#' [See here for more details](https://pbreheny.github.io/grpreg/articles/models.html).
#' 
#' Presently, ties are not handled by \code{grpsurv} in a particularly
#' sophisticated manner.  This will be improved upon in a future release of
#' \code{grpreg}.
#' 
#' @param X The design matrix.
#' @param y The time-to-event outcome, as a two-column matrix or
#' \code{\link[survival]{Surv}} object.  The first column should be time on
#' study (follow up time); the second column should be a binary variable with 1
#' indicating that the event has occurred and 0 indicating (right) censoring.
#' @param group A vector describing the grouping of the coefficients.  For
#' greatest efficiency and least ambiguity (see details), it is best if
#' \code{group} is a factor or vector of consecutive integers, although
#' unordered groups and character vectors are also allowed.  If there are
#' coefficients to be included in the model without being penalized, assign
#' them to group 0 (or \code{"0"}).
#' @param penalty The penalty to be applied to the model.  For group selection,
#' one of \code{grLasso}, \code{grMCP}, or \code{grSCAD}.  For bi-level
#' selection, one of \code{gel} or \code{cMCP}.  See below for details.
#' @param gamma Tuning parameter of the group or composite MCP/SCAD penalty
#' (see details).  Default is 3 for MCP and 4 for SCAD.
#' @param alpha \code{grpsurv} allows for both a group penalty and an L2
#' (ridge) penalty; \code{alpha} controls the proportional weight of the
#' regularization parameters of these two penalties.  The group penalties'
#' regularization parameter is \code{lambda*alpha}, while the regularization
#' parameter of the ridge penalty is \code{lambda*(1-alpha)}.  Default is 1: no
#' ridge penalty.
#' @param nlambda The number of lambda values.  Default is 100.
#' @param lambda.min The smallest value for lambda, as a fraction of
#' lambda.max.  Default is .001 if the number of observations is larger than
#' the number of covariates and .05 otherwise.
#' @param lambda A user-specified sequence of lambda values.  By default, a
#' sequence of values of length \code{nlambda} is computed automatically,
#' equally spaced on the log scale.
#' @param eps Convergence threshhold.  The algorithm iterates until the RMSD
#' for the change in linear predictors for each coefficient is less than
#' \code{eps}.  Default is \code{0.001}.
#' @param max.iter Maximum number of iterations (total across entire path).
#' Default is 10000.
#' @param dfmax Limit on the number of parameters allowed to be nonzero.  If
#' this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' @param gmax Limit on the number of groups allowed to have nonzero elements.
#' If this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' @param tau Tuning parameter for the group exponential lasso; defaults to
#' 1/3.
#' @param group.multiplier A vector of values representing multiplicative
#' factors by which each group's penalty is to be multiplied.  Often, this is a
#' function (such as the square root) of the number of predictors in each
#' group.  The default is to use the square root of group size for the group
#' selection methods, and a vector of 1's (i.e., no adjustment for group size)
#' for bi-level selection.
#' @param warn Return warning messages for failures to converge and model
#' saturation?  Default is TRUE.
#' @param returnX Return the standardized design matrix?  Default is FALSE.
#' @param ... Not used.
#' 
#' @returns An object with S3 class `"grpsurv"` containing:
#' \item{beta}{The fitted matrix of coefficients. The number of rows is equal to
#' the number of coefficients, and the number of columns is equal to `nlambda`.}
#' \item{group}{Same as above.}
#' \item{lambda}{The sequence of `lambda` values in the path.}
#' \item{penalty}{Same as above.}
#' \item{gamma}{Same as above.}
#' \item{alpha}{Same as above.}
#' \item{deviance}{The deviance of the fitted model at each value of `lambda`.}
#' \item{n}{The number of observations.}
#' \item{df}{A vector of length `nlambda` containing estimates of effective
#' number of model parameters all the points along the regularization path. For
#' details on how this is calculated, see Breheny and Huang (2009).}
#' \item{iter}{A vector of length `nlambda` containing the number of iterations
#' until convergence at each value of `lambda`.}
#' \item{group.multiplier}{A named vector containing the multiplicative constant
#' applied to each group's penalty.}
#' 
#' For Cox models, the following objects are also returned (and are necessary
#' to estimate baseline survival conditional on the estimated regression
#' coefficients), all of which are ordered by time on study (i.e., the ith row
#' of `W` does not correspond to the ith row of `X`):
#' 
#' \item{W}{Matrix of `exp(beta)` values for each subject over all `lambda`
#' values.}
#' \item{time}{Times on study.}
#' \item{fail}{Failure event indicator.}
#' 
#' @author Patrick Breheny
#' 
#' @seealso [plot.grpreg()], [predict.grpsurv()], [cv.grpsurv()]
#' 
#' @references
#' \itemize{
#' \item Breheny P and Huang J. (2009) Penalized methods for bi-level variable
#' selection. *Statistics and its interface*, **2**: 369-380.
#' \doi{10.4310/sii.2009.v2.n3.a10}
#' 
#' \item Huang J, Breheny P, and Ma S. (2012). A selective review of group
#' selection in high dimensional models. *Statistical Science*, **27**: 481-499.
#' \doi{10.1214/12-sts392}
#' 
#' \item Breheny P and Huang J. (2015) Group descent algorithms for nonconvex
#' penalized linear and logistic regression models with grouped predictors.
#' *Statistics and Computing*, **25**: 173-187. \doi{10.1007/s11222-013-9424-2}
#' 
#' \item Breheny P. (2015) The group exponential lasso for bi-level variable
#' selection. *Biometrics*, **71**: 731-740. \doi{10.1111/biom.12300}
#' 
#' \item Simon N, Friedman JH, Hastie T, and Tibshirani R. (2011)
#' Regularization Paths for Cox's Proportional Hazards Model via Coordinate
#' Descent. *Journal of Statistical Software*, **39**: 1-13.
#' \doi{10.18637/jss.v039.i05}
#' }
#' 
#' @examples
#' data(Lung)
#' X <- Lung$X
#' y <- Lung$y
#' group <- Lung$group
#' 
#' fit <- grpsurv(X, y, group)
#' plot(fit)
#' 
#' S <- predict(fit, X, type='survival', lambda=0.05)
#' plot(S, xlim=c(0,200))
#' @export

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
  loss <- -res[[4]]
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
                        deviance = 2 * loss,
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
