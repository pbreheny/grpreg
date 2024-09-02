#' Fit a group penalized regression path
#' 
#' Fit regularization paths for models with grouped penalties over a grid of
#' values for the regularization parameter lambda. Fits linear and logistic
#' regression models.
#' 
#' There are two general classes of methods involving grouped penalties: those
#' that carry out bi-level selection and those that carry out group selection.
#' Bi-level means carrying out variable selection at the group level as well as
#' the level of individual covariates (i.e., selecting important groups as well
#' as important members of those groups).  Group selection selects important
#' groups, and not members within the group -- i.e., within a group,
#' coefficients will either all be zero or all nonzero.  The `grLasso`,
#' `grMCP`, and `grSCAD` penalties carry out group selection, while
#' the `gel` and `cMCP` penalties carry out bi-level selection.  For
#' bi-level selection, see also the [gBridge()] function.  For
#' historical reasons and backwards compatibility, some of these penalties have
#' aliases; e.g., `gLasso` will do the same thing as `grLasso`, but
#' users are encouraged to use `grLasso`.
#' 
#' Please note the distinction between \code{grMCP} and \code{cMCP}.  The
#' former involves an MCP penalty being applied to an L2-norm of each group.
#' The latter involves a hierarchical penalty which places an outer MCP penalty
#' on a sum of inner MCP penalties for each group, as proposed in Breheny &
#' Huang, 2009.  Either penalty may be referred to as the "group MCP",
#' depending on the publication.  To resolve this confusion, Huang et al.
#' (2012) proposed the name "composite MCP" for the \code{cMCP} penalty.
#' 
#' For more information about the penalties and their properties, please
#' consult the references below, many of which contain discussion, case
#' studies, and simulation studies comparing the methods.  If you use
#' \code{grpreg} for an analysis, please cite the appropriate reference.
#' 
#' In keeping with the notation from the original MCP paper, the tuning
#' parameter of the MCP penalty is denoted 'gamma'.  Note, however, that in
#' Breheny and Huang (2009), \code{gamma} is denoted 'a'.
#' 
#' The objective function for \code{grpreg} optimization is defined to be
#' \deqn{Q(\beta|X, y) = \frac{1}{n} L(\beta|X, y) + }{Q(\beta|X, y) =
#' (1/n)*L(\beta|X, y) + P(\beta, \lambda),}\deqn{ P_\lambda(\beta)}{Q(\beta|X,
#' y) = (1/n)*L(\beta|X, y) + P(\beta, \lambda),} where the loss function L is
#' the negative log-likelihood (half the deviance) for the specified outcome
#' distribution (gaussian/binomial/poisson). For more details, refer to the
#' following:
#' \itemize{
#'   \item [Models and loss functions](https://pbreheny.github.io/grpreg/articles/models.html)
#'   \item [Penalties](https://pbreheny.github.io/grpreg/articles/penalties.html)
#' }
#'   
#' For the bi-level selection methods, a locally approximated coordinate
#' descent algorithm is employed.  For the group selection methods, group
#' descent algorithms are employed.
#' 
#' The algorithms employed by \code{grpreg} are stable and generally converge
#' quite rapidly to values close to the solution.  However, especially when p
#' is large compared with n, \code{grpreg} may fail to converge at low values
#' of \code{lambda}, where models are nonidentifiable or nearly singular.
#' Often, this is not the region of the coefficient path that is most
#' interesting.  The default behavior warning the user when convergence
#' criteria are not met may be distracting in these cases, and can be modified
#' with \code{warn} (convergence can always be checked later by inspecting the
#' value of \code{iter}).
#' 
#' If models are not converging, increasing \code{max.iter} may not be the most
#' efficient way to correct this problem.  Consider increasing \code{n.lambda}
#' or \code{lambda.min} in addition to increasing \code{max.iter}.
#' 
#' Although \code{grpreg} allows groups to be unordered and given arbitary
#' names, it is recommended that you specify groups as consecutive integers.
#' The first reason is efficiency: if groups are out of order, \code{X} must be
#' reordered prior to fitting, then this process reversed to return
#' coefficients according to the original order of \code{X}.  This is
#' inefficient if \code{X} is very large.  The second reason is ambiguity with
#' respect to other arguments such as \code{group.multiplier}.  With
#' consecutive integers, \code{group=3} unambiguously denotes the third element
#' of \code{group.multiplier}.
#' 
#' Seemingly unrelated regressions/multitask learning can be carried out using
#' \code{grpreg} by passing a matrix to \code{y}.  In this case, \code{X} will
#' be used in separate regressions for each column of \code{y}, with the
#' coefficients grouped across the responses.  In other words, each column of
#' \code{X} will form a group with m members, where m is the number of columns
#' of \code{y}.  For multiple Gaussian responses, it is recommended to
#' standardize the columns of \code{y} prior to fitting, in order to apply the
#' penalization equally across columns.
#' 
#' \code{grpreg} requires groups to be non-overlapping.
#' 
#' @param X The design matrix, without an intercept.  \code{grpreg}
#' standardizes the data and includes an intercept by default.
#' @param y The response vector, or a matrix in the case of multitask learning
#' (see details).
#' @param group A vector describing the grouping of the coefficients.  For
#' greatest efficiency and least ambiguity (see details), it is best if
#' \code{group} is a factor or vector of consecutive integers, although
#' unordered groups and character vectors are also allowed.  If there are
#' coefficients to be included in the model without being penalized, assign
#' them to group 0 (or \code{"0"}).
#' @param penalty The penalty to be applied to the model.  For group selection,
#' one of \code{grLasso}, \code{grMCP}, or \code{grSCAD}.  For bi-level
#' selection, one of \code{gel} or \code{cMCP}.  See below for details.
#' @param family Either "gaussian" or "binomial", depending on the response.
#' @param nlambda The number of \code{lambda} values.  Default is 100.
#' @param lambda A user supplied sequence of \code{lambda} values.  Typically,
#' this is left unspecified, and the function automatically computes a grid of
#' lambda values that ranges uniformly on the log scale over the relevant range
#' of lambda values.
#' @param lambda.min The smallest value for \code{lambda}, as a fraction of
#' \code{lambda.max}.  Default is .0001 if the number of observations is larger
#' than the number of covariates and .05 otherwise.
#' @param log.lambda Whether compute the grid values of lambda on log scale
#' (default) or linear scale.
#' @param alpha \code{grpreg} allows for both a group penalty and an L2 (ridge)
#' penalty; \code{alpha} controls the proportional weight of the regularization
#' parameters of these two penalties.  The group penalties' regularization
#' parameter is \code{lambda*alpha}, while the regularization parameter of the
#' ridge penalty is \code{lambda*(1-alpha)}.  Default is 1: no ridge penalty.
#' @param eps Convergence threshhold.  The algorithm iterates until the RMSD
#' for the change in linear predictors for each coefficient is less than
#' \code{eps}.  Default is \code{1e-4}.  See details.
#' @param max.iter Maximum number of iterations (total across entire path).
#' Default is 10000.  See details.
#' @param dfmax Limit on the number of parameters allowed to be nonzero.  If
#' this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' @param gmax Limit on the number of groups allowed to have nonzero elements.
#' If this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' @param gamma Tuning parameter of the group or composite MCP/SCAD penalty
#' (see details).  Default is 3 for MCP and 4 for SCAD.
#' @param tau Tuning parameter for the group exponential lasso; defaults to
#' 1/3.
#' @param group.multiplier A vector of values representing multiplicative
#' factors by which each group's penalty is to be multiplied.  Often, this is a
#' function (such as the square root) of the number of predictors in each
#' group.  The default is to use the square root of group size for the group
#' selection methods, and a vector of 1's (i.e., no adjustment for group size)
#' for bi-level selection.
#' @param warn Should the function give a warning if it fails to converge?
#' Default is TRUE.  See details.
#' @param returnX Return the standardized design matrix (and associated group
#' structure information)?  Default is FALSE.
#' @param ... Arguments passed to other functions (such as gBridge).
#' 
#' @return An object with S3 class \code{"grpreg"} containing:
#' \describe{
#' \item{beta}{The fitted matrix of coefficients.  The number of rows is equal
#' to the number of coefficients, and the number of columns is equal to
#' \code{nlambda}.}
#' \item{family}{Same as above.}
#' \item{group}{Same as above.}
#' \item{lambda}{The sequence of \code{lambda} values in the path.}
#' \item{alpha}{Same as above.}
#' \item{deviance}{A vector containing the deviance of the fitted model at each
#' value of `lambda`.}
#' \item{n}{Number of observations.}
#' \item{penalty}{Same as above.}
#' \item{df}{A vector of length `nlambda` containing estimates of effective number of model parameters all the points along the regularization path.  For details on how this is calculated, see Breheny and Huang (2009).}
#' \item{iter}{A vector of length `nlambda` containing the number of iterations until convergence at each value of `lambda`.}
#' \item{group.multiplier}{A named vector containing the multiplicative constant applied to each group's penalty.}
#' }
#' 
#' @author Patrick Breheny
#' 
#' @seealso [cv.grpreg()], as well as [plot.grpreg()] and [select.grpreg()] methods.
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
#' }
#' 
#' @examples
#' # Birthweight data
#' data(Birthwt)
#' X <- Birthwt$X
#' group <- Birthwt$group
#' 
#' # Linear regression
#' y <- Birthwt$bwt
#' fit <- grpreg(X, y, group, penalty="grLasso")
#' plot(fit)
#' fit <- grpreg(X, y, group, penalty="grMCP")
#' plot(fit)
#' fit <- grpreg(X, y, group, penalty="grSCAD")
#' plot(fit)
#' fit <- grpreg(X, y, group, penalty="gel")
#' plot(fit)
#' fit <- grpreg(X, y, group, penalty="cMCP")
#' plot(fit)
#' select(fit, "AIC")
#' 
#' # Logistic regression
#' y <- Birthwt$low
#' fit <- grpreg(X, y, group, penalty="grLasso", family="binomial")
#' plot(fit)
#' fit <- grpreg(X, y, group, penalty="grMCP", family="binomial")
#' plot(fit)
#' fit <- grpreg(X, y, group, penalty="grSCAD", family="binomial")
#' plot(fit)
#' fit <- grpreg(X, y, group, penalty="gel", family="binomial")
#' plot(fit)
#' fit <- grpreg(X, y, group, penalty="cMCP", family="binomial")
#' plot(fit)
#' select(fit, "BIC")
#' 
#' # Multitask learning (simulated example)
#' set.seed(1)
#' n <- 50
#' p <- 10
#' k <- 5
#' X <- matrix(runif(n*p), n, p)
#' y <- matrix(rnorm(n*k, X[,1] + X[,2]), n, k)
#' fit <- grpreg(X, y)
#' # Note that group is set up automatically
#' fit$group
#' plot(fit)
#' @export

grpreg <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),
                   family=c("gaussian","binomial", "poisson"), nlambda=100, lambda,
                   lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, log.lambda = TRUE,
                   alpha=1, eps=1e-4, max.iter=10000, dfmax=p, gmax=length(unique(group)),
                   gamma=ifelse(penalty=="grSCAD", 4, 3), tau=1/3, group.multiplier,
                   warn=TRUE, returnX=FALSE, ...) {

  # Deprecation support / error checking
  if (!missing(penalty)) {
    if (penalty[1]=="gBridge") stop("gBridge has been divorced from the grpreg function; use the gBridge() function instead", call.=FALSE)
    if (penalty[1]=="gMCP") {
      writeLines(strwrap("penalty='gMCP' is deprecated and may not be supported in future versions.  Use penalty='cMCP' instead."))
      penalty <- "cMCP"
    }
    if (penalty[1]=="gLasso") {
      writeLines(strwrap("You have specified penalty='gLasso'; grpreg is assuming you mean group lasso (penalty='grLasso')"))
      penalty <- "grLasso"
    }
  }
  family <- match.arg(family)
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
  
  # Construct XG, yy
  bilevel <- strtrim(penalty, 2) != "gr"
  yy <- newY(y, family)
  XG <- newXG(X, group, group.multiplier, attr(yy, 'm'), bilevel)
  if (nrow(XG$X) != length(yy)) stop("X and y do not have the same number of observations", call.=FALSE)

  # Setup lambda
  if (missing(lambda)) {
    lambda <- setupLambda(XG$X, yy, XG$g, family, penalty, alpha, lambda.min, log.lambda, nlambda, XG$m)
    lam.max <- lambda[1]
    user.lambda <- FALSE
  } else {
    lam.max <- -1
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  # Fit
  n <- length(yy)
  p <- ncol(XG$X)
  K <- as.integer(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  if (K0) {
    lambda[1] <- lambda[1] + 1e-5
    user.lambda <- TRUE
  }
  if (family=="gaussian") {
    if (bilevel) fit <- .Call("lcdfit_gaussian", XG$X, yy, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), XG$m, as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    else fit <- .Call("gdfit_gaussian", XG$X, yy, penalty, K1, K0, lambda, lam.max, alpha, eps, as.integer(max.iter), gamma, XG$m, as.integer(dfmax), as.integer(gmax), as.integer(user.lambda))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    dev <- fit[[2]]
    Eta <- matrix(fit[[3]], nrow=n) + mean(y)
    df <- fit[[4]] + 1 # Intercept
    iter <- fit[[5]]
  } else {
    if (bilevel) fit <- .Call("lcdfit_glm", XG$X, yy, family, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter), XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    else fit <- .Call("gdfit_glm", XG$X, yy, family, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter), gamma, XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    dev <- fit[[3]]
    Eta <- matrix(fit[[4]], nrow=n)
    df <- fit[[5]]
    iter <- fit[[6]]
  }

  # Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  dev <- dev[ind]
  Eta <- Eta[, ind, drop=FALSE]
  if (iter[1] == max.iter) stop("Algorithm failed to converge for any values of lambda.  This indicates a combination of (a) an ill-conditioned feature matrix X and (b) insufficient penalization.  You must fix one or the other for your model to be identifiable.", call.=FALSE)
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda", call.=FALSE)

  # Unstandardize
  if (strtrim(penalty,2)=="gr") b <- unorthogonalize(b, XG$X, XG$g)
  if (XG$reorder) b[-1,] <- b[1+XG$ord.inv,]
  beta <- unstandardize(b, XG)

  # Names
  varnames <- c("(Intercept)", XG$names)
  ncolY <- attr(yy, 'm')
  if (ncolY > 1) {
    beta[2:ncolY,] <- sweep(beta[2:ncolY, , drop=FALSE], 2, beta[1,], FUN="+")
    beta <- array(beta, dim=c(ncolY, nrow(beta)/ncolY, ncol(beta)))
    group <- group[-(1:(ncolY-1))]
    dimnames(beta) <- list(colnames(yy), varnames, round(lambda, digits=4))
  } else {
    dimnames(beta) <- list(varnames, round(lambda, digits=4))
  }
  colnames(Eta) <- round(lambda, digits=4)

  val <- structure(list(beta = beta,
                        family = family,
                        group = factor(group),
                        lambda = lambda,
                        alpha = alpha,
                        deviance = dev,
                        linear.predictors = Eta,
                        n = n,
                        penalty = penalty,
                        df = df,
                        iter = iter,
                        group.multiplier = XG$m),
                   class = "grpreg")
  if (family == 'gaussian') {
    val$y <- yy + attr(yy, 'mean')
  } else {
    val$y <- yy
  }
  if (returnX) val$XG <- XG
  if (expanded) {
    val$meta <- list(knots = knots,
                     boundary = boundary,
                     degree = degree,
                     originalx = originalx,
                     type = type,
                     X = X)
    attr(val, "class") <- c("grpreg", "expanded")
  }
  val
}
