#' Fit a group bridge regression path
#' 
#' Fit regularization paths for linear and logistic group bridge-penalized
#' regression models over a grid of values for the regularization parameter
#' lambda.
#' 
#' This method fits the group bridge method of Huang et al. (2009).  Unlike the
#' penalties in \code{grpreg}, the group bridge is not differentiable at zero;
#' because of this, a number of changes must be made to the algorithm, which is
#' why it has its own function.  Most notably, the method is unable to start at
#' \code{lambda.max}; it must start at \code{lambda.min} and proceed in the
#' opposite direction.
#' 
#' In other respects, the usage and behavior of the function is similar to the
#' rest of the \code{grpreg} package.
#' 
#' @param X The design matrix, as in \code{grpreg}.
#' @param y The response vector (or matrix), as in \code{grpreg}.
#' @param group The grouping vector, as in \code{grpreg}.
#' @param family Either "gaussian" or "binomial", depending on the response.
#' @param nlambda The number of \code{lambda} values, as in \code{grpreg}.
#' @param lambda A user supplied sequence of `lambda values, as in `grpreg()`.
#' @param lambda.min The smallest value for \code{lambda}, as in \code{grpreg}.
#' @param lambda.max The maximum value for \code{lambda}.  Unlike the penalties
#' in \code{grpreg}, it is not possible to solve for \code{lambda.max} directly
#' with group bridge models.  Thus, it must be specified by the user.  If it is
#' not specified, \code{gBridge} will attempt to guess \code{lambda.max}, but
#' this is not particularly accurate.
#' @param alpha Tuning parameter for the balance between the group penalty and
#' the L2 penalty, as in \code{grpreg}.
#' @param eps Convergence threshhold, as in \code{grpreg}.
#' @param delta The group bridge penalty is not differentiable at zero, and
#' requires a small number \code{delta} to bound it away from zero.  There is
#' typically no need to change this value.
#' @param max.iter Maximum number of iterations, as in \code{grpreg}.
#' @param gamma Tuning parameter of the group bridge penalty (the exponent to
#' which the L1 norm of the coefficients in the group are raised).  Default is
#' 0.5, the square root.
#' @param group.multiplier The multiplicative factor by which each group's
#' penalty is to be multiplied, as in \code{grpreg}.
#' @param warn Should the function give a warning if it fails to converge?  As
#' in \code{grpreg}.
#' @param returnX Return the standardized design matrix (and associated group
#' structure information)?  Default is FALSE.
#' @param ... Not used.
#' 
#' @return An object with S3 class \code{"grpreg"}, as in \code{grpreg}.
#' 
#' @seealso [grpreg()]
#' 
#' @references
#' \itemize{
#' \item Huang J, Ma S, Xie H, and Zhang C. (2009) A group bridge approach for
#' variable selection. *Biometrika*, **96**: 339-355. \doi{10.1093/biomet/asp020}
#' 
#' \item Breheny P and Huang J. (2009) Penalized methods for bi-level variable
#' selection. *Statistics and its interface*, **2**: 369-380.
#' \doi{10.4310/sii.2009.v2.n3.a10}
#' }
#' 
#' @examples
#' data(Birthwt)
#' X <- Birthwt$X
#' group <- Birthwt$group
#' 
#' ## Linear regression
#' y <- Birthwt$bwt
#' fit <- gBridge(X, y, group, lambda.max=0.08)
#' plot(fit)
#' select(fit)$beta
#' 
#' ## Logistic regression
#' y <- Birthwt$low
#' fit <- gBridge(X, y, group, family="binomial", lambda.max=0.17)
#' plot(fit)
#' select(fit)$beta
#' @export

gBridge <- function(X, y, group=1:ncol(X), family=c("gaussian", "binomial", "poisson"), nlambda=100, lambda,
                    lambda.min={if (nrow(X) > ncol(X)) .001 else .05}, lambda.max, alpha=1, eps=.001, delta=1e-7,
                    max.iter=10000, gamma=0.5, group.multiplier, warn=TRUE, returnX=FALSE, ...) {
  # Error checking
  family <- match.arg(family)
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0, 1]", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to gBridge", call.=FALSE)
  if (length(group)!=ncol(X)) stop("group does not match X", call.=FALSE)
  if (delta <= 0) stop("Delta must be a positive number", call.=FALSE)

  # Construct XG, yy
  yy <- newY(y, family)
  m <- attr(yy, "m")
  XG <- newXG(X, group, group.multiplier, m, TRUE)
  if (nrow(XG$X) != length(yy)) stop("X and y do not have the same number of observations", call.=FALSE)

  # Set up lambda
  if (missing(lambda)) {
    lambda <- setupLambda.gBridge(XG$X, yy, XG$g, family, alpha, lambda.min, lambda.max, nlambda, gamma, XG$m)
  } else {
    nlambda <- length(lambda)
  }

  # Fit
  n <- length(yy)
  p <- ncol(XG$X)
  K <- as.integer(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  if (family=="gaussian") {
    fit <- .Call("lcdfit_gaussian", XG$X, yy, "gBridge", K1, K0, lambda, alpha, eps, delta, gamma, 0, as.integer(max.iter), as.double(XG$m), as.integer(p), as.integer(max(XG$g)), as.integer(TRUE))
    b <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    loss <- fit[[2]]
    Eta <- matrix(fit[[3]], nrow=n) + mean(y)
    df <- fit[[4]] + 1 # Intercept
    iter <- fit[[5]]
  } else {
    fit <- .Call("lcdfit_glm", XG$X, yy, family, "gBridge", K1, K0, lambda, alpha, eps, delta, gamma, 0, as.integer(max.iter), as.double(XG$m), as.integer(p), as.integer(max(XG$g)), as.integer(warn), as.integer(TRUE))
    b <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    loss <- fit[[3]]
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
  loss <- loss[ind]
  if (iter[1] == max.iter) stop("Algorithm failed to converge for any values of lambda.  This indicates a combination of (a) an ill-conditioned feature matrix X and (b) insufficient penalization.  You must fix one or the other for your model to be identifiable.", call.=FALSE)
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda", call.=FALSE)

  # Unstandardize
  if (XG$reorder) b[-1,] <- b[1+XG$ord.inv,]
  beta <- unstandardize(b, XG)

  # Names
  varnames <- c("(Intercept)", XG$names)
  if (m > 1) {
    beta[2:m,] <- sweep(beta[2:m, , drop=FALSE], 2, beta[1,], FUN="+")
    beta <- array(beta, dim=c(m, nrow(beta)/m, ncol(beta)))
    group <- group[-(1:(m-1))]
    dimnames(beta) <- list(colnames(yy), varnames, round(lambda, digits=4))
  } else {
    dimnames(beta) <- list(varnames, round(lambda, digits=4))
  }

  val <- structure(list(beta = beta,
                        family = family,
                        group = group,
                        lambda = lambda,
                        alpha = alpha,
                        deviance = 2 * loss,
                        n = n,
                        penalty = "gBridge",
                        df = df,
                        iter = iter,
                        group.multiplier = XG$m),
                   class = "grpreg")
  if (returnX) {
    val$XG = XG
    val$y = yy
  } else if (family=="poisson") {
    val$y <- y
  }
  val
}
