#' @rdname select
#' @export
select <- function(obj,...) UseMethod("select")

#' Select an value of lambda along a grpreg path
#' 
#' Selects a point along the regularization path of a fitted grpreg object
#' according to the AIC, BIC, or GCV criteria.
#' 
#' The criteria are defined as follows, where \eqn{L}{L} is the deviance (i.e,
#' -2 times the log-likelihood), \eqn{\nu}{df} is the degrees of freedom, and
#' \eqn{n}{n} is the sample size:
#' 
#' \deqn{AIC = L + 2\nu}{AIC = L + 2*df} \deqn{BIC = L + \log(n)\nu}{BIC = L +
#' log(n)*df} \deqn{GCV = \frac{L}{(1-\nu/n)^2}}{GCV= L/((1-df/n)^2)}
#' \deqn{AICc = AIC + 2\frac{\nu(\nu+1)}{n-\nu-1}}{AICc = AIC +
#' 2*df*(df+1)/(n-df-1)} \deqn{EBIC = BIC + 2 \log{p \choose \nu}}{EBIC = BIC +
#' 2*log(p choose df)}
#' 
#' @rdname select
#' 
#' @param obj A fitted grpreg object.
#' @param criterion The criterion by which to select the regularization
#' parameter.  One of \code{"AIC"}, \code{"BIC"}, \code{"GCV"}, \code{"AICc"},
#' or \code{"EBIC"}; default is \code{"BIC"}.
#' @param df.method How should effective model parameters be calculated?  One
#' of: \code{"active"}, which counts the number of nonzero coefficients; or
#' \code{"default"}, which uses the calculated \code{df} returned by
#' \code{grpreg}.  Default is \code{"default"}.
#' @param smooth Applies a smoother to the information criteria before
#' selecting the optimal value.
#' @param \dots For S3 method compatibility.
#' 
#' @return A list containing:
#' \describe{
#' \item{lambda}{The selected value of the regularization parameter, `lambda`.}
#' \item{beta}{The vector of coefficients at the chosen value of `lambda`.}
#' \item{df}{The effective number of model parameters at the chosen value of `lambda`.}
#' \item{IC}{A vector of the calculated model selection criteria for each point on the regularization path.}
#' }
#' 
#' @seealso [grpreg()]
#' 
#' @examples
#' data(Birthwt)
#' X <- Birthwt$X
#' y <- Birthwt$bwt
#' group <- Birthwt$group
#' fit <- grpreg(X, y, group, penalty="grLasso")
#' select(fit)
#' select(fit,crit="AIC",df="active")
#' plot(fit)
#' abline(v=select(fit)$lambda)
#' par(mfrow=c(1,3))
#' l <- fit$lambda
#' xlim <- rev(range(l))
#' plot(l, select(fit)$IC, xlim=xlim, pch=19, type="o", ylab="BIC")
#' plot(l, select(fit,"AIC")$IC, xlim=xlim, pch=19, type="o",ylab="AIC")
#' plot(l, select(fit,"GCV")$IC, xlim=xlim, pch=19, type="o",ylab="GCV")
#' @export

select.grpreg <- function(obj, criterion=c("BIC","AIC","GCV","AICc","EBIC"), df.method=c("default","active"), smooth=FALSE, ...) {
  criterion <- match.arg(criterion)
  df.method <- match.arg(df.method)
  ll <- logLik(obj, df.method=df.method, ...)
  df <- as.double(attr(ll,"df"))
  d <- dim(obj$beta)
  p <- if (length(d)==2) d[1] - 1 else d[2] - 1
  j <- if(obj$family=="gaussian") df - 2 else df - 1
  
  IC <- switch(criterion,
               AIC = AIC(ll),
               BIC = BIC(ll),
               GCV = (1/obj$n) * (-2) * as.double(ll) / (1-df/obj$n)^2,
               AICc = AIC(ll) + 2*df*(df+1)/(obj$n-df-1),
               EBIC = BIC(ll) + 2*(lgamma(p+1) - lgamma(j+1) - lgamma(p-j+1)))
  n.l <- length(obj$lambda)
  if (smooth & (n.l < 4)) {
    smooth <- FALSE
    warning("Need at least 4 points to use smooth=TRUE", call.=FALSE)
  }
  if (smooth) {
    fit.ss <- smooth.spline(IC[is.finite(IC)])
    d <- diff(fit.ss$y)
    if (all(d<0)) i <- n.l
    else i <- min(which(d>0))-1
    if (i==0) i <- 1
  } else i <- which.min(IC)
  
  if (min(obj$lambda) == obj$lambda[i]) {
    warning(paste("minimum lambda selected for", obj$penalty), call.=FALSE)
  } else if ((max(obj$lambda) == obj$lambda[i]) & obj$penalty=="gBridge") {
    warning("maximum lambda selected", call.=FALSE)
  }
  return(list(beta=obj$beta[,i],
              lambda=obj$lambda[i],
              df=df[i],
              IC=IC))
}

