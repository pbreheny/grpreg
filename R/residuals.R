#' Extract residuals from a grpreg or grpsurv fit
#' 
#' Currently, only residuals on the linear predictor scale are supported.
#' 
#' @param object   Object of class `grpreg` or `grpsurv`.
#' @param lambda   Values of the regularization parameter at which residuals are requested. For values of lambda not in the sequence of fitted models, linear interpolation is used.
#' @param which    Indices of the penalty parameter at which residuals are requested. Default: all indices. If lambda is specified, this will override `which`.
#' @param ...      Not used.
#' 
#' @examples
#' data(Birthwt)
#' X <- Birthwt$X
#' y <- Birthwt$bwt
#' group <- Birthwt$group
#' fit <- grpreg(X, y, group, returnX=TRUE)
#' residuals(fit)
#' residuals(fit, lambda=0.1)
#' @export

residuals.grpreg <- function(object, lambda, which=1:length(object$lambda), ...) {
  NULL
}
