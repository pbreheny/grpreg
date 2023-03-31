#' @rdname predict.grpreg
#' @export

predict.cv.grpreg <- function(object, X, lambda=object$lambda.min, which=object$min, type=c("link", "response", "class", "coefficients", "vars", "groups", "nvars", "ngroups", "norm"), ...) {
  type <- match.arg(type)
  if (inherits(object, 'cv.grpsurv')) {
    return(predict.grpsurv(object$fit, X=X, lambda=lambda, which=which, type=type, ...))
  } else {
    return(predict.grpreg(object$fit, X=X, lambda=lambda, which=which, type=type, ...))
  }
}

#' @rdname predict.grpreg
#' @export

coef.cv.grpreg <- function(object, lambda=object$lambda.min, which=object$min, ...) {
  coef.grpreg(object$fit, lambda=lambda, which=which, ...)
}
