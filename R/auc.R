#' Calculates AUC for cv.grpsurv objects
#' 
#' Calculates the cross-validated AUC (concordance) from a "cv.grpsurv" object.
#' 
#' The area under the curve (AUC), or equivalently, the concordance statistic
#' (C), is calculated according to the procedure described in van Houwelingen
#' and Putter (2011). The function calls `survival::concordancefit()`, except
#' cross-validated linear predictors are used to guard against overfitting.
#' Thus, the values returned by `AUC.cv.grpsurv()` will be lower than those you
#' would obtain with `concordancefit()` if you fit the full (unpenalized) model.
#' 
#' @aliases AUC
#' 
#' @param obj     A `cv.grpsurv` object. You must run `cv.grpsurv()` with the option `returnY=TRUE` in order for \code{AUC} to work.
#' @param \dots   For S3 method compatibility.
#' 
#' @seealso [cv.grpsurv()], [survival::survConcordance()]
#' 
#' @references van Houwelingen H, Putter H (2011). *Dynamic Prediction in Clinical Survival Analysis*. CRC Press.
#' 
#' @examples
#' \dontshow{set.seed(1)}
#' data(Lung)
#' X <- Lung$X
#' y <- Lung$y
#' group <- Lung$group
#' 
#' cvfit <- cv.grpsurv(X, y, group, returnY=TRUE)
#' head(AUC(cvfit))
#' ll <- log(cvfit$fit$lambda)
#' plot(ll, AUC(cvfit), xlim=rev(range(ll)), lwd=3, type='l',
#'      xlab=expression(log(lambda)), ylab='AUC', las=1)
#' @export

AUC.cv.grpsurv <- function(obj, ...) {
  if (!("Y" %in% names(obj))) stop("Must run cv.grpsurv with 'returnY=TRUE' in order to calculate AUC", call.=FALSE)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is needed for AUC() to work. Please install it.", call. = FALSE)
  }
  if (utils::packageVersion("survival") < "3.2.10") stop("AUC.cv.grpsurv requires version 3.2.10 of 'survival' package or higher", call.=FALSE)
  S <- survival::Surv(obj$fit$time, obj$fit$fail)
  apply(obj$Y, 2, concord, y = S)
}

#' @rdname AUC.cv.grpsurv
#' @export

AUC <- function(obj, ...) UseMethod("AUC")

concord <- function(x, y) {
  survival::concordancefit(y, -x)$concordance
}
