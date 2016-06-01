AUC.cv.grpsurv <- function(obj, ...) {
  if (!("Y" %in% names(obj))) stop("Must run cv.grpsurv with 'returnY=TRUE' in order to calculate AUC")
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is needed for AUC() to work. Please install it.", call. = FALSE)
  }  
  S <- survival::Surv(obj$fit$time, obj$fit$fail)
  res <- apply(obj$Y, 2, survival::survConcordance.fit, y = S)
  num <- res['concordant',] + 0.5*res['tied.risk',] + 0.5*res['tied.time',]
  num/sum(res[1:4,1])
}
AUC <- function(obj, ...) UseMethod("AUC")
