summary.cv.grpsurv <- function(object, ...) {
  S <- pmax(object$null.dev - object$cve, 0)
  rsq <- S/object$null.dev
  snr <- S/object$cve
  nvars <- predict(object$fit, type="nvars")
  ngroups <- predict(object$fit, type="ngroups")
  model <- "Cox"
  val <- list(penalty=object$fit$penalty,
              model=model,
              n=object$fit$n,
              p=nrow(object$fit$beta),
              min=object$min,
              lambda=object$lambda,
              cve=object$cve,
              r.squared=rsq,
              snr=snr,
              nvars=nvars,
              ngroups=ngroups,
              d=2)
  structure(val, class="summary.cv.grpreg")
}
