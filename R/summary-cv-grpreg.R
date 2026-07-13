#' Summarizing inferences based on cross-validation
#'
#' Summary method for `cv.grpreg` or `cv.grpsurv` objects
#'
#' @param object A `cv.grpreg` object.
#' @param x A `summary.cv.grpreg` object.
#' @param digits Number of digits past the decimal point to print out. Can be a vector specifying
#'   different display digits for each of the five non-integer printed values.
#' @param ... Further arguments passed to or from other methods.
#'
#' @returns `summary(cvfit)` produces an object with S3 class `"summary.cv.grpreg"`. The class has
#'   its own print method and contains the following list elements:
#'
#' \describe{
#'   \item{penalty}{The penalty used by `grpreg` /`grpsurv`.}
#'   \item{model}{The type of model: `"linear"`, `"logistic"`, `"Poisson"`, `"Cox"`, etc.}
#'   \item{n}{Number of observations}
#'   \item{p}{Number of regression coefficients (not including the intercept).}
#'   \item{min}{The index of `lambda` with the smallest cross-validation error.}
#'   \item{lambda}{The sequence of `lambda` values used by `cv.grpreg` /`cv.grpsurv`.}
#'   \item{cve}{Cross-validation error (deviance).}
#'   \item{r.squared}{
#'     Proportion of variance explained by the model, as estimated by cross-validation.
#'   }
#'   \item{snr}{Signal to noise ratio, as estimated by cross-validation.}
#'   \item{sigma}{For linear regression models, the scale parameter estimate.}
#'   \item{pe}{For logistic regression models, the prediction error (misclassification error).}
#' }
#'
#' @seealso [grpreg()], [cv.grpreg()], [cv.grpsurv()], [plot.cv.grpreg()]
#'
#' @examples
#' # Birthweight data
#' data(Birthwt)
#' X <- Birthwt$X
#' group <- Birthwt$group
#'
#' # Linear regression
#' y <- Birthwt$bwt
#' cvfit <- cv.grpreg(X, y, group)
#' summary(cvfit)
#'
#' # Logistic regression
#' y <- Birthwt$low
#' cvfit <- cv.grpreg(X, y, group, family = "binomial")
#' summary(cvfit)
#'
#' # Cox regression
#' data(Lung)
#' cvfit <- with(Lung, cv.grpsurv(X, y, group))
#' summary(cvfit)
#'
#' @export
summary.cv.grpreg <- function(object, ...) {
  S <- pmax(object$null.dev - object$cve, 0)
  if (!inherits(object, "cv.grpsurv") && object$fit$family == "gaussian") {
    rsq <- pmin(pmax(1 - object$cve / object$null.dev, 0), 1)
  } else {
    rsq <- pmin(pmax(1 - exp(object$cve - object$null.dev), 0), 1)
  }
  snr <- S / object$cve
  nvars <- predict(object$fit, type = "nvars")
  ngroups <- predict(object$fit, type = "ngroups")
  if (inherits(object, "cv.grpsurv")) {
    model <- "Cox"
  } else {
    model <- switch(
      object$fit$family,
      gaussian = "linear",
      binomial = "logistic",
      poisson = "Poisson"
    )
  }
  d <- dim(object$fit$beta)
  if (length(d) == 3) {
    p <- d[2] - 1
  } else {
    if (model == "Cox") {
      p <- d[1]
    } else {
      p <- d[1] - 1
    }
  }
  val <- list(
    penalty = object$fit$penalty,
    model = model,
    n = object$fit$n,
    p = p,
    min = object$min,
    lambda = object$lambda,
    cve = object$cve,
    r.squared = rsq,
    snr = snr,
    nvars = nvars,
    ngroups = ngroups,
    d = d
  )
  if (!inherits(object, "cv.grpsurv") && object$fit$family == "gaussian") {
    val$sigma <- sqrt(object$cve)
  }
  if (!inherits(object, "cv.grpsurv") && object$fit$family == "binomial") {
    val$pe <- object$pe
  }
  structure(val, class = "summary.cv.grpreg")
}

#' @rdname summary.cv.grpreg
#' @export
print.summary.cv.grpreg <- function(x, digits, ...) {
  digits <- if (missing(digits)) digits <- c(2, 4, 2, 2, 3) else rep(digits, length.out = 5)
  if (length(x$d) == 3) {
    cat(
      x$penalty,
      "-penalized multivariate ",
      x$model,
      " regression with m=",
      x$d[1],
      ", n=",
      x$n / x$d[1],
      ", p=",
      x$p,
      "\n",
      sep = ""
    )
  } else {
    cat(x$penalty, "-penalized ", x$model, " regression with n=", x$n, ", p=", x$p, "\n", sep = "")
  }
  cat(
    "At minimum cross-validation error (lambda=",
    formatC(x$lambda[x$min], digits[2], format = "f"),
    "):\n",
    sep = ""
  )
  cat("-------------------------------------------------\n")
  cat("  Nonzero coefficients: ", x$nvars[x$min], "\n", sep = "")
  cat("  Nonzero groups: ", x$ngroups[x$min], "\n", sep = "")
  cat("  Cross-validation error of ", formatC(min(x$cve), digits[1], format = "f"), "\n", sep = "")
  cat("  Maximum R-squared: ", formatC(max(x$r.squared), digits[3], format = "f"), "\n", sep = "")
  cat(
    "  Maximum signal-to-noise ratio: ",
    formatC(max(x$snr), digits[4], format = "f"),
    "\n",
    sep = ""
  )
  if (x$model == "logistic") {
    cat(
      "  Prediction error at lambda.min: ",
      formatC(x$pe[x$min], digits[5], format = "f"),
      "\n",
      sep = ""
    )
  }
  if (x$model == "linear") {
    cat(
      "  Scale estimate (sigma) at lambda.min: ",
      formatC(sqrt(x$cve[x$min]), digits[5], format = "f"),
      "\n",
      sep = ""
    )
  }
}
