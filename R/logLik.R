#' logLik method for grpreg
#'
#' Calculates the log likelihood and degrees of freedom for a fitted grpreg object.
#'
#' Exists mainly for use with [stats::AIC()] and [stats::BIC()].
#'
#' @param object A fitted `grpreg` or `grpsurv` object, as obtained from [grpreg()] or [grpsurv()]
#' @param df.method How should effective model parameters be calculated? One of: `"active"`, which
#'   counts the number of nonzero coefficients; or `"default"`, which uses the calculated `df`
#'   returned by `grpreg`. Default is `"default"`.
#' @param REML Use restricted MLE for estimation of the scale parameter in a gaussian model? Default
#'   is FALSE.
#' @param ... For S3 method compatibility.
#'
#' @returns
#' An object with S3 class `logLik`. This is a vector of numbers (the log-likelihood at each value
#' of the regularization parameter) with the following attributes:
#'
#' \describe{
#'   \item{df}{Degrees of freedom}
#'   \item{nobs}{Number of observations}
#' }
#'
#' The 'print' method for 'logLik' objects is not intended to handle vectors; consequently, the
#' value of the function does not necessarily display correctly. However, it works with 'AIC' and
#' 'BIC' without any glitches and returns the expected vectorized output.
#'
#' @seealso [stats::logLik()]
#'
#' @examples
#' data(Birthwt)
#' X <- Birthwt$X
#' y <- Birthwt$bwt
#' group <- Birthwt$group
#' fit <- grpreg(X, y, group, penalty = "cMCP")
#' logLik(fit) ## Display is glitchy for vectors
#' AIC(fit)
#' BIC(fit)
#'
#' @rdname logLik.grpreg
#' @export
logLik.grpreg <- function(object, df.method = c("default", "active"), REML = FALSE, ...) {
  df.method <- match.arg(df.method)
  n <- as.integer(object$n)
  df <- if (df.method == "active") apply(coef(object) != 0, 2, sum) else object$df
  if (object$family == "gaussian") {
    rdf <- if (REML) n - df else n
    RSS <- object$deviance
    l <- -n / 2 * (log(2 * pi) + log(RSS) - log(rdf)) - rdf / 2
    df <- df + 1
  } else if (object$family == "poisson") {
    y <- object$y
    ind <- y != 0
    l <- -object$deviance / 2 + sum(y[ind] * log(y[ind])) - sum(y) - sum(lfactorial(y))
  } else {
    l <- -object$deviance / 2
  }
  structure(l, df = df, nobs = n, class = "logLik")
}

#' @rdname logLik.grpreg
#' @export
logLik.grpsurv <- function(object, df.method = c("default", "active"), ...) {
  df.method <- match.arg(df.method)
  n <- as.integer(object$n)
  df <- if (df.method == "active") apply(coef(object) != 0, 2, sum) else object$df
  structure(-object$deviance / 2, df = df, nobs = n, class = "logLik")
}
