#' Marginal false discovery rates
#'
#' Estimates the marginal false discovery rate (mFDR) of a group penalized
#' regression model.
#'
#' The function estimates the marginal false discovery rate (mFDR) for groups in
#' a group lasso or group MCP penalized regression model. The estimate tends to
#' be accurate in most settings, but will be somewhat conservative if predictors
#' are highly correlated.
#'
#' @param fit A `grpreg` or `grpsurv` object.
#' @param X The model matrix corresponding to `fit`. This is not necessary
#'   for linear regression, but in logistic and Cox regression, the mFDR depends
#'   on X. It is not necessary to supply `X` if it is already contained in
#'   `fit`; i.e., if `ncvreg`/`ncvsurv` was run with `returnX = TRUE`.
#'
#' @returns An object with S3 class `mfdr` inheriting from `data.frame`, containing:
#' \describe{
#'   \item{ef}{The number of variables selected at each value of `lambda`,
#'     averaged over the permutation fits.}
#'   \item{s}{The actual number of selected variables for the non-permuted data.}
#'   \item{mfdr}{The estimated marginal false discovery rate (`ef/s`).}
#' }
#'
#' @seealso [grpreg()], [grpsurv()]
#'
#' @examples
#' # Birthweight data ---------------------------
#' data(Birthwt)
#' x <- Birthwt$X
#' group <- Birthwt$group
#'
#' # Linear regression --------------------------
#' y <- Birthwt$bwt
#' fit <- grpreg(x, y, group)
#' obj <- mfdr(fit)
#' head(obj)
#'
#' # Logistic regression ------------------------------
#' y <- Birthwt$low
#' fit <- grpreg(x, y, group, penalty="grMCP", family="binomial", returnX = TRUE)
#' obj <- mfdr(fit)
#' # If returnX is not TRUE, user must supply X
#' fit <- grpreg(x, y, group, penalty="grMCP", family="binomial")
#' obj <- mfdr(fit, x)
#' head(obj)
#'
#' # Cox regression -----------------------------------
#' data(Lung)
#' x_lung <- Lung$X
#' y_lung <- Lung$y
#' g_lung <- Lung$group
#' fit <- grpsurv(x_lung, y_lung, g_lung, penalty = "grSCAD")
#' obj <- mfdr(fit, x_lung)
#' head(obj)
#' @export mfdr

mfdr <- function(fit, X) {
  # Initial checks
  if (!inherits(fit, "grpreg")) stop('"fit" must be an grpreg object', call. = FALSE)
  if (!(fit$penalty %in% c("grLasso", "grMCP", "grSCAD")))
    stop('mfdr() is only avaiable for "grLasso", "grMCP", and "grSCAD" penalties', call. = FALSE)
  if (!missing(X)) {
    if (inherits(fit, "grpsurv")) {
      m <- 1
    } else {
      m <- attr(fit$y, "m")
    }
    fit$XG <- newXG(X, fit$group, fit$group.multiplier, m, FALSE)
  }
  needs_x <- inherits(fit, "grpsurv") || fit$family == "binomial"
  if (needs_x && !("XG" %in% names(fit))) {
    stop(
      paste(
        "For GLM/Cox models, you must either:",
        "- Supply X or",
        "- Specify 'returnX = TRUE' in the call to grpreg()",
        "",
        "This is required to calculate mFDR.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  ## Setup
  s0 <- sum(fit$group.multiplier == 0)
  s <- predict(fit, type = "ngroups") - s0
  fit$gl <- as.vector(table(fit$group))

  # Call C functions
  if (inherits(fit, "grpsurv")) {
    fit$XX <- fit$XG$X
    ef <- .Call("mfdr_cox", fit)
  } else if (fit$family == "binomial") {
    fit$P <- binomial()$linkinv(fit$linear.predictors)
    fit$XX <- fit$XG$X
    ef <- .Call("mfdr_binomial", fit)
  } else {
    ef <- .Call("mfdr_gaussian", fit)
  }

  # Calculate rate, return
  ef <- pmin(ef - s0, s)
  mfdr <- ef / s
  mfdr[s == 0] <- 0
  out <- data.frame(ef = ef, s = s, mfdr = mfdr)
  rownames(out) <- lamNames(fit$lambda)
  structure(out, class = c("mfdr", "data.frame"))
}
