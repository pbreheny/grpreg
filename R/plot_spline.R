#' Plot spline curve for a fitted additive model
#'
#' Plots a spline curve for a single variable using a `grpreg` or `cv.grpreg` object for which an additive model was fit.
#'
#' `plot_spline()` takes a model fit using both the [grpreg()] and [expand_spline()] functions and plots a spline curve for a given variable.
#'
#' @param fit         A `grpreg` object. The model must have been fit using a `expand_spline` object.
#' @param variable    The name of the variable which will be plotted (character).
#' @param lambda      Values of the regularization parameter `lambda` which will be used for the plot. If a vector is passed, a curve will be drawn for each value of lambda (numeric vector; if a `cv.grpreg` object is passed, the `lambda` value minimizing cross-validation error will be used as a default; otherwise, there is no default value)
#' @param which       Index of penalty parameter `lambda` which will be used for the plot. If both `lambda` and `which` are specified, `lambda` takes precedence (integer vector).
#' @param partial     If `TRUE`, a scatter plot of the partial residuals is superimposed on the curve (logical; default = `FALSE`). If multiple lambdas are specified, the largest value is used to calculate the residuals.
#' @param type        Type of plot to be produced (default = `"contrast"`). The following options are supported:
#'   * If `"conditional"`, the plot returned shows the value of the variable on the x-axis and the change in linear predictor on the y-axis, holding all other variables constant at their mean value.
#'   * If `"contrast"`, the plot returned shows the effect on the linear predictor by moving the x variable away from its mean.
#' @param warnings    If `FALSE`, warnings will be suppressed (default = `TRUE`).
#' @param points.par  List of parameters (see [par()] to pass to [points()] when `partial=TRUE`.
#' @param ...         Further arguments to be passed to `plot()`. Note that these arguments also control the appearance of the lines.
#'
#' @examples
#' Data <- gen_nonlinear_data(n=1000)
#' X <- expand_spline(Data$X)
#' fit <- grpreg(X, Data$y)
#' plot_spline(fit, "V02", lambda = 0.03)
#' plot_spline(fit, "V02", which = c(10, 90))
#' plot_spline(fit, "V02", lambda = 0.03, partial=TRUE)
#' plot_spline(fit, "V02", lambda = 0.03, partial=TRUE, type='conditional')
#' plot_spline(fit, "V02", lambda = 0.03, partial=TRUE, lwd=6, col='yellow',
#'             points.par=list(pch=9, col='blue'))
#' 
#' op <- par(mfrow=c(3,2), mar=c(4.5, 4.5, 0.25, 0.25))
#' for (i in 1:6) plot_spline(fit, sprintf("V%02d", i), lambda = 0.03, partial=TRUE)
#' par(op)
#' 
#' cvfit <- cv.grpreg(X, Data$y)
#' plot_spline(cvfit, "V02")
#' plot_spline(cvfit, "V02", partial=TRUE)

plot_spline <- function(fit, variable, lambda, which = NULL, partial = FALSE, 
                           type = "contrast", warnings = TRUE, points.par = NULL, ...){
  if (inherits(fit, "cv.grpreg")) {
    if (missing(lambda) & missing(which)) lambda <- fit$lambda.min
    fit <- fit$fit
  }
  if (!inherits(fit, "expanded")) stop("Model must have been fit using a expand_spline object", call. = FALSE)
  if (!inherits(fit, "grpreg")) stop("Model must have been fit using the grpreg function", call. = FALSE)
  if (missing(lambda)) {
    if (missing(which)) stop("lambda not specified", call. = FALSE)
    if (!all(which %in% 1:length(fit$lambda))) stop("'which' must be an integer between 1 and length(fit$lambda)", call.=FALSE)
    lambda <- fit$lambda[which]
  }
  if (length(which(fit$group == variable)) == 0) stop(paste("Cannot find variable", variable), call. = FALSE)
  if (partial == TRUE && length(lambda) > 1 && warnings == TRUE) warning("Scatter plot represents largest value of lambda imputed")
  
  meta <- fit$meta
  if (meta$type == 'bs') {
    df <- length(meta$knots[[1]]) + meta$degree
  } else if (meta$type == 'ns') {
    df <- length(meta$knots[[1]]) + 1
  }
  j <- which(fit$group == variable)
  i <- j[df]/df
  l <- length(lambda)
  p <- ncol(meta$originalx)
  n <- nrow(meta$originalx)
  
  #create sequence and basis and calculate y's and residuals
  min <- meta$boundary[[i]][1]
  max <- meta$boundary[[i]][2]
  newx <- seq(min, max, length.out = 200)
  if (type == "conditional"){
    xmeans <- matrix(colMeans(meta$originalx), 200, p, byrow = TRUE)
    xmeans[,i] <- newx
    y <- predict(fit, xmeans, lambda = lambda)
    if (partial == TRUE) {
      xmeans <- matrix(colMeans(meta$originalx), n, p, byrow = TRUE)
      xmeans[,i] <- meta$originalx[,i]
      const <- predict(fit, xmeans, lambda = max(lambda))
      r <- residuals(fit, lambda=max(lambda))
      if (inherits(fit, 'grpsurv')) r <- r[fit$order]
      parresid <- residuals(fit, lambda=max(lambda)) + const
    }
  } else if (type == "contrast") {
    mat <- fit$meta$X[,j]
    attr(mat, "degree") <- meta$degree
    attr(mat, "knots") <- meta$knots[[i]]
    attr(mat, "Boundary.knots") <- meta$boundary[[i]]
    attr(mat, "intercept") <- FALSE
    attr(mat, "class") <- c(meta$type, "basis", "matrix")
    newxbs <- predict(mat, newx)
    betas <- matrix(coef.grpreg(fit, lambda = lambda), ncol = l)
    if (!inherits(fit, 'grpsurv')) betas <- betas[-1,,drop=FALSE]
    newxmean <- predict(mat, mean(meta$originalx[,i]))
    y <- newxbs%*%betas[j,] - matrix(newxmean%*%betas[j,],200,l, byrow = TRUE)
    if (partial == TRUE) {
      betas <- coef.grpreg(fit, lambda = max(lambda))
      if (!inherits(fit, 'grpsurv')) betas <- betas[-1]
      r <- residuals(fit, lambda=max(lambda))
      if (inherits(fit, 'grpsurv')) r <- r[fit$order]
      offset <- rep(newxmean %*% betas[j], length = fit$n)
      parresid <- r + fit$meta$X[,j] %*% betas[j] - offset
    }
  } else {
    stop(paste(type, "is not a valid type"), call. = FALSE)
  }
  
  #check for selection
  if (warnings == TRUE){
    notselected <- NULL
    for(q in 1:length(lambda)) {
      if (all(coef.grpreg(fit, lambda = lambda[q])[j+1]==0)){
        notselected <- c(notselected, lambda[q])
      }
    }
    if (length(notselected) > 0){
      warning(paste(variable, "was not selected at lambda =", notselected, "\n"))
    }
  }
  
  #plot
  cols <- hcl(h=seq(15, 240, len=l), l=60, c=150, alpha=1)
  if (partial == TRUE){
    ymax <- max(max(y), max(parresid))
    ymin <- min(min(y), min(parresid))
  } else {
    ymax <- max(y)
    ymin <- min(y)
  }
  plot.args <- list(x=newx, y=y, xlab=variable, ylab = "y", type="l", las = 1, bty='n',
                    col = cols, lty = 1, ylim = c(ymin, ymax), lwd=2)
  new.args <- list(...)
  if (length(new.args)) {
    new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
    plot.args[names(new.plot.args)] <- new.plot.args
  }
  do.call("matplot", plot.args)
  if (partial == TRUE) {
    points.args <- list(x=meta$originalx[,i], y=parresid, pch=19, cex=0.8, col='gray')
    if (length(points.par)) {
      new.points.args <- points.par[names(points.par) %in% c(names(par()), names(formals(points)))]
      points.args[names(new.points.args)] <- new.points.args
    }
    do.call("points", points.args)
  }
  # Re-plot line so that it stays on top
  line.args <- list(x=newx, y=y, col = cols, lty = 1, lwd=2)
  if (length(new.args)) {
    new.line.args <- new.args[names(new.args) %in% c(names(par()), names(formals(matlines)))]
    line.args[names(new.line.args)] <- new.line.args
  }
  do.call("matlines", line.args)
}
