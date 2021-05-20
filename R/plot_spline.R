#' Plot spline curve
#'
#' Plots a spline curve for a single variable using a grpreg object for which an additive model was fit.
#'
#' `plot_spline()` takes a model fit using both the [grpreg()] and [expand_spline()] functions and plots a spline curve for a given variable.
#'
#' @param fit        A \code{grpreg} object. The model must have been fit using a \code{expand_spline} object.
#' @param variable   The name of the variable which will be plotted.
#' @param lambda     Values of the regularization parameter \code{lambda} which will be used for the plot. Each value of lambda will produce a different curve.
#' @param which      Indices of the penalty parameter \code{lambda} which will be used for the plot. If both `lambda` and `which` are specified, `lambda` takes precedence.
#' @param scatter    If \code{TRUE}, a scatter plot of the partial residuals will also be produced. Default is \code{FALSE}. If multiple lambdas are specified, the largest value will be used to calculate the residuals. 
#' @param type       Type of plot to be produced. Default is \code{contrast}.  The following options are supported:
#'                     * If \code{conditional} is selected, the plot returned shows the value of the variable on the x-axis and the change in response on the y-axis, holding all other variables constant at their mean value.
#'                     * If \code{contrast} is selected, the plot returned shows the effect on the expected value of the response by moving the x variable away from the mean on the x-axis. 
#' @param warnings   If \code{FALSE}, warnings will be suppressed. Default is `TRUE`.
#' @param ...        Further arguments to be passed to `plot()`
#'
#' @examples
#' X <- expand_spline(attitude[-1], df = 3)
#' fit <- grpreg(X, attitude$rating, penalty="grLasso")
#' plot_spline(fit, "complaints", which = c(5, 90))

plot_spline <- function(fit, variable, lambda, which = NULL, scatter = FALSE, 
                           type = "contrast", warnings = TRUE, ...){
  if (inherits(fit, "cv.grpreg")) {
    fit <- fit$fit
  }
  if (!inherits(fit, "expanded")) stop("Model must have been fit using a expand_spline object", call. = FALSE)
  if (!inherits(fit, "grpreg")) stop("Model must have been fit using the grpreg function", call. = FALSE)
  if (missing(lambda)) {
    if(missing(which)) stop("Lambda not specified",call. = FALSE)
    lambda <- fit$lambda[which]
  }
  if (length(which(fit$group == variable)) == 0) stop(paste("Cannot find variable", variable), call. = FALSE)
  if (scatter == TRUE && length(lambda) > 1 && warnings == TRUE) warning("Scatter plot represents largest value of lambda imputed")
  
  meta <- fit$meta
  fit$y <- fit$y + attr(fit$y, "mean")
  df <- length(meta$knots[[1]]) + meta$degree
  j <- which(fit$group == variable)
  i <- j[df]/df
  l <- length(lambda)
  p <- ncol(meta$originalx)
  n <- nrow(meta$originalx)
  mat <- fit$XG$X[,j]
  attr(mat, "degree") <- meta$degree
  attr(mat, "knots") <- meta$knots[[i]]
  attr(mat, "Boundary.knots") <- meta$boundary[[i]]
  attr(mat, "intercept") <- FALSE
  attr(mat, "class") <- c(meta$type, "basis", "matrix")
  
  #create sequence and basis and calculate y's and residuals
  min <- meta$boundary[[i]][1]
  max <- meta$boundary[[i]][2]
  newx <- seq(min, max, length.out = 200) 
  if(type == "conditional"){
    xmeans <- matrix(colMeans(meta$originalx), 200, p, byrow = TRUE)
    xmeans[,i] <- newx
    y <- predict(fit, xmeans, lambda = lambda)
    if(scatter == TRUE){
      xmeans <- matrix(colMeans(meta$originalx), n, p, byrow = TRUE)
      xmeans[,i] <- meta$originalx[,i]
      const <- predict(fit, xmeans, lambda = max(lambda))
      betas <- coef.grpreg(fit, lambda = max(lambda))
      parresid <- fit$y - cbind(1, fit$X)%*%betas + const
    }
  } else if(type == "contrast"){
    newxbs <- predict(mat, newx)
    betas <- matrix(coef.grpreg(fit, lambda = lambda), ncol = l)
    newxmean <- predict(mat, mean(meta$originalx[,i]))
    y <- newxbs%*%betas[j+1,] - matrix(newxmean%*%betas[j+1,],200,l, byrow = TRUE)
    if(scatter == TRUE){
      betas <- coef.grpreg(fit, lambda = max(lambda))
      parresid <- fit$y - cbind(1, fit$X)%*%betas + fit$X[,j]%*%betas[j+1] - rep(newxmean%*%betas[j+1], length = length(fit$y))
    }
  } else {
    stop(paste(type, "is not a valid type"), call. = FALSE)
  }
  
  #check for selection
  if(warnings == TRUE){
    notselected <- NULL
    for(q in 1:length(lambda)){
      if (all(coef.grpreg(fit, lambda = lambda[q])[j+1]==0)){
        notselected <- c(notselected, lambda[q])
      }
    }
    if(length(notselected) > 0){
      warning(paste(variable, "was not selected at lambda =", notselected, "\n"))
    }
  }
  
  #plot
  cols <- hcl(h=seq(15, 240, len=l), l=60, c=150, alpha=1)
  if(scatter == TRUE){
    ymax <- max(max(y), max(parresid))
    ymin <- min(min(y), min(parresid))
  } else {
    ymax <- max(y)
    ymin <- min(y)
  }
  plot.args <- list(x=newx, y=y, xlab=variable, ylab = "y", type="l", 
                    col = cols, lty = 1, ylim = c(ymin, ymax))
  new.args <- list(...)
  if (length(new.args)) {
    new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
    plot.args[names(new.plot.args)] <- new.plot.args
  }
  do.call("matplot", plot.args)
  if(scatter == TRUE){
    points(meta$originalx[,i], parresid)
  }
}
