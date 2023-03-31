#' Expand feature matrix using basis splines
#'
#' Performs a basis expansion for many features at once, returning output that is compatible
#' for use with the `grpreg()` function. Returns an expanded matrix along with a vector
#' that describes its grouping.
#'
#' `expand_spline()` uses the function [splines::bs()] or [splines::ns()] to generate a basis
#' matrix for each column of `x`. These matrices represent the spline basis for piecewise
#' polynomials with specified degree evaluated separately for each original column of `x`.
#' These matrices are then column-bound to form a single grouped matrix of derived features. A vector
#' that describes the grouping present in the resulting matrix is also generated. The resulting 
#' object can be passed to [grpreg()].
#' 
#' This methodology was originally proposed by Ravikumar et al. (2009), who named it SPAM (SParse Additive Modeling).
#'
#' @param x        Features to be expanded (numeric matrix).
#' @param df       Degrees of freedom (numeric; default = 3).
#' @param degree   Degree of the piecewise polynomial (integer; default = 3 (cubic splines)).
#' @param type     Type of spline, either B-spline (`"bs"`) or natural cubic spline (`"ns"`; default).
#' 
#' @return
#' An object of class `expandedMatrix` consisting of:
#' * `X`: A matrix of dimension `nrow(x)` by `df*ncol(x)`
#' * `group`: A vector of length `df*ncol(x)` that describes the grouping structure
#' * Additional metadata on the splines, such as knot locations, required in order to evaluate spline at new feature values (e.g., for prediction)
#' 
#' @references
#'   * Ravikumar P, Lafferty J, Liu H and Wasserman L (2009). Sparse additive models. *Journal of the Royal Statistical Society Series B*, **71**: 1009-1030.
#'   
#' @seealso [plot_spline()] to visualize the resulting nonlinear fits
#'
#' @examples
#' \dontshow{set.seed(1)}
#' Data <- gen_nonlinear_data(n=1000)
#' X <- expand_spline(Data$X)
#' fit <- grpreg(X, Data$y)
#' plot_spline(fit, "V02", lambda = 0.03)
#' @export

expand_spline <- function(x, df = 3, degree = 3, type = c("ns", "bs")) {
  type <- match.arg(type)
  if(type == "ns" && degree != 3) {
    warning("Degree has been set to 3 for natural splines")
  }
  n <- nrow(x)
  p <- ncol(x)
  finalx <- matrix(NA, n, (p*df))
  knots <- rep(list(rep(NA, (df-degree))), p)
  boundary <- rep(list(rep(NA, 2)), p)
  
  if(type == "bs"){
    for(i in 0:(p-1)){
      bs <- splines::bs(x[,i+1], df = df, degree = degree)
      finalx[,(df*i+1):(df*i+df)] <- bs
      boundary[[i+1]] <- attr(bs, "Boundary.knots")
      knots[[i+1]] <- attr(bs, "knots")
    }
  }
  else if (type == "ns") {
    for (i in 0:(p-1)) {
      ns <- splines::ns(x[,i+1], df = df)
      finalx[,(df*i+1):(df*i+df)] <- ns
      boundary[[i+1]] <- attr(ns, "Boundary.knots")
      knots[[i+1]] <- attr(ns, "knots")
    }
  }
  else {
    stop(paste(type, "is not a valid type"), call. = FALSE)
  }
  
  if (length(colnames(x)) == p) {
    groups <- rep(colnames(x), each = df)
    colnames(finalx) <- paste(groups, 1:df, sep='_')
  }
  else {
    groups <- rep(paste0("V", 1:p), each = df)
    colnames(finalx) <- paste(groups, 1:df, sep="_")
  }
  
  return(structure(list(X = finalx, 
                        group = groups, 
                        knots = knots,
                        boundary = boundary,
                        degree = degree,
                        originalx = x,
                        type = type), class='expandedMatrix'))
}
