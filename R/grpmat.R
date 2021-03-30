#' Spline Basis for Multiple Variables
#'
#' \code{grpmat} Performs a basis expansion for many variables at once for 
#' later use with the grpreg function. Fuction returns one large matrix 
#' and a vector that describes its grouping.
#'
#' \code{grpmat} is based on the function [splines::bs()]. It takes each
#' column of \code{x} and generates a basis matrix. These basis 
#' matrices represent the family of piecewise polynomials with the specified 
#' degree evaluated at the column values of \code{x}. These matrices are then 
#' column-bound to form a single grouped design matrix. A vector that describes 
#' the grouping present in the resulting matrix is also generated. The resulting 
#' object can be passed to [grpreg()].
#'
#' @param x the design matrix. Columns must represent numeric variables.
#' @param df degrees of freedom; default is 4.
#' @param degree degree of the piecewise polynomial; default is 3 for cubic 
#' splines.
#' @param type specifies type of splines, default is \code{"bs"}
#' @return A \code{grouped_hat} object composed of a matrix of dimension 
#' \code{c(nrow(x), df*ncol(x))} and a vector of length 
#' \code{df*ncol(x)}
#'
#' @examples
#' 
#' X <- grpmat(attitude[-1], df = 3)
#' fit <- grpreg(X, attitude$rating, penalty="grLasso")
#' plot(fit)
#'
#'

grpmat <- function(x, df = 4, degree = 3, type = "bs"){
  if(type == "ns" && degree != 3){
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
  else if(type == "ns"){
    for(i in 0:(p-1)){
      ns <- splines::ns(x[,i+1], df = df)
      finalx[,(df*i+1):(df*i+df)] <- ns
      boundary[[i+1]] <- attr(ns, "Boundary.knots")
      knots[[i+1]] <- attr(ns, "knots")
    }
  }
  else{
    stop(paste(type, "is not a valid type"))
  }
  
  if(length(colnames(x)) == p){
    groups <- rep(colnames(x), each = df)
    colnames(finalx) <- paste0(groups, 1:df)
  }
  else{
    groups <- rep(paste0("V", 1:p), each = df)
    colnames(finalx) <- paste(groups, 1:df, sep = "_")
  }
  
  return(structure(list(x = finalx, 
                        groups = groups, 
                        knots = knots,
                        boundary = boundary,
                        degree = degree,
                        originalx = x,
                        type = type), class='grouped_mat'))
}
  
