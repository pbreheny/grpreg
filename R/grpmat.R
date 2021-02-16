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
#' @return A \code{grouped_hat} object composed of a matrix of dimension 
#' \code{c(nrow(x), df*ncol(x))} and a vector of length 
#' \code{df*ncol(x)}
#'
#' @examples
#' grpmat(examplestuff)
#'
#'

grpmat <- function(x, df = 4, degree = 3){
  n <- nrow(x)
  p <- ncol(x)
  finalx <- matrix(1:40, n, (p*df))
  
  for(i in 0:(p-1)){
    finalx[,(4*i+1):(4*i+df)] <- splines::bs(x[,i+1], df = df, degree = degree)
  }
  
  if(length(colnames(x)) == p){
    groups <- rep(colnames(x), each = df)
    colnames(finalx) <- paste0(groups, 1:df)
  }
  else{
    groups <- rep(paste0("V", 1:p), each = df)
    colnames(finalx) <- paste(groups, 1:df, sep = "_")
  }
  
  return(structure(list(x = finalx, groups = groups), class='grouped_hat'))
}
