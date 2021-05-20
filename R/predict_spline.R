predict_spline <- function(object, X){
  p <- length(unique(object$group))
  meta <- object$meta
  df <- length(meta$knots[[1]]) + meta$degree
  bsX <- matrix(NA, nrow(X), (p*df)) #make this work for vectors
  i <- 1
  for(i in 0:(p-1)){
    mat <- object$X[,(df*i+1):(df*i+df)]
    attr(mat, "degree") <- meta$degree
    attr(mat, "knots") <- meta$knots[[i+1]]
    attr(mat, "Boundary.knots") <- meta$boundary[[i+1]]
    attr(mat, "intercept") <- FALSE
    attr(mat, "class") <- c(meta$type, "basis", "matrix")
    bsX[,(df*i+1):(df*i+df)] <- predict(mat, X[,i+1])
  }
  X <- bsX
}