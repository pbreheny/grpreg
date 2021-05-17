predict.mat <- function(object, X, indic){
  p <- sum(indic)
  grpmat <- object$grpmat
  df <- length(grpmat$knots[[1]])+grpmat$degree
  bsX <- matrix(NA, nrow(X), (p*df)) #make this work for vectors
  i <- 1
  for(i in 0:(p-1)){
    mat <- object$X[,(df*i+1):(df*i+df)]
    attr(mat, "degree") <- grpmat$degree
    attr(mat, "knots") <- grpmat$knots[[i+1]]
    attr(mat, "Boundary.knots") <- grpmat$boundary[[i+1]]
    attr(mat, "intercept") <- FALSE
    attr(mat, "class") <- c(grpmat$type, "basis", "matrix")
    bsX[,(df*i+1):(df*i+df)] <- predict(mat, X[,i+1])
  }
  X <- bsX
}