predict_spline <- function(object, X) {
  p <- length(unique(object$group))
  meta <- object$meta
  if (meta$type == 'bs') {
    df <- length(meta$knots[[1]]) + meta$degree
  } else if (meta$type == 'ns') {
    df <- length(meta$knots[[1]]) + 1
  }
  bsX <- matrix(NA, nrow(X), (p*df)) # make this work for vectors
  for(i in 1:p) {
    mat <- object$meta$X[,(df*(i-1)+1):(df*(i-1)+df)]
    attr(mat, "degree") <- meta$degree
    attr(mat, "knots") <- meta$knots[[i]]
    attr(mat, "Boundary.knots") <- meta$boundary[[i]]
    attr(mat, "intercept") <- FALSE
    attr(mat, "class") <- c(meta$type, "basis", "matrix")
    bsX[,(df*(i-1)+1):(df*(i-1)+df)] <- predict(mat, X[,i])
  }
  bsX
}
