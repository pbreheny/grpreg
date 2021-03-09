plot.sp.grpreg <- function(x, fit, variable, lambda, xlab = variable, ...){
  #reproduce bs object - will need to check if bs or ns
  #check if x is grped mat
  df <- length(x$knots[[1]]+x$degree)
  j <- which(x$groups == variable)
  i <- j[df]/df
  mat <- x$x[,j]
  attr(mat, "degree") <- x$degree
  attr(mat, "knots") <- x$knots[[i]]
  attr(mat, "Boundary.knots") <- x$boundary[[i]]
  attr(mat, "intercept") <- FALSE
  attr(mat, "class") <- c("bs", "basis", "matrix")
  
  #create sequence and basis
  min <- attr(mat, "Boundary.knots")[1]
  max <- attr(mat, "Boundary.knots")[2]
  newx <- seq(min, max, length.out = 200) 
  newxbs <- cbind(1, predict(mat, newx))
  
  #plot
  cols <- hcl(h=seq(15, 240, len=length(lambda)), l=60, c=150, alpha=1)
  notsel <- NULL
  y <- newxbs%*%fit$beta[c(1,j+1), lambda[1]]
  plot(newx, y, xlab = xlab, type = "l", col = cols[1])
  if(length(lambda > 1)){
    for(i in 2:(length(lambda))){
      y <- newxbs%*%fit$beta[c(1,j+1), lambda[i]]
      lines(newx, y, col = cols[i])
      if(beta[j+1, lambda[i]] == rep(0, length(j))){
        notsel <- c(notsel, lambda[i])
      }
    }
  }
  #points(x$originalx[,i], fit$y)
}
