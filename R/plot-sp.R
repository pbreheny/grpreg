plot.sp.grpreg <- function(fit, variable, lambda, which, scatter = FALSE, type = "contrast", ...){
  #check if x is grped mat
  if(scatter == TRUE && length(which) > 1){
    warning("Scatter plot represents largest value of lambda imputed")
  }
  if(missing(lambda)){
    lambda <- fit$lambda[which]
  }
  
  grpmat <- fit$grpmat
  df <- length(grpmat$knots[[1]])+grpmat$degree
  j <- which(fit$group == variable)
  i <- j[df]/df
  l <- length(which)
  mat <- fit$X[,j]
  attr(mat, "degree") <- grpmat$degree
  attr(mat, "knots") <- grpmat$knots[[i]]
  attr(mat, "Boundary.knots") <- grpmat$boundary[[i]]
  attr(mat, "intercept") <- FALSE
  attr(mat, "class") <- c(grpmat$type, "basis", "matrix")
  
  #create sequence and basis and calculate y's and residuals
  min <- grpmat$boundary[[i]][1]
  max <- grpmat$boundary[[i]][2]
  newx <- seq(min, max, length.out = 200) 
  if(type == "conditional"){
    xmeans <- matrix(colMeans(grpmat$originalx), 200, ncol(grpmat$originalx), byrow = TRUE)
    xmeans[,i] <- newx
    y <- predict(fit, xmeans, lambda = lambda)
    #parresid <- fit$y - cbind(1, fit$X[,-j])%*%fit$beta[-(j+1), max(which)]
  } else if(type == "contrast"){
    newxbs <- predict(mat, newx)
    newxmean <- predict(mat, mean(grpmat$originalx[,i]))
    y <- newxbs%*%fit$beta[j+1, which] - matrix(newxmean%*%fit$beta[j+1, which],200,length(which), byrow = TRUE)
    #parresid <- fit$y - cbind(1, fit$X[,-j])%*%fit$beta[-(j+1), max(which)]
  } else {
    stop(paste(type, "is not a valid type"))
  }
  
  #check for selection
  #notselected <- NULL
  #for(i in 1:length(which)){
  #  if (all(fit$beta[(j+1), which[i]]==0)){
  #    notselected <- c(notselected, fit$lambda[which[i]])
  #  }
  #}
  #if(length(notselected) > 0){
  #  warning(paste(variable, "was not selected at lambda =", notselected))
  #}
  
  #plot
  cols <- hcl(h=seq(15, 240, len=length(which)), l=60, c=150, alpha=1)
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
    #points(fit$originalx[,i], parresid)
  }
}
