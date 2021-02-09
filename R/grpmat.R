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
