# library(splines)
# 
notTheFinalName <- function(x, df = 4, degree = 3){
  n <- nrow(x)
  p <- ncol(x)
  finalx <- matrix(1:40, n, (p*df))
  
  for(i in 0:(p-1)){
    finalx[,(4*i+1):(4*i+df)] <- bs(x[,i+1], df = df, degree = degree)
  }
  
  if(length(colnames(x)) == p){
    colnames(finalx) <- paste0(colnames(x), rep(1:df, times = p))
  }
  
  groups <- rep(paste0("x", 1:p), each = df)
  
  return(structure(list(x = finalx, groups = groups), class='grouped_hat'))
}
# 
# n <- 100
# p <- 1000
# x <- matrix(runif(100), n, p, byrow = TRUE)
# hat <- notTheFinalName(x)
# hat
# 
# x <- matrix(1:20, 10,2)
# 
# system.time(x_list <- split(x, col(x, as.factor = TRUE))
# x_list <- lapply(x_list, FUN = bs, df = 4, degree = 3)
# matrix(unlist(x_list), ncol = 2*4))
# 
# system.time(for(i in 0:(p-1)){
#   finalx[,(4*i+1):(4*i+df)] <- bs(x[,i+1], df = df, knots = knots, degree = degree)
# })
