newS <- function(y) {
  ind <- order(y[,1])
  list(time = as.numeric(y[ind,1]),
       fail = as.numeric(y[ind,2]),
       ind  = ind)
}
