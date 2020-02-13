newS <- function(y) {
  ind <- order(y[,1])
  list(time = as.double(y[ind,1]),
       fail = as.double(y[ind,2]),
       ind  = ind)
}
