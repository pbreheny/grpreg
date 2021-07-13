newS <- function(y) {
  if(ncol(y) == 2){
    ind <- order(y[,1])
    list(start_time = as.double(rep(0,nrow(y))),
         stop_time = as.double(y[ind,1]),
         fail = as.double(y[ind,2]),
         ind  = ind)
  }else if(ncol(y)==3){
    ind <- order(y[,2])
    list(start_time = as.double(y[ind,1]),
         stop_time = as.double(y[ind,2]),
         fail = as.double(y[ind,3]),
         ind  = ind)
  }else{
    stop('y must be specified as Surv(time,fail) or Surv(start,stop,fail)')
  }
}
