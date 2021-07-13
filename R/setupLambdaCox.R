setupLambdaCox <- function(X, y, group, penalty, alpha, lambda.min, nlambda, group.multiplier) {
  n <- nrow(X)
  p <- ncol(X)

  ## Fit to unpenalized covariates
  K <- table(group)
  K1 <- as.integer(if (min(group)==0) cumsum(K) else c(0, cumsum(K)))
  if (K1[1]!=0) {
    SURV <- get("Surv", asNamespace("survival"))
    COXPH <- get("coxph", asNamespace("survival"))
    nullFit <- COXPH(SURV(y$start_time, y$stop_time, y$fail) ~ X[, group==0, drop=FALSE])
    eta <- nullFit$linear.predictors
    rsk <- rev(cumsum(rev(exp(eta))))
    current_sum = 0
    start_idx = n
    stop_idx = n
    start_o = order(y[,1])
    while(start_idx>0 && stop_idx >0){
      if(y[start_o[start_idx],1]<y[ind[stop_idx],2]){
        rsk[stop_idx] = rsk[stop_idx] - current_sum
        stop_idx = stop_idx - 1
      }else{
        current_sum = current_sum + exp(eta[start_o[start_idx]])
        start_idx = start_idx - 1
      }
    }
    s <- y$fail - exp(eta)*cumsum(y$fail/rsk)
  } else {
    w <- 1/(n-(1:n)+1)
    s <- y$fail - cumsum(y$fail*w)
  }

  ## Determine lambda.max
  if (strtrim(penalty, 2) == "gr") {
    zmax <- .Call("maxgrad", X, s, K1, as.double(group.multiplier)) / n
  } else {
    zmax <- .Call("maxprod", X, s, K1, as.double(group.multiplier)) / n
  }
  lambda.max <- zmax/alpha

  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
  else lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
  lambda
}
