# Bootstrap method for calculating CVSE for Cox models
se_grpsurv <-function (y, eta, B = 100) {
  cve <- matrix(NA, B, ncol(eta))
  for (b in 1:B) {
    ind <- sample(1:nrow(eta), replace=TRUE)
    cve[b,] <- deviance_grpsurv(y[ind,], eta[ind,])
  }
  apply(cve, 2, sd)
}
