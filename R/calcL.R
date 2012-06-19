calcL <- function(X,y,beta,family)
  {
    L <- ncol(beta)
    val <- numeric(L)
    if (family=="gaussian")
      {
        for (l in 1:L)
          {
            val[l] <- apply((y - X %*% beta[,l])^2,2,sum)/2
          }
      }
    if (family=="binomial")
      {
        for (l in 1:L)
          {
            eta <- X %*% beta[,l]
            pi. <- exp(eta)/(1+exp(eta))
            val[l] <- -sum(y*log(pi.)+(1-y)*log(1-pi.))
          }
      }
    return(val)
  }
