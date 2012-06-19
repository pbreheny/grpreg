unstandardize <- function(beta,meanx,normx)
  {
    val <- beta
    val[1,] <- beta[1,]-apply(meanx*beta[-1,,drop=F]/normx,2,sum)
    val[-1,] <- beta[-1,]/normx
    return(val)
  }
