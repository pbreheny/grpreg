unstandardize <- function(b, XG) {
  beta <- matrix(0, nrow=1+length(XG$scale), ncol=ncol(b))
  beta[1 + XG$nz,] <- b[-1,] / XG$scale[XG$nz]
  beta[1,] <- b[1,] - crossprod(XG$center, beta[-1,,drop=FALSE])
  beta
}
