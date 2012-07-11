standardize <- function(X)
{
  n <- nrow(X)
  center <- colMeans(X)
  X.c <- sweep(X, 2, center)
  scale <- sqrt(apply(X.c,2,crossprod)/n)
  val <- sweep(X.c,2,scale,"/")
  attr(val,"center") <- center
  attr(val,"scale") <- center
  val
}
