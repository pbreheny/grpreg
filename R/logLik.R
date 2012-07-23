logLik.grpreg <- function(object, df.method=c("default","active"), REML=FALSE, ...)
{
  df.method <- match.arg(df.method)
  n <- as.numeric(object$n)
  df <- if (df.method=="active") apply(coef(object)!=0, 2, sum) else object$df
  if (object$family=="gaussian")
  {
    rdf <- if (REML) n-df else n
    RSS <- object$loss
    l <- -n/2 * (log(2*pi) + log(RSS) - log(rdf)) - rdf/2
    df <- df + 1
  }
  if (object$family=="binomial") l <- -1*object$loss
  
  val <- l
  attr(val,"df") <- df
  attr(val,"nobs") <- n
  class(val) <- "logLik"
  return(val)    
}
