logLik.grpreg <- function(object, df.method=c("default","active"), ...)
  {
    df.method <- match.arg(df.method)
    if (df.method=="default") df <- object$df
    else if (df.method=="active") df <- apply(object$beta!=0,2,sum)

    if (object$family=="gaussian")
      {
        n <- object$n
        RSS <- 2*object$loss
        l <- -n/2 * (log(2*pi) + 1 - log(n) + log(RSS))
        df <- df + 1
      }
    if (object$family=="binomial") l <- -1*object$loss

    val <- l
    attr(val,"df") <- as.numeric(df)
    attr(val,"nobs") <- as.numeric(object$n)
    class(val) <- "logLik"
    return(val)
  }
