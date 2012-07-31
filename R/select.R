select.grpreg <- function(obj, criterion=c("BIC","AIC","GCV"), df.method=c("default","active"), smooth=FALSE, ...)
{
  criterion <- match.arg(criterion)
  df.method <- match.arg(df.method)
  ll <- logLik(obj, df.method=df.method, ...)
  df <- as.numeric(attr(ll,"df"))
  
  if (criterion=="AIC") IC <- AIC(ll)
  if (criterion=="BIC") IC <- BIC(ll)
  if (criterion=="GCV") IC <- (1/obj$n) * (-2) * as.numeric(ll) / (1-df/obj$n)^2
  n.l <- length(obj$lambda)
  if (smooth & (n.l < 4)) {
    smooth <- FALSE
    warning("Need at least 4 points to use smooth=TRUE")
  }
  if (smooth) {
    fit.ss <- smooth.spline(IC[is.finite(IC)])
    ##plot(obj$lambda,IC,pch=19,xlim=rev(range(obj$lambda)))
    ##lines(obj$lambda,fit.ss$y)
    d <- diff(fit.ss$y)
    if (all(d<0)) i <- n.l
    else i <- min(which(d>0))-1
    if (i==0) i <- 1
  } else i <- which.min(IC)
  
  if (min(obj$lambda) == obj$lambda[i]) warning(paste("minimum lambda selected for",obj$penalty))
  else if ((max(obj$lambda) == obj$lambda[i]) & obj$penalty=="gBridge") warning("maximum lambda selected")
  return(list(beta=obj$beta[,i],
              lambda=obj$lambda[i],
              df=df[i],
              IC=IC))
}
select <- function(obj,...) UseMethod("select")
