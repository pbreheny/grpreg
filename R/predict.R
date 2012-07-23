predict.grpreg <- function(object, X, lambda, which=1:length(object$lambda), type=c("link","response","class","coefficients"),...)
{
  type <- match.arg(type)
  beta <- coef.grpreg(object, lambda=lambda, which=which)
  if (type=="coefficients") return(drop(object$beta[,which]))
  if (missing(X)) stop("Must supply X")
  eta <- sweep(X %*% beta[-1,,drop=FALSE], 2, beta[1,], "+")
  if (object$family=="gaussian" & type=="class") stop("type='class' is not applicable for family='gaussian'")
  if (object$family=="gaussian" | type=="link") return(drop(eta))
  pihat <- exp(eta)/(1+exp(eta))
  if (type=="response") return(drop(pihat))
  if (type=="class") return(drop(eta > 0))
}
coef.grpreg <- function(object, lambda, which=1:length(object$lambda), ...)
{
  if (!missing(lambda)) {
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    beta <- (1-w)*object$beta[,l,drop=FALSE] + w*object$beta[,r,drop=FALSE]
    if (length(lambda) > 1) colnames(beta) <- round(lambda,4)
  }
  else beta <- object$beta[, which, drop=FALSE]
  return(beta)
}
