predict.grpreg <- function(object,X,which=1:length(object$lambda),type=c("link","response","class","coefficients"),...)
  {
    type <- match.arg(type)
    if (type=="coefficients") return(object$beta[,which])
    eta <- t(object$beta[1,which] + t(X%*%object$beta[-1,which]))
    if (object$family=="gaussian" | type=="link") return(eta)
    pihat <- exp(eta)/(1+exp(eta))
    if (type=="response") return(pihat)
    if (type=="class") return(eta>0)
  }
coef.grpreg <- function(object,which=1:length(object$lambda),...)
  {
    return(object$beta[,which])
  }
