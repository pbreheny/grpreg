setupLambda <- function(X, y, group, family, penalty, alpha, lambda.min, nlambda, gamma, group.multiplier)
{
  ## Fit to unpenalized covariates
  n <- length(y)
  ind <- which(group!=0)
  if (length(ind)!=length(group)) {
    fit <- glm(y~X[, group==0], family=family)
  } else fit <- glm(y~1, family=family)
  
  ## Determine lambda.max
  if (family=="gaussian") {
    z <- crossprod(X[,ind],fit$residuals)
  } else z <- crossprod(X[,ind], fit$weights*residuals(fit, "working"))
  if (strtrim(penalty,2)=="gr") maxGradient <- sqrt(tapply(z^2,group[ind],sum))/n
  if (strtrim(penalty,2)=="ge") maxGradient <- tapply(abs(z),group[ind],max)/n
  if (penalty=="gMCP") maxGradient <- sqrt(tapply(abs(z),group[ind],max)/n)
  lambda.max <- max(maxGradient/group.multiplier) / alpha
  
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 0)
  } else lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  lambda
}

setupLambda.gBridge <- function(X, y, group, family, penalty, lambda.max, lambda.min, nlambda, gamma, group.multiplier)
{
  n <- length(y)
  K <- as.numeric(table(group))[-1]
  if (missing(lambda.max))
  {
    fit <- glm(y~0+X[,group==0],family=family)
    if (family=="gaussian")
    {
      r <- fit$residuals
      if (penalty=="gLasso") maxGradient <- sqrt(tapply(crossprod(X[,group!=0],r)^2,group[group!=0],sum))/n
      if (penalty=="gBridge") maxGradient <- tapply(abs(crossprod(X[,group!=0],r)),group[group!=0],max)*.35^(1-gamma)/(n*gamma*K^gamma)
      if (penalty=="gMCP") maxGradient <- sqrt(tapply(abs(crossprod(X[,group!=0],r)),group[group!=0],max)/n)
      if (penalty=="geLasso" | penalty=="geMCP") maxGradient <- tapply(abs(crossprod(X[,group!=0],r)),group[group!=0],max)/n
    }
    if (family=="binomial")
    {
      eta <- X[,group==0]%*%matrix(fit$coef,ncol=1)
      pi. <- exp(eta)/(1+exp(eta))
      w <- as.numeric(sqrt(pi.*(1-pi.)))
      r = (y - pi.)/w;
      if (penalty=="gLasso") maxGradient <- sqrt(tapply(crossprod(w*X[,group!=0],r)^2,group[group!=0],sum))/n
      if (penalty=="gBridge") maxGradient <- tapply(abs(crossprod(w*X[,group!=0],r)),group[group!=0],max)*.2^(1-gamma)/(n*gamma*K^gamma)
      if (penalty=="gMCP") maxGradient <- sqrt(tapply(abs(crossprod(w*X[,group!=0],r)),group[group!=0],max)/n)
      if (penalty=="geLasso" | penalty=="geMCP") maxGradient <- tapply(abs(crossprod(w*X[,group!=0],r)),group[group!=0],max)/n
    }
    lambda.max <- max(maxGradient/group.multiplier)+1e-6
  }
  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  if (penalty=="gBridge") lambda <- rev(lambda)
  return(lambda)
}
