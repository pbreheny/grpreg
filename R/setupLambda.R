setupLambda <- function(X, y, group, family, penalty, alpha, lambda.min, nlambda, group.multiplier)
{
  ## Fit to unpenalized covariates
  n <- length(y)
  ind <- which(group!=0)
  if (length(ind)!=length(group)) {
    fit <- glm(y~X[, group==0], family=family)
  } else fit <- glm(y~1, family=family)
  
  ## Determine lambda.max
  if (family=="gaussian") {
    z <- crossprod(X[,ind], fit$residuals) / n
  } else {
    z <- crossprod(X[,ind], fit$weights * residuals(fit, "working")) / n
    ## v <- apply(X[ ,ind, drop=FALSE] * fit$weights * X[ ,ind, drop=FALSE], 2, sum) / n
    ## z <- u/v
  }
  if (strtrim(penalty,2)=="gr") maxGradient <- sqrt(tapply(z^2, group[ind], sum))
  if (strtrim(penalty,2)=="ge") maxGradient <- tapply(abs(z),group[ind],max)
  if (penalty=="gMCP") maxGradient <- sqrt(tapply(abs(z),group[ind],max))
  lambda.max <- max(maxGradient/group.multiplier) / alpha
  
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 0)
  } else lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  lambda
}

setupLambda.gBridge <- function(X, y, group, family, alpha, lambda.min, lambda.max, nlambda, gamma, group.multiplier)
{
  ## Fit to unpenalized covariates
  n <- length(y)
  ind <- which(group!=0)
  if (length(ind)!=length(group)) {
    fit <- glm(y~X[, group==0], family=family)
  } else {
    fit <- glm(y~1, family=family)
  }
  
  ## Guess lambda.max
  if (missing(lambda.max)) {
    if (family=="gaussian") {
      z <- crossprod(X[,ind], fit$residuals) / n
      a <- .35
    } else {
      z <- crossprod(X[,ind], fit$weights * residuals(fit, "working")) / n
      a <- .2
    }
    maxGradient <- tapply(abs(z), group[ind],max)*a^(1-gamma)/gamma
    lambda.max <- max(maxGradient/group.multiplier) / alpha
  }
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)),0)                  
  } else {
    lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  }
  return(rev(lambda))
}
