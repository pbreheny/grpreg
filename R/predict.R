predict.grpreg <- function(object, X, lambda, which=1:length(object$lambda), type=c("link","response","class","coefficients","vars","groups","norm"), ...) {
  type <- match.arg(type)
  beta <- coef.grpreg(object, lambda=lambda, which=which, drop=FALSE)
  if (type=="coefficients") return(beta)
  if (length(dim(object$beta)) == 2) {
    if (type=="vars") return(drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which)))
    if (type=="groups") return(drop(apply(beta[-1, , drop=FALSE]!=0, 2, function(x) unique(object$group[x]))))
    if (type=="norm") return(drop(apply(beta[-1, , drop=FALSE], 2, function(x) tapply(x, object$group, l2))))
    if (missing(X)) stop("Must supply X")
    eta <- sweep(X %*% beta[-1,,drop=FALSE], 2, beta[1,], "+")
    if (object$family=="gaussian" & type=="class") stop("type='class' is not applicable for family='gaussian'")
    if (object$family=="gaussian" | type=="link") return(drop(eta))
    pihat <- exp(eta)/(1+exp(eta))
    if (type=="response") return(drop(pihat))
    if (type=="class") return(drop(eta > 0))
  } else {
    if (type=="vars") return(apply(beta[,-1, , drop=FALSE], 1, function(x){apply(x!=0, 2, FUN=which)}))
    if (type=="groups") return(drop(apply(beta[,-1, , drop=FALSE], 3, function(x){which(apply(x!=0, 2, any))})))
    if (type=="norm") return(drop(apply(beta[, -1, , drop=FALSE], 3, function(x) apply(x, 2, l2))))
    if (missing(X)) stop("Must supply X")
    eta <- apply(beta[,-1,,drop=FALSE],1,function(b){X%*%b})
    eta <- array(eta, dim=c(nrow(X), dim(beta)[1], dim(beta)[3]), dimnames=list(NULL, dimnames(beta)[[1]], dimnames(beta)[[3]]))
    eta <- sweep(eta, 2:3, beta[,1,], "+")
    if (object$family=="gaussian" & type=="class") stop("type='class' is not applicable for family='gaussian'")
    if (object$family=="gaussian" | type=="link") return(drop(eta))
    pihat <- exp(eta)/(1+exp(eta))
    if (type=="response") return(drop(pihat))
    if (type=="class") return(drop(eta > 0))    
  }
}
coef.grpreg <- function(object, lambda, which=1:length(object$lambda), drop=TRUE, ...) {
  if (!missing(lambda)) {
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    if (length(dim(object$beta)) == 3) {
      beta <- (1-w)*object$beta[,,l,drop=FALSE] + w*object$beta[,,r,drop=FALSE]
      dimnames(beta)[[3]] <- round(lambda,4)      
    } else {
      beta <- (1-w)*object$beta[,l,drop=FALSE] + w*object$beta[,r,drop=FALSE]
      colnames(beta) <- round(lambda,4)
    }
  } else {
    if (length(dim(object$beta)) == 3) {
      beta <- object$beta[,,which, drop=FALSE]
    } else {
      beta <- object$beta[, which, drop=FALSE]
    }
  }
  if (drop) return(drop(beta)) else return(beta)
}
l2 <- function(x) sqrt(sum(x^2))