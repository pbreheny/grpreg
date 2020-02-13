predict.grpreg <- function(object, X, type=c("link", "response", "class", "coefficients", "vars", "groups", "nvars", "ngroups", "norm"), lambda, which=1:length(object$lambda), ...) {
  if (!missing(X) && is.character(X)) {
    type <- X
    X <- NULL
  }
  type <- match.arg(type)
  beta <- coef.grpreg(object, lambda=lambda, which=which, drop=FALSE)
  if (type=="coefficients") return(beta)
  d <- length(dim(beta))
  if (inherits(object, 'grpsurv')) {
    beta <- beta
  } else {
    if (d==3) {
      alpha <- beta[,1,]
      beta <- beta[, -1, , drop=FALSE]
    } else {
      alpha <- beta[1,]
      beta <- beta[-1, , drop=FALSE]
    }
  }
  if (d == 2) {
    if (type=="vars") return(drop(apply(beta!=0, 2, FUN=which)))
    if (type=="groups") return(drop(apply(beta!=0, 2, function(x) unique(object$group[x]))))
    if (type=="nvars") {
      v <- drop(apply(beta!=0, 2, FUN=which))
      if (is.list(v)) {
        res <- sapply(v, length)
      } else {
        res <- length(v)
      }
      return(res)
    }
    if (type=="ngroups") {
      g <- drop(apply(beta!=0, 2, function(x) unique(object$group[x])))
      if (is.list(g)) {
        res <- sapply(g, length)
      } else {
        res <- length(g)
      }
      return(res)
    }
    if (type=="norm") return(drop(apply(beta, 2, function(x) tapply(x, object$group, function(x){sqrt(sum(x^2))}))))
    if (missing(X) | is.null(X)) stop("Must supply X", call.=FALSE)
    eta <- sweep(X %*% beta, 2, alpha, "+")
    if (object$family=="gaussian" & type=="class") stop("type='class' is not applicable for family='gaussian'", call.=FALSE)
    if (object$family=="gaussian" | type=="link") return(drop(eta))
    resp <- switch(object$family,
                   binomial = exp(eta)/(1+exp(eta)),
                   poisson = exp(eta))
    if (type=="response") return(drop(resp))
    if (type=="class") {
      if (object$family=="binomial") {
        return(drop(1*(eta>0)))
      } else {
        stop("type='class' can only be used with family='binomial'", call.=FALSE)
      }
    }
  } else {
    if (type=="vars") stop("Predicting type 'vars' not implemented with multivariate outcomes", call.=FALSE)
    if (type=="groups") return(drop(apply(beta, 3, function(x){which(apply(x!=0, 2, any))})))
    if (type=="norm") return(drop(apply(beta, 3, function(x) apply(x, 2, function(x){sqrt(sum(x^2))}))))
    if (type=="nvars") {
      return(drop(apply(beta!=0, 3, FUN=sum)))
    }
    if (type=="ngroups") {
      return(apply(apply(beta!=0, c(2,3), any), 2, sum))
    }
    if (missing(X)) stop("Must supply X", call.=FALSE)
    eta <- apply(beta, 1, function(b){X%*%b})
    eta <- array(eta, dim=c(nrow(X), dim(beta)[1], dim(beta)[3]), dimnames=list(NULL, dimnames(beta)[[1]], dimnames(beta)[[3]]))
    eta <- sweep(eta, 2:3, alpha, "+")
    if (object$family=="gaussian" & type=="class") stop("type='class' is not applicable for family='gaussian'", call.=FALSE)
    if (object$family=="gaussian" | type=="link") return(drop(eta))
    resp <- switch(object$family,
                   binomial = exp(eta)/(1+exp(eta)),
                   poisson = exp(eta))
    if (type=="response") return(drop(resp))
    if (type=="class") {
      if (object$family=="binomial") {
        return(drop(1*(eta>0)))
      } else {
        stop("type='class' can only be used with family='binomial'", call.=FALSE)
      }
    }
  }
}
coef.grpreg <- function(object, lambda, which=1:length(object$lambda), drop=TRUE, ...) {
  if (!missing(lambda)) {
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    if (length(dim(object$beta)) == 3) {
      beta <- (1-w)*object$beta[, , l, drop=FALSE] + w*object$beta[, , r, drop=FALSE]
      dimnames(beta)[[3]] <- round(lambda, 4)
    } else {
      beta <- (1-w)*object$beta[, l, drop=FALSE] + w*object$beta[, r, drop=FALSE]
      colnames(beta) <- round(lambda, 4)
    }
  } else {
    if (length(dim(object$beta)) == 3) {
      beta <- object$beta[, , which, drop=FALSE]
    } else {
      beta <- object$beta[, which, drop=FALSE]
    }
  }
  if (drop) return(drop(beta)) else return(beta)
}
