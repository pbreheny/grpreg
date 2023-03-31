#' Model predictions based on a fitted \code{grpreg} object
#' 
#' Similar to other predict methods, this function returns predictions from a
#' fitted \code{"grpreg"} object.
#' 
#' \code{coef} and \code{predict} methods are provided for \code{"cv.grpreg"}
#' options as a convenience.  They simply call \code{coef.grpreg} and
#' \code{predict.grpreg} with \code{lambda} set to the value that minimizes the
#' cross-validation error.
#' 
#' @aliases predict.grpreg coef.grpreg predict.cv.grpreg coef.cv.grpreg
#' 
#' @param object Fitted \code{"grpreg"} or \code{"cv.grpreg"} model object.
#' @param X Matrix of values at which predictions are to be made.  Not used for
#' \code{type="coefficients"}
#' @param lambda Values of the regularization parameter \code{lambda} at which
#' predictions are requested.  For values of \code{lambda} not in the sequence
#' of fitted models, linear interpolation is used.
#' @param which Indices of the penalty parameter \code{lambda} at which
#' predictions are required.  By default, all indices are returned.  If
#' \code{lambda} is specified, this will override \code{which}.
#' @param type Type of prediction: \code{"link"} returns the linear predictors;
#' \code{"response"} gives the fitted values; \code{"class"} returns the
#' binomial outcome with the highest probability; \code{"coefficients"} returns
#' the coefficients; \code{"vars"} returns the indices for the nonzero
#' coefficients; \code{"groups"} returns the indices for the groups with at
#' least one nonzero coefficient; \code{"nvars"} returns the number of nonzero
#' coefficients; \code{"ngroups"} returns the number of groups with at least
#' one nonzero coefficient; \code{"norm"} returns the L2 norm of the
#' coefficients in each group.
#' @param drop By default, if a single value of \code{lambda} is supplied, a
#' vector of coefficients is returned.  Set \code{drop=FALSE} if you wish to
#' have \code{coef} always return a matrix (see \code{\link{drop}}).
#' @param \dots Not used.
#' @return The object returned depends on type.
#' @author Patrick Breheny
#' @seealso \code{grpreg}
#' @examples
#' # Fit penalized logistic regression model to birthweight data
#' data(Birthwt)
#' X <- Birthwt$X
#' y <- Birthwt$low
#' group <- Birthwt$group
#' fit <- grpreg(X, y, group, penalty="grLasso", family="binomial")
#' 
#' # Coef and predict methods
#' coef(fit, lambda=.001)
#' predict(fit, X, type="link", lambda=.07)[1:10]
#' predict(fit, X, type="response", lambda=.07)[1:10]
#' predict(fit, X, type="class", lambda=.01)[1:15]
#' predict(fit, type="vars", lambda=.07)
#' predict(fit, type="groups", lambda=.07)
#' predict(fit, type="norm", lambda=.07)
#' 
#' # Coef and predict methods for cross-validation
#' cvfit <- cv.grpreg(X, y, group, family="binomial", penalty="grMCP")
#' coef(cvfit)
#' predict(cvfit, X)[1:10]
#' predict(cvfit, X, type="response")[1:10]
#' predict(cvfit, type="groups")
#' @export

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
    if (type=="groups") {
      if (ncol(beta) == 1) return(unique(object$group[beta != 0]))
      else return(drop(apply(beta!=0, 2, function(x) unique(object$group[x]))))
    }
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
    if (inherits(object, "expanded")) X <- predict_spline(object, X)
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
    if (type=="groups") return(drop(apply(beta, 3, function(x) which(apply(x!=0, 2, any)))))
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

#' @rdname predict.grpreg
#' @export

coef.grpreg <- function(object, lambda, which=1:length(object$lambda), drop=TRUE, ...) {
  if (!missing(lambda)) {
    if (any(lambda > max(object$lambda) | lambda < min(object$lambda))) stop('lambda must lie within the range of the fitted coefficient path', call.=FALSE)
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
