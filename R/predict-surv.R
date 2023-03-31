#' Model predictions for grpsurv objects
#' 
#' Similar to other predict methods, this function returns predictions from a fitted `grpsurv` object.
#' 
#' Estimation of baseline survival function conditional on the estimated values of `beta` is carried out according to the method described in Chapter 4.3 of Kalbfleisch and Prentice.
#' 
#' @param object   Fitted `grpsurv` model object.
#' @param X        Matrix of values at which predictions are to be made. Not required for some `type` values.
#' @param lambda   Regularization parameter at which predictions are requested. For values of `lambda` not in the sequence of fitted models, linear interpolation is used.
#' @param which    Indices of the penalty parameter `lambda` at which predictions are required. Default: all indices. If `lambda` is specified, this will override `which`.
#' @param type     Type of prediction:
#'   * `link`: linear predictors
#'   * `response`: risk (i.e., `exp(link)`)
#'   * `survival`: the estimated survival function
#'   * `hazard`: the estimated cumulative hazard function
#'   * `median`: median survival time
#'   * The other options are all identical to their [grpreg()] counterparts
#' @param ...      Not used.
#' 
#' @return The object returned depends on type.
#' 
#' @references
#'   * Kalbfleish JD and Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd edition. Wiley.
#'   
#' @author Patrick Breheny
#' 
#' @seealso [grpsurv()]
#' 
#' @examples
#' data(Lung)
#' X <- Lung$X
#' 
#' y <- Lung$y
#' group <- Lung$group
#'  
#' fit <- grpsurv(X, y, group)
#' coef(fit, lambda=0.05)
#' head(predict(fit, X, type="link", lambda=0.05))
#' head(predict(fit, X, type="response", lambda=0.05))
#'  
#' # Survival function
#' S <- predict(fit, X[1,], type="survival", lambda=0.05)
#' S(100)
#' S <- predict(fit, X, type="survival", lambda=0.05)
#' plot(S, xlim=c(0,200))
#'  
#' # Medians
#' predict(fit, X[1,], type="median", lambda=0.05)
#' M <- predict(fit, X, type="median")
#' M[1:10, 1:10]
#'  
#' # Nonzero coefficients
#' predict(fit, type="vars", lambda=c(0.1, 0.01))
#' predict(fit, type="nvars", lambda=c(0.1, 0.01))
#' @export

predict.grpsurv <- function(object, X,
                            type=c("link", "response", "survival", "hazard", "median", "norm", 
                                   "coefficients", "vars", "nvars", "groups", "ngroups"),
                            lambda, which=1:length(object$lambda), ...) {
  type <- match.arg(type)
  if (type %in% c("norm", "coefficients", "vars", "nvars", "groups", "ngroups")) {
    return(predict.grpreg(object=object, X=X, type=type, lambda=lambda, which=which, ...))
  }
  if (!missing(lambda)) {
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    x <- ind %% 1
    beta <- (1-x)*object$beta[, l, drop=FALSE] + x*object$beta[, r, drop=FALSE]
    colnames(beta) <- round(lambda, 4)
  } else {
    beta <- object$beta[, which, drop=FALSE]
  }

  if (missing(X)) {
    eta <- matrix(0, 1, 1)
    warning('Returning "baseline" prediction; supply X for more interesting prediction')
  } else {
    if (inherits(object, "expanded")) X <- predict_spline(object, X)
    eta <- X %*% beta
  }
  if (type=='link') return(drop(eta))
  if (type=='response') return(drop(exp(eta)))
  if (!missing(lambda)) {
    W <- (1-x)*exp(object$linear.predictors)[, l, drop=FALSE] + x*exp(object$linear.predictors)[, r, drop=FALSE]
  } else {
    W <- exp(object$linear.predictors)[, which, drop=FALSE]
  }

  if (type %in% c('survival', 'hazard') & ncol(W) > 1) stop('Can only return type="survival" for a single lambda value', call.=FALSE)
  if (type %in% c('survival', 'hazard')) val <- vector('list', length(eta))
  if (type == 'median') val <- matrix(NA, nrow(eta), ncol(eta))
  for (j in 1:ncol(eta)) {
    # Estimate baseline hazard
    w <- W[,j]
    r <- rev(cumsum(rev(w)))
    a <- ifelse(object$fail, (1-w/r)^(1/w), 1)
    S0 <- c(1, cumprod(a))
    H0 <- c(0, cumsum(1-a))
    x <- c(0, object$time)
    for (i in 1:nrow(eta)) {
      S <- S0^exp(eta[i,j])
      H <- H0*exp(eta[i,j])
      if (type == 'survival') val[[i]] <- approxfun(x, S, method='constant', ties="ordered")
      else if (type == 'hazard') val[[i]] <- approxfun(x, H, method='constant', ties="ordered")
      else if (type == 'median') {
        if (any(S < 0.5)) {
          val[i,j] <- x[min(which(S < .5))]
        }
      }
    }
  }
  if (type %in% c('survival', 'hazard')) {
    if (nrow(eta)==1) val <- val[[1]]
    class(val) <- c('grpsurv.func', class(val))
    attr(val, 'time') <- object$time
  } else if (type == 'median') {
    val <- drop(val)
  }
  val
}
