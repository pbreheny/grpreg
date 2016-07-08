\name{predict.grpreg}
\alias{predict.grpreg}
\alias{coef.grpreg}
\alias{predict.cv.grpreg}
\alias{coef.cv.grpreg}
\title{Model predictions based on a fitted \code{grpreg} object}
\description{Similar to other predict methods, this function returns
  predictions from a fitted \code{"grpreg"} object.}
\usage{
\method{predict}{grpreg}(object, X, type=c("link", "response", "class",
"coefficients", "vars", "groups", "nvars", "ngroups", "norm"), lambda,
which=1:length(object$lambda), ...) 
\method{coef}{grpreg}(object, lambda, which=1:length(object$lambda),
drop=TRUE, ...)
\method{predict}{cv.grpreg}(object, X, lambda=object$lambda.min,
which=object$min, type=c("link", "response", "class", "coefficients",
"vars", "groups", "nvars", "ngroups", "norm"), ...) 
\method{coef}{cv.grpreg}(object, lambda=object$lambda.min,
which=object$min, ...)
}
\arguments{
  \item{object}{Fitted \code{"grpreg"} or \code{"cv.grpreg"} model object.}
  \item{X}{Matrix of values at which predictions are to be made.  Not
    used for \code{type="coefficients"}}
  \item{lambda}{Values of the regularization parameter \code{lambda} at
    which predictions are requested.  For values of \code{lambda} not in
    the sequence of fitted models, linear interpolation is used.}
  \item{which}{Indices of the penalty parameter \code{lambda} at which
    predictions are required.  By default, all indices are returned.  If
    \code{lambda} is specified, this will override \code{which}.}
  \item{type}{Type of prediction: \code{"link"} returns the
    linear predictors; \code{"response"} gives the fitted values;
    \code{"class"} returns the binomial outcome with the highest
    probability; \code{"coefficients"} returns the coefficients;
    \code{"vars"} returns the indices for the nonzero coefficients;
    \code{"groups"} returns the indices for the groups with at least one
    nonzero coefficient; \code{"nvars"} returns the number of nonzero
    coefficients; \code{"ngroups"} returns the number of groups with at
    least one nonzero coefficient; \code{"norm"} returns the L2 norm of
    the coefficients in each group.}
  \item{drop}{By default, if a single value of \code{lambda} is
    supplied, a vector of coefficients is returned.  Set
    \code{drop=FALSE} if you wish to have \code{coef} always return a
    matrix (see \code{\link{drop}}).}
  \item{\dots}{Not used.}}
\details{
  \code{coef} and \code{predict} methods are provided for
  \code{"cv.grpreg"} options as a convenience.  They simply call
  \code{coef.grpreg} and \code{predict.grpreg} with \code{lambda} set to
  the value that minimizes the cross-validation error.}
\value{The object returned depends on type.}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{grpreg}}
\examples{
# Fit penalized logistic regression model to birthweight data
data(Birthwt)
X <- Birthwt$X
y <- Birthwt$low
group <- Birthwt$group
fit <- grpreg(X, y, group, penalty="grLasso", family="binomial")

# Coef and predict methods
coef(fit, lambda=.001)
predict(fit, X, type="link", lambda=.001)
predict(fit, X, type="response", lambda=.001)
predict(fit, X, type="class", lambda=.001)
predict(fit, type="vars", lambda=.07)
predict(fit, type="groups", lambda=.07)
predict(fit, type="norm", lambda=.07)

# Coef and predict methods for cross-validation
cvfit <- cv.grpreg(X, y, group, family="binomial", penalty="grMCP")
coef(cvfit)
predict(cvfit, X)
predict(cvfit, X, type="response")
predict(cvfit, type="groups")
}
