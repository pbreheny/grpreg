\name{predict.grpreg}
\alias{predict.grpreg}
\alias{coef.grpreg}
\title{Model predictions based on a fitted \code{grpreg} object}
\description{Similar to other predict methods, this function returns
  predictions from a fitted \code{"grpreg"} object.}
\usage{
\method{predict}{grpreg}(object, X, lambda, which=1:length(object$lambda),
type=c("link", "response", "class", "coefficients", "vars", "groups",
"norm"), ...)
\method{coef}{grpreg}(object, lambda, which=1:length(object$lambda),...)
}
\arguments{
  \item{object}{Fitted \code{"grpreg"} model object.}
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
    nonzero coefficient; \code{"norm"} returns the L2 norm of the
    coefficients in each group.}
  \item{\dots}{Not used.}}
\value{The object returned depends on type.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{grpreg}}
\examples{
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$low
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X, y, group, penalty="grLasso", family="binomial")
coef(fit, lambda=.001)
predict(fit, X, type="link", lambda=.001)
predict(fit, X, type="response", lambda=.001)
predict(fit, X, type="class", lambda=.001)
predict(fit, X, type="vars", lambda=.07)
predict(fit, X, type="groups", lambda=.07)
predict(fit, X, type="norm", lambda=.07)
}
