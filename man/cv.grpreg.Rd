\name{cv.grpreg}
\alias{cv.grpreg}
\title{Cross-validation for grpreg}
\description{Performs k-fold cross validation for penalized regression
  models with grouped covariates over a grid of values for the
  regularization parameter lambda.}
\usage{
cv.grpreg(X, y, ..., nfolds=10, seed, trace=FALSE)
}
\arguments{
  \item{X}{The design matrix, as in \code{grpreg}.}
  \item{y}{The response vector, as in \code{grpreg}.}
  \item{...}{Additional arguments to \code{grpreg}.}
  \item{nfolds}{The number of cross-validation folds.  Default is 10.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{trace}{If set to TRUE, cv.grpreg will inform the user of its
    progress by announcing the beginning of each CV fold.  Default is
    FALSE.}
  }
\details{
  The function calls \code{grpreg} \code{nfolds} times, each time
  leaving out 1/\code{nfolds} of the data.  The cross-validation
  error is based on the loss function used to fit the model; see
  \code{grpreg} for details.}
\value{
  An object with S3 class \code{"cv.grpreg"} containing:
  \item{cve}{The error for each value of \code{lambda}, averaged
    across the cross-validation folds.}
  \item{cvse}{The estimated standard error associated with each value of
    for \code{cve}.}
  \item{lambda}{The sequence of regularization parameter values along
    which the cross-validation error was calculated.}
  \item{fit}{The fitted \code{grpreg} object for the whole data.}
  \item{min}{The index of \code{lambda} corresponding
    \code{lambda.min}.}
  \item{lambda.min}{The value of \code{lambda} with the minimum
    cross-validation error.}
  }
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{\link{grpreg}}}
\examples{
## Linear regression
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)

fit <- cv.grpreg(X, y, group, penalty="grLasso")
cv.fit <- cv.grpreg(X, y, group, penalty="grLasso")
plot(cv.fit)
coef(fit, cv.fit$lambda.min)
}
