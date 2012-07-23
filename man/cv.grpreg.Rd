\name{cv.grpreg}
\alias{cv.grpreg}
\title{Cross-validation for grpreg}
\description{Performs k-fold cross validation for penalized regression
  models with grouped covariates over a grid of values for the
  regularization parameter lambda.}
\usage{
cv.grpreg(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD",
"gMCP", "gLasso"), family=c("gaussian","binomial"), nlambda=100, lambda,
lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, alpha=1, nfolds=10,
seed, trace=FALSE, group.multiplier=rep(1,J),...)
}
\arguments{
  \item{X}{The design matrix, as in \code{grpreg}.}
  \item{y}{The response vector.}
  \item{group}{The grouping vector, as in \code{grpreg}.}
  \item{penalty}{The penalty to be applied to the model, as in
    \code{grpreg}.}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{nlambda}{The number of \code{lambda} values, as in
    \code{grpreg}.} 
  \item{lambda}{A user supplied sequence of \code{lambda} values, as in
    \code{grpreg}.} 
  \item{lambda.min}{The smallest value for \code{lambda}, as in
    \code{grpreg}.} 
  \item{alpha}{Tuning parameter for the balance between the group
    penalty and the L2 penalty, as in \code{grpreg}.}
  \item{nfolds}{The number of cross-validation folds.  Default is 10.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{trace}{If set to TRUE, cv.grpreg will inform the user of its
    progress by announcing the beginning of each CV fold.  Default is
    FALSE.}
  \item{group.multiplier}{The multiplicative factor by which each
    group's penalty is to be multiplied, as in \code{grpreg}.}
  \item{...}{Additional arguments to \code{grpreg}.}}
\details{
  The function calls \code{grpreg} \code{nfolds} times, each time
  leaving out 1/\code{nfolds} of the data.  The cross-validation
  error is based on the loss function used to fit the model; see
  \code{grpreg} for details.}
\value{
  An object with S3 class \code{"cv.grpreg"} containing:
  \item{E}{A matrix with \code{nfolds} rows and \code{nlambda}
    columns, containing the cross-validation error for each fit.}
  \item{cve}{The error for each value of \code{lambda}, averaged
    across the cross-validation folds.}
  \item{lambda}{The sequence of regularization parameter values along
    which the cross-validation error was calculated.}
  \item{lambda.min}{The value of \code{lambda} with the minimum
    cross-validation error.}
  \item{min}{The index of \code{lambda} corresponding
    \code{lambda.min}.}
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
