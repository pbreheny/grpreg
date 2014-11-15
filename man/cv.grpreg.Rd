\name{cv.grpreg}
\alias{cv.grpreg}
\title{Cross-validation for grpreg}
\description{Performs k-fold cross validation for penalized regression
  models with grouped covariates over a grid of values for the
  regularization parameter lambda.}
\usage{
cv.grpreg(X, y, group=1:ncol(X), ..., nfolds=10, seed, trace=FALSE)
}
\arguments{
  \item{X}{The design matrix, as in \code{grpreg}.}
  \item{y}{The response vector (or matrix), as in \code{grpreg}.}
  \item{group}{The grouping vector, as in \code{grpreg}.}
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
  error is based on the residual sum of squares when
  \code{family="gaussian"} and the deviance when
  \code{family="binomial"} or \code{family="poisson"}.

  For Gaussian and Poisson responses, the folds are chosen according to
  simple random sampling.  For binomial responses, the numbers for each
  outcome class are balanced across the folds; i.e., the number of
  outcomes in which \code{y} is equal to 1 is the same for each fold, or
  possibly off by 1 if the numbers do not divide evenly.

  As in \code{grpreg}, seemingly unrelated regressions/multitask
  learning can be carried out by setting \code{y} to be a matrix, in
  which case groups are set up automatically (see \code{\link{grpreg}}
  for details), and cross-validation is carried out with respect to rows
  of \code{y}.  As mentioned in the details there, it is recommended to
  standardize the responses prior to fitting.
}
\value{
  An object with S3 class \code{"cv.grpreg"} containing:
  \item{cve}{The error for each value of \code{lambda}, averaged
    across the cross-validation folds.}
  \item{cvse}{The estimated standard error associated with each value of
    for \code{cve}.}
  \item{lambda}{The sequence of regularization parameter values along
    which the cross-validation error was calculated.}
  \item{fit}{The fitted \code{grpreg} object for the whole data.}
  \item{min}{The index of \code{lambda} corresponding to
    \code{lambda.min}.}
  \item{lambda.min}{The value of \code{lambda} with the minimum
    cross-validation error.}
  \item{null.dev}{The deviance for the intercept-only model.}
  \item{pe}{If \code{family="binomial"}, the cross-validation prediction
    error for each value of \code{lambda}.}
}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{grpreg}}, \code{\link{plot.cv.grpreg}},
  \code{\link{summary.cv.grpreg}}, \code{\link{predict.cv.grpreg}}}
\examples{
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)

cvfit <- cv.grpreg(X, y, group)
plot(cvfit)
summary(cvfit)
coef(cvfit) ## Beta at minimum CVE

cvfit <- cv.grpreg(X, y, group, penalty="gel")
plot(cvfit)
summary(cvfit)
}
