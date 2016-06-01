\name{summary.cv.grpreg}
\alias{summary.cv.grpreg}
\alias{print.summary.cv.grpreg}
\title{Summarizing inferences based on cross-validation}
\description{Summary method for \code{cv.grpreg} objects}
\usage{
\method{summary}{cv.grpreg}(object, ...)
\method{print}{summary.cv.grpreg}(x, digits, ...)
}
\arguments{
  \item{object}{A \code{"cv.grpreg"} object.}
  \item{x}{A \code{"summary.cv.grpreg"} object.}
  \item{digits}{Number of digits past the decimal point to print out.
    Can be a vector specifying different display digits for each of the
    five non-integer printed values.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\value{
  \code{summary.cv.grpreg} produces an object with S3 class
  \code{"summary.cv.grpreg"}.  The class has its own print method and
  contains the following list elements:
  \item{penalty}{The penalty used by \code{grpreg}.}
  \item{model}{Either \code{"linear"} or \code{"logistic"}, depending on
    the \code{family} option in \code{grpreg}.}
  \item{n}{Number of observations}
  \item{p}{Number of regression coefficients (not including the
    intercept).}
  \item{min}{The index of \code{lambda} with the smallest
    cross-validation error.}
  \item{lambda}{The sequence of \code{lambda} values used by
    \code{cv.grpreg}.}
  \item{cve}{Cross-validation error (deviance).}
  \item{r.squared}{Proportion of variance explained by the model, as
    estimated by cross-validation.}
  \item{snr}{Signal to noise ratio, as estimated by cross-validation.}
  \item{sigma}{For linear regression models, the scale parameter
    estimate.}
  \item{pe}{For logistic regression models, the prediction error
    (misclassification error).}}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{grpreg}}, \code{\link{cv.grpreg}}, \code{\link{plot.cv.grpreg}}}
\examples{
# Birthweight data
data(Birthwt)
X <- Birthwt$X
group <- Birthwt$group

# Linear regression
y <- Birthwt$bwt
cvfit <- cv.grpreg(X, y, group)
summary(cvfit)

# Logistic regression
y <- Birthwt$low
cvfit <- cv.grpreg(X, y, group, family="binomial")
summary(cvfit)

# Cox regression
data(Lung)
cvfit <- with(Lung, cv.grpsurv(X, y, group))
summary(cvfit)
}
