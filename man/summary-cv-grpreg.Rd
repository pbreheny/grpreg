\name{summary.cv.grpreg}
\alias{summary.cv.grpreg}
\alias{print.summary.cv.grpreg}
\title{Summarizing inferences based on cross-validation}
\description{Summary method for \code{cv.grpreg} or \code{cv.grpsurv}
  objects}
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
  \code{summary(cvfit)} produces an object with S3 class
  \code{"summary.cv.grpreg"}.  The class has its own print method and
  contains the following list elements:
  \item{penalty}{The penalty used by \code{grpreg}/\code{grpsurv}.}
  \item{model}{The type of model: \code{"linear"}, \code{"logistic"},
    \code{"Poisson"}, \code{"Cox"}, etc.}
  \item{n}{Number of observations}
  \item{p}{Number of regression coefficients (not including the
    intercept).}
  \item{min}{The index of \code{lambda} with the smallest
    cross-validation error.}
  \item{lambda}{The sequence of \code{lambda} values used by
    \code{cv.grpreg}/\code{cv.grpsurv}.}
  \item{cve}{Cross-validation error (deviance).}
  \item{r.squared}{Proportion of variance explained by the model, as
    estimated by cross-validation.}
  \item{snr}{Signal to noise ratio, as estimated by cross-validation.}
  \item{sigma}{For linear regression models, the scale parameter
    estimate.}
  \item{pe}{For logistic regression models, the prediction error
    (misclassification error).}}
\author{Patrick Breheny}
\seealso{\code{\link{grpreg}}, \code{\link{cv.grpreg}},
  \code{\link{cv.grpsurv}}, \code{\link{plot.cv.grpreg}}}
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
