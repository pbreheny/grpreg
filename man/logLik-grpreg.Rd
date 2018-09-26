\name{logLik.grpreg}
\alias{logLik}
\alias{logLik.grpreg}
\title{logLik method for grpreg}
\description{Calculates the log likelihood and degrees of freedom for a
  fitted grpreg object.}
\usage{
\method{logLik}{grpreg}(object, df.method=c("default","active"),
REML=FALSE, ...)
}
\arguments{
  \item{object}{A fitted grpreg object.}
  \item{df.method}{How should effective model parameters be calculated?
    One of: \code{"active"}, which counts the number of nonzero
    coefficients; or \code{"default"}, which uses the calculated
    \code{df} returned by \code{grpreg}.  Default is \code{"default"}.}
  \item{REML}{Use restricted MLE for estimation of the scale parameter
    in a gaussian model?  Default is FALSE.}
  \item{\dots}{For S3 method compatibility.}}
\details{
  Exists mainly for use with \code{'AIC'} and \code{'BIC'}.
  }
\value{
  Returns an object of class 'logLik', in this case consisting of a
  number (or vector of numbers) with two attributes: 'df' (the estimated
  degrees of freedom in the model) and 'nobs' (number of observations).

  The 'print' method for 'logLik' objects is not intended to handle
  vectors; consequently, the value of the function does not necessarily
  display correctly.  However, it works with 'AIC' and 'BIC' without any
  glitches and returns the expected vectorized output.
  }
\author{Patrick Breheny}
\seealso{\code{grpreg}}
\examples{
data(Birthwt)
X <- Birthwt$X
y <- Birthwt$bwt
group <- Birthwt$group
fit <- grpreg(X,y,group,penalty="cMCP")
logLik(fit) ## Display is glitchy for vectors
AIC(fit)
BIC(fit)
}
