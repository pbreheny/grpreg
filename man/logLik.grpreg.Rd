\name{logLik.grpreg}
\alias{logLik}
\alias{logLik.grpreg}
\title{logLik method for grpreg}
\description{Calculates the log likelihood and degrees of freedom for a
  fitted grpreg object.}
\usage{
\method{logLik}{grpreg}(object, df.method=c("default","active"), ...)
}
\arguments{
  \item{object}{A fitted grpreg object.}
  \item{df.method}{How should effective model parameters be calculated?
    One of: \code{"active"}, which counts the number of nonzero
    coefficients; or \code{"default"}, which uses the calculated
    \code{df} returned by \code{grpreg}.  Default is \code{"default"}.}
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
\references{Breheny, P. and Huang, J. (2009) Penalized methods for
  bi-level variable selection.  Statistics and its interface, 2:
  369-380.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{grpreg}}
\examples{
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X,y,group,penalty="gMCP")
logLik(fit) ## Display is glitchy
AIC(fit) ## Works fine
\dontrun{BIC(fit) ## Works on R > 2.13-0, or if certain packages are loaded}
}
\keyword{models}
\keyword{regression}
