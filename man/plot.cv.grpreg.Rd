\name{plot.cv.grpreg}
\alias{plot.cv.grpreg}
\title{Plots the cross-validation curve from a \code{cv.grpreg} object}
\description{Plots the cross-validation curve from a \code{cv.grpreg}
  object, along with standard error bars.}
\usage{
\method{plot}{cv.grpreg}(x, log.l=TRUE, ...)
}
\arguments{
  \item{x}{A \code{cv.grpreg} object.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    TRUE.} 
  \item{\dots}{Other graphical parameters to \code{plot}}
}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{grpreg}, \code{cv.grpreg}}
\examples{
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
cv.fit <- cv.grpreg(X, y, group, penalty="grLasso")
plot(cv.fit)
cv.fit <- cv.grpreg(X, y, group, penalty="gMCP")
plot(cv.fit)
}
\keyword{graphics}
