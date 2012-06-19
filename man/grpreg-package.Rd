\name{grpreg-package}
\alias{grpreg-package}
\docType{package}
\title{Regularization paths for regression models with grouped covariates}
\description{
  This package fits regularization paths for linear or logistic
  regression models penalized by the group lasso, group bridge, or group
  MCP methods.  The algorithm is based on the idea of a locally
  approximated coordinate descent, and is stable and very fast.
  }
\details{
\tabular{ll}{
Package: \tab grpreg\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2009-11-12\cr
License: \tab GPL-2\cr
}
Accepts a design matrix \code{X}, a vector of responses \code{y}, and
covariate grouping \code{group}, and produces the regularization path
over a grid of values for the tuning parameter \code{lambda}.  Also
provides methods for selecting the optimal point along the path using a
variety of information criteria and for plotting the paths.
}
\references{Breheny, P. and Huang, J. (2009) Penalized methods for
  bi-level variable selection.  Statistics and its interface, 2:
  369-380.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\keyword{models}
\keyword{regression}
\keyword{package}
\examples{
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)

fit <- grpreg(X,y,group,penalty="gMCP")
plot(fit)
select(fit,"AIC")
}
