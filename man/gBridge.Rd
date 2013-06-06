\name{gBridge}
\alias{gBridge}
\title{Fit a group bridge regression path}
\description{Fit regularization paths for linear and logistic group
  bridge-penalized regression models over a grid of values for the
  regularization parameter lambda.}
\usage{
gBridge(X, y, group=1:ncol(X), family=c("gaussian","binomial"),
nlambda=100, lambda, lambda.min={if (nrow(X) > ncol(X)) .001 else .05},
lambda.max, alpha=1, eps=.005, delta=1e-7, max.iter=1000, gamma=0.5,
group.multiplier=rep(1,J), warn=TRUE) 
}
\arguments{
  \item{X}{The design matrix, as in \code{grpreg}.}
  \item{y}{The response vector (or matrix), as in \code{grpreg}.}
  \item{group}{The grouping vector, as in \code{grpreg}.}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{nlambda}{The number of \code{lambda} values, as in
    \code{grpreg}.} 
  \item{lambda}{A user supplied sequence of \code{lambda} values, as in
    \code{grpreg}.} 
  \item{lambda.min}{The smallest value for \code{lambda}, as in
    \code{grpreg}.} 
  \item{lambda.max}{The maximum value for \code{lambda}.  Unlike the
    penalties in \code{grpreg}, it is not possible to solve for
    \code{lambda.max} directly with group bridge models.  Thus, it must
    be specified by the user.  If it is not specified, \code{gBridge}
    will attempt a guess at \code{lambda.max}, but this is not
    particularly accurate.}
  \item{alpha}{Tuning parameter for the balance between the group
    penalty and the L2 penalty, as in \code{grpreg}.}
  \item{eps}{Convergence threshhold, as in \code{grpreg}.}
  \item{delta}{The group bridge penalty is not differentiable at zero,
    and requires a small number \code{delta} to bound it away from
    zero.  There is typically no need to change this value.}
  \item{max.iter}{Maximum number of iterations, as in
    \code{grpreg}.}
  \item{gamma}{Tuning parameter of the group bridge penalty (the
    exponent to which the L1 norm of the coefficients in the group are
    raised).  Default is 0.5, the square root.}
  \item{group.multiplier}{The multiplicative factor by which each
    group's penalty is to be multiplied, as in \code{grpreg}.}
  \item{warn}{Should the function give a warning if it fails to
    converge?  As in \code{grpreg}.}
}
\details{
  This method fits the group bridge method of Huang et al. (2009).
  Unlike the penalties in \code{grpreg}, the group bridge is not
  differentiable at zero; because of this, a number of changes must be 
  made to the algorithm, which is why it has its own function.  Most
  notably, the method is unable to start at \code{lambda.max}; it must
  start at \code{lambda.max} and proceed in the opposite direction.

  In other respects, the usage and behavior of the function is similar
  to the rest of the \code{grpreg} package.
  }
\value{An object with S3 class \code{"grpreg"}, as in
  \code{grpreg}.}
\references{
  \itemize{
    \item Huang, J., Ma, S., Xie, H., and Zhang, C. (2009) A group bridge
    approach for variable selection.  \emph{Biometrika}, \strong{96}:
    339-355. 
    \item Breheny, P. and Huang, J. (2009) Penalized methods for
    bi-level variable selection.  \emph{Statistics and its interface},
    \strong{2}: 369-380.
    }}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{\link{grpreg}}}
\examples{
## Linear regression
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- gBridge(X, y, group)
plot(fit)
select(fit)

## Logistic regression
y <- birthwt.grpreg$low
fit <- gBridge(X, y, group, family="binomial")
plot(fit)
select(fit)
}
\keyword{models}

