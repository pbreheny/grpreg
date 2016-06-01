\name{cv.grpsurv}
\alias{cv.grpsurv}
\title{Cross-validation for grpsurv}
\description{Performs k-fold cross validation for penalized Cox
  regression models with grouped covariates over a grid of values for
  the regularization parameter lambda.}
\usage{
cv.grpsurv(X, y, group, ..., nfolds=10, seed, cv.ind, returnY=FALSE,
trace=FALSE)
}
\arguments{
  \item{X}{The design matrix, as in \code{grpsurv}.}
  \item{y}{The response matrix, as in \code{grpsurv}.}
  \item{group}{The grouping vector, as in \code{grpsurv}.}
  \item{...}{Additional arguments to \code{grpsurv}.}
  \item{nfolds}{The number of cross-validation folds.  Default is 10.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{cv.ind}{Which fold each observation belongs to.  By default the
    observations are randomly assigned by \code{cv.grpsurv}.}
  \item{returnY}{Should \code{cv.grpsurv} return the linear predictors
    from the cross-validation folds?  Default is FALSE; if TRUE, this
    will return a matrix in which the element for row i, column j is the
    fitted value for observation i from the fold in which observation i
    was excluded from the fit, at the jth value of lambda.  NOTE: The
    rows of \code{Y} are ordered by time on study, and therefore do not
    correspond to the original order of observations pased to
    \code{cv.grpsurv}.}
  \item{trace}{If set to TRUE, cv.grpsurv will inform the user of its
    progress by announcing the beginning of each CV fold.  Default is
    FALSE.}
  }
\details{
  The function calls \code{grpsurv} \code{nfolds} times, each time
  leaving out 1/\code{nfolds} of the data.  Because of the
  semiparametric nature of Cox regression, cross-validation is not
  clearly defined.  \code{cv.grpsurv} uses the approach of calculating
  the full Cox partial likelihood using the cross-validated set of
  linear predictors.  Unfortunately, using this approach there is no
  clear way (yet) of determining standard errors, so \code{cv.grpsurv},
  unlike \code{cv.grpreg}, does not provide any.

  Other approaches to cross-validation for the Cox regression model have
  been proposed; the strenghts and weaknesses of the various methods for
  penalized regression in the Cox model are not well understood.
  Because of this, the approach used by \code{cv.grpsurv} may change in
  the future as additional research is carried out.
  }
\value{
  An object with S3 class \code{"cv.grpsurv"} inheriting from
  \code{"cv.grpreg"} and containing:
  \item{cve}{The error for each value of \code{lambda}, averaged
    across the cross-validation folds.}
  \item{lambda}{The sequence of regularization parameter values along
    which the cross-validation error was calculated.}
  \item{fit}{The fitted \code{grpsurv} object for the whole data.}
  \item{min}{The index of \code{lambda} corresponding to
    \code{lambda.min}.}
  \item{lambda.min}{The value of \code{lambda} with the minimum
    cross-validation error.}
  \item{null.dev}{The cross-validated deviance for the first model along
    the grid of \code{lambda} (i.e., the cross-validated deviance for
    max(lambda), unless you have supplied your own \code{lambda}
    sequence, in which case this quantity is probably not meaningful).
    Although the actual null deviance can be calculated, it cannot be
    compared with the cross-validated deviance due to the manner in
    which deviance must be calculated for Cox regression models (see
    details).}
  }
\references{
  \itemize{
    \item Verweij PJ and van Houwelingen HC. (1993) Cross-validation in
    survival analysis.  \emph{Statistics in Medicine}, \strong{12}:
    2305-2314.
  }
}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{\link{grpsurv}},
  \code{\link{plot.cv.grpreg}},
  \code{\link{summary.cv.grpreg}}
}
\examples{
data(Lung)
X <- Lung$X
y <- Lung$y
group <- Lung$group

cvfit <- cv.grpsurv(X, y, group)
plot(cvfit)
coef(cvfit)
plot(cvfit$fit)
plot(cvfit, type="rsq")
}
