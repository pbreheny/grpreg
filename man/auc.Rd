\name{AUC.cv.grpsurv}
\alias{AUC}
\alias{AUC.cv.grpsurv}
\title{Calculates AUC for cv.grpsurv objects}
\description{Calculates the cross-validated AUC (concordance) from a
  "cv.grpsurv" object.}
\usage{
\method{AUC}{cv.grpsurv}(obj, ...)
}
\arguments{
  \item{obj}{A \code{cv.grpsurv} object.  You must run \code{cv.grpsurv}
  with the option \code{returnY=TRUE} in order for \code{AUC} to work.}
  \item{\dots}{For S3 method compatibility.}
}
\details{
  The area under the curve (AUC), or equivalently, the concordance
  statistic (C), is calculated according to the procedure outlined in
  the reference below.  This calls the \code{survConcordance} function
  in the \code{survival} package, except the cross-validated linear
  predictors are used to guard against overfitting.  Thus, the values
  returned by \code{AUC.cv.grpsurv} will be lower than those
  you would obtain with \code{survConcordance} if you fit the full
  (unpenalized) model.
}
\references{van Houwelingen H, Putter H (2011). Dynamic Prediction in
  Clinical Survival Analysis.  CRC Press.}
\author{Patrick Breheny}
\seealso{\code{\link{cv.grpsurv}},
  \code{\link[survival]{survConcordance}}} 
\examples{
\dontshow{set.seed(1)}
data(Lung)
X <- Lung$X
y <- Lung$y
group <- Lung$group

cvfit <- cv.grpsurv(X, y, group, returnY=TRUE)
head(AUC(cvfit))
ll <- log(cvfit$fit$lambda)
plot(ll, AUC(cvfit), xlim=rev(range(ll)), lwd=3, type='l',
     xlab=expression(log(lambda)), ylab='AUC', las=1)
}
