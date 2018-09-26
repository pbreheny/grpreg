\name{plot.cv.grpreg}
\alias{plot.cv.grpreg}
\title{Plots the cross-validation curve from a \code{cv.grpreg} object}
\description{Plots the cross-validation curve from a \code{cv.grpreg}
  object, along with standard error bars.}
\usage{
\method{plot}{cv.grpreg}(x, log.l=TRUE, type=c("cve", "rsq", "scale",
"snr", "pred", "all"),  selected=TRUE, vertical.line=TRUE, col="red",
...)}
\arguments{
  \item{x}{A \code{cv.grpreg} object.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    TRUE.}
  \item{type}{What to plot on the vertical axis.  \code{cve} plots the
    cross-validation error (deviance); \code{rsq} plots an estimate of
    the fraction of the deviance explained by the model (R-squared);
    \code{snr} plots an estimate of the signal-to-noise ratio;
    \code{scale} plots, for \code{family="gaussian"}, an estimate of the
    scale parameter (standard deviation); \code{pred} plots, for
    \code{family="binomial"}, the estimated prediction error; \code{all}
    produces all of the above.}
  \item{selected}{If \code{TRUE} (the default), places an axis on top of
    the plot denoting the number of groups in the model (i.e., that
    contain a nonzero regression coefficient) at that value of
    \code{lambda}.}  
  \item{vertical.line}{If \code{TRUE} (the default), draws a vertical
    line at the value where cross-validaton error is minimized.}
  \item{col}{Controls the color of the dots (CV estimates).}
  \item{\dots}{Other graphical parameters to \code{plot}}
}
\details{
  Error bars representing approximate +/- 1 SE (68\% confidence
  intervals) are plotted along with the estimates at value of
  \code{lambda}.  For \code{rsq} and \code{snr}, these confidence
  intervals are quite crude, especially near zero, and will hopefully be
  improved upon in later versions of \code{grpreg}.}
\author{Patrick Breheny}
\seealso{\code{\link{grpreg}}, \code{\link{cv.grpreg}}}
\examples{
# Birthweight data
data(Birthwt)
X <- Birthwt$X
group <- Birthwt$group

# Linear regression
y <- Birthwt$bwt
cvfit <- cv.grpreg(X, y, group)
plot(cvfit)
par(mfrow=c(2,2))
plot(cvfit, type="all")

## Logistic regression
y <- Birthwt$low
cvfit <- cv.grpreg(X, y, group, family="binomial")
plot(cvfit)
par(mfrow=c(2,2))
plot(cvfit, type="all")
}
