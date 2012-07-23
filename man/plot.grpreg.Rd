\name{plot.grpreg}
\alias{plot.grpreg}
\title{Plot coefficients from a "grpreg" object}
\description{Produces a plot of the coefficient paths for a fitted \code{grpreg} object.}
\usage{
\method{plot}{grpreg}(x, alpha=1, legend.loc, log.l = FALSE, ...)
}
\arguments{
  \item{x}{Fitted \code{"grpreg"} model.}
  \item{alpha}{Controls alpha-blending.  Default is alpha=1.}
  \item{legend.loc}{Where should the legend go?  If left unspecified, no
    legend is drawn.  See \code{legend} for details.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    FALSE.} 
  \item{\dots}{Other graphical parameters to \code{plot},
    \code{matlines}, or \code{legend}}
}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{grpreg}}
\examples{
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X, y, group, penalty="grLasso")
plot(fit)
myColors <- c("black", "red", "green", "blue", "yellow", "purple",
"orange", "brown")
plot(fit, legend.loc="topleft", col=myColors)
labs <- c("Mother's Age", "Mother's weight", "Race", "Smoking",
"# Premature", "Hypertension", "Uterine irritability", "# Phys. visits")
plot(fit, legend.loc="topleft", lwd=6, alpha=0.5, legend=labs)
}
\keyword{graphics}
