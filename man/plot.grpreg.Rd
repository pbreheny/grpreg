\name{plot.grpreg}
\alias{plot.grpreg}
\title{Plot coefficients from a "grpreg" object}
\description{Produces a plot of the coefficient paths for a fitted \code{grpreg} object.}
\usage{
\method{plot}{grpreg}(x, alpha=1, legend.loc, label=FALSE, log.l=FALSE,
norm=FALSE, ...)
}
\arguments{
  \item{x}{Fitted \code{"grpreg"} model.}
  \item{alpha}{Controls alpha-blending.  Default is alpha=1.}
  \item{legend.loc}{Where should the legend go?  If left unspecified, no
    legend is drawn.  See \code{\link[graphics]{legend}} for details.}
  \item{label}{If TRUE, annotates the plot with text labels in the right
    margin describing which variable/group the corresponding line
    belongs to.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    FALSE.}
  \item{norm}{If \code{TRUE}, plot the norm of each group, rather than
    the individual coefficients.}
  \item{\dots}{Other graphical parameters to \code{plot},
    \code{matlines}, or \code{legend}}
}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\seealso{\code{grpreg}}
\examples{
# Fit model to birthweight data
data(Birthwt)
X <- Birthwt$X
y <- Birthwt$bwt
group <- Birthwt$group
fit <- grpreg(X, y, group, penalty="grLasso")

# Plot (basic)
plot(fit)

# Plot group norms, with labels in right margin
plot(fit, norm=TRUE, label=TRUE)

# Plot (miscellaneous options)
myColors <- c("black", "red", "green", "blue", "yellow", "purple",
"orange", "brown")
plot(fit, legend.loc="topleft", col=myColors)
labs <- c("Mother's Age", "# Phys. visits", "Hypertension", "Mother's weight",
          "# Premature", "Race", "Smoking", "Uterine irritability")
plot(fit, legend.loc="topleft", lwd=6, alpha=0.5, legend=labs)
plot(fit, norm=TRUE, legend.loc="topleft", lwd=6, alpha=0.5, legend=labs)
}
