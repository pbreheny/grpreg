% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expand_spline.R
\name{expand_spline}
\alias{expand_spline}
\title{Expand feature matrix using basis splines}
\usage{
expand_spline(x, df = 3, degree = 3, type = c("ns", "bs"))
}
\arguments{
\item{x}{Features to be expanded (numeric matrix).}

\item{df}{Degrees of freedom (numeric; default = 3).}

\item{degree}{Degree of the piecewise polynomial (integer; default = 3 (cubic splines)).}

\item{type}{Type of spline, either B-spline (\code{"bs"}) or natural cubic spline (\code{"ns"}; default).}
}
\value{
An object of class \code{expandedMatrix} consisting of:
\itemize{
\item \code{X}: A matrix of dimension \code{nrow(x)} by \code{df*ncol(x)}
\item \code{group}: A vector of length \code{df*ncol(x)} that describes the grouping structure
\item Additional metadata on the splines, such as knot locations, required in order to evaluate spline at new feature values (e.g., for prediction)
}
}
\description{
Performs a basis expansion for many features at once, returning output that is compatible
for use with the \code{grpreg()} function. Returns an expanded matrix along with a vector
that describes its grouping.
}
\details{
\code{expand_spline()} uses the function \code{\link[splines:bs]{splines::bs()}} or \code{\link[splines:ns]{splines::ns()}} to generate a basis
matrix for each column of \code{x}. These matrices represent the spline basis for piecewise
polynomials with specified degree evaluated separately for each original column of \code{x}.
These matrices are then column-bound to form a single grouped matrix of derived features. A vector
that describes the grouping present in the resulting matrix is also generated. The resulting
object can be passed to \code{\link[=grpreg]{grpreg()}}.

This methodology was originally proposed by Ravikumar et al. (2009), who named it SPAM (SParse Additive Modeling).
}
\examples{
\dontshow{set.seed(1)}
Data <- gen_nonlinear_data(n=1000)
X <- expand_spline(Data$X)
fit <- grpreg(X, Data$y)
plot_spline(fit, "V02", lambda = 0.03)
}
\references{
\itemize{
\item Ravikumar P, Lafferty J, Liu H and Wasserman L (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society Series B}, \strong{71}: 1009-1030.
}
}
\seealso{
\code{\link[=plot_spline]{plot_spline()}} to visualize the resulting nonlinear fits
}
