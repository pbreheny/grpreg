% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_nonlinear_data.R
\name{gen_nonlinear_data}
\alias{gen_nonlinear_data}
\title{Generate nonlinear example data}
\usage{
gen_nonlinear_data(n = 100, p = 16, seed)
}
\arguments{
\item{n}{Sample size (numeric; default = 100).}

\item{p}{Number of features (numeric; default = 16).}

\item{seed}{Set to get different random data sets, passed to \code{\link[=set.seed]{set.seed()}}}
}
\description{
Mainly intended to demonstrate the use of basis expansion models for sparse additive modeling; intended for use with \code{\link[=expand_spline]{expand_spline()}}.
}
\examples{
Data <- gen_nonlinear_data()
}
