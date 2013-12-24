\name{grpreg-internal}
\title{Internal grpreg functions}
\alias{calcL}
\alias{criteria.grpreg}
\alias{setupLambda}
\alias{unstandardize}
\description{Internal grpreg functions}
\usage{
setupLambda(X, y, group, family, penalty, alpha, lambda.min, nlambda,
group.multiplier) 
setupLambda.gBridge(X, y, group, family, alpha, lambda.min, lambda.max,
nlambda, gamma, group.multiplier) 
standardize(X)
unstandardize(b, center, scale)
orthogonalize(X, group)
unorthogonalize(b, XX, group.full, group)
}
\author{Patrick Breheny <patrick-breheny@uiowa.edu>}
\details{These are not intended for use by users. \code{setupLambda} and
  \code{setupLambda.gBridge} create appropriate vectors of
  regularization parameter values. The functions of \code{standardize},
  \code{unstandardize}, \code{orthogonalize}, and \code{unorthogonalize}
  are self-explanatory.
  }
\keyword{internal}
