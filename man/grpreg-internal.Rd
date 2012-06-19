\name{grpreg-internal}
\title{Internal grpreg functions}
\alias{calcL}
\alias{criteria.grpreg}
\alias{setupLambda}
\alias{unstandardize}
\description{Internal grpreg functions}
\usage{
calcL(X,y,beta,family)
setupLambda(X,y,group,family,penalty,lambda.max,lambda.min,nlambda,gamma,group.multiplier)
unstandardize(beta,meanx,normx)
}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\details{These are not intended for use by users. \code{calcL}
  calculates the loss function of a model given a formatted Data object
  and a vector of coefficients beta.  \code{setupLambda} creates an
  appropriate vector of regularization parameter values.
  \code{unstandardize} returns standardized coefficients to
  their original scale.
  }
\keyword{internal}
