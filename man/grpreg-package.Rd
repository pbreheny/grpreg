\name{grpreg-package}
\alias{grpreg-package}
\docType{package}
\title{Regularization paths for regression models with grouped covariates}
\description{
  This package fits regularization paths for linear or logistic
  regression models with grouped penalties, such as the group lasso,
  group MCP, group SCAD, and group bridge.  The algorithms are based on
  the idea of either locally approximated coordinate descent or group
  descent, depending on the penalty.  All of the algorithms (with the
  exception of group bridge) are stable and fast.
}
\details{The following penalties are available:
  \itemize{
    \item \code{grLasso}: Group lasso (Yuan and Lin, 2006)
    \item \code{grMCP}: Group MCP; like the group lasso, but with an MCP
    penalty on the norm of each group
    \item \code{grSCAD}: Group SCAD; like the group lasso, but with a SCAD
    penalty on the norm of each group
    \item \code{gMCP}: A hierarchical penalty which places an outer MCP
    penalty on a sum of inner MCP penalties for each group (Breheny &
    Huang, 2009)
    \item \code{gBridge}: A penalty which places a bridge penalty on the
    L1-norm of each group (Huang et al., 2009)}
  The \code{gMCP} and \code{gBridge} penalties carry out bi-level
  selection, meaning that they carry out variable selection at the group
  level and at the level of individual covariates (i.e., they select
  important groups as well as important members of those groups).  The
  \code{grLasso}, \code{grMCP}, and \code{grSCAD} penalties carry out
  group selection, meaning that within a group, coefficients will either
  all be zero or all nonzero.

  A variety of supporting methods for selecting lambda and plotting the
  paths are provided also.

  See further details at \code{\link{grpreg}} and \code{\link{gBridge}}.

  The following functions are documented in the help pages:
  \itemize{
    \item \code{\link{grpreg}}
    \item \code{\link{gBridge}}
    \item \code{\link{cv.grpreg}}
    \item \code{\link{logLik.grpreg}}
    \item \code{\link{plot.grpreg}}
    \item \code{\link{plot.cv.grpreg}}
    \item \code{\link{predict.grpreg}}
    \item \code{\link{select.grpreg}}
  }

  The following data sets are documented in the help pages:
  \itemize{
    \item \code{\link{birthwt.grpreg}}
    }
}
\references{
  \itemize{
    \item Yuan, M. and Lin, Y. (2006) Model selection and estimation in
    regression with grouped variables.  \emph{Journal of the Royal Statistical
    Society Series B}, \strong{68}: 49-67. 
    \item Huang, J., Ma, S., Xie, H., and Zhang, C. (2009) A group bridge
    approach for variable selection.  \emph{Biometrika}, \strong{96}:
    339-355. 
    \item Breheny, P. and Huang, J. (2009) Penalized methods for
    bi-level variable selection.  \emph{Statistics and its interface},
    \strong{2}: 369-380.
    \item Huang J., Breheny, P. and Ma, S. (2012). A selective
    review of group selection in high dimensional
    models. \emph{Statistical Science}, \strong{27}: 481-499.}
  }
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\examples{## See examples in grpreg}
