\name{grpreg}
\alias{grpreg}
\title{Fit a group penalized regression path}
\description{Fit regularization paths for models with grouped penalties
  over a grid of values for the regularization parameter lambda. Fits
  linear and logistic regression models.}
\usage{
grpreg(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", 
"gMCP", "gBridge", "gLasso"), family=c("gaussian","binomial"), 
nlambda=100, lambda, lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, 
alpha=1, eps=.005, max.iter=1000, dfmax=p, gamma = 3, 
group.multiplier={if (strtrim(penalty,2)=="gr") 
sqrt(table(group[group!=0])) else rep(1,J)}, warn=TRUE, ...) 
}
\arguments{
  \item{X}{The design matrix, without an intercept.  \code{grpreg}
    standardizes the data and includes an intercept by default.}
  \item{y}{The response vector, or a matrix in the case of multitask
    learning (see details).}
  \item{group}{A vector of consecutive integers describing the grouping
    of the coefficients.  Integers for penalized groups must run from 1
    to the number of penalized groups (see example below).  If there are
    coefficients to be included in the model without being penalized,
    assign them to group 0.}
  \item{penalty}{The penalty to be applied to the model.  One of
    \code{grLasso}, \code{grMCP}, \code{grSCAD}, or \code{gMCP}.  See
    below for details.}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{nlambda}{The number of \code{lambda} values.  Default is 100.}
  \item{lambda}{A user supplied sequence of \code{lambda} values.
    Typically, this is left unspecified, and the function automatically
    computes a grid of lambda values that ranges uniformly on the log 
    scale over the relevant range of lambda values.}
  \item{lambda.min}{The smallest value for \code{lambda}, as a fraction
    of \code{lambda.max}.  Default is .0001 if the number of observations
    is larger than the number of covariates and .05 otherwise.} 
  \item{alpha}{\code{grpreg} allows for both a group penalty and an L2
    (ridge) penalty; \code{alpha} controls the proportional weight of
    the regularization parameters of these two penalties.  The group
    penalties' regularization parameter is \code{lambda*alpha}, while
    the regularization parameter of the ridge penalty is
    \code{lambda*(1-alpha)}.  Default is 1: no ridge penalty.}
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{.005}.  See details.}
  \item{max.iter}{Maximum number of iterations.  Default is 1000.  See
    details.}
  \item{dfmax}{Limit on the number of parameters allowed to be nonzero.
    If this limit is exceeded, the algorithm will exit early from the
    regularization path.}
  \item{gamma}{Tuning parameter of the MCP penalty; defaults to 3.}
  \item{group.multiplier}{A vector of values representing multiplicative
    factors by which each group's penalty is to be multiplied.
    Often, this is a function (such as the square root) of the number of
    predictors in each group.  The default is a vector of 1's (i.e., no
    correction for group size).}
  \item{warn}{Should the function give a warning if it fails to
    converge?  Default is TRUE.  See details.}
  \item{...}{Arguments passed to other functions (such as gBridge).}
}
\details{
  There are two general classes of methods involving grouped penalties:
  those that carry out bi-level selection and those that carry out group
  selection.  Bi-level means  carrying out variable selection at the group
  level as well as the level of individual covariates (i.e., selecting
  important groups as well as important members of those groups).  Group
  selection selects important groups, not members within the group --
  i.e., within a group, coefficients will either all be zero or all
  nonzero.  The \code{grLasso}, \code{grMCP}, and \code{grSCAD}
  penalties carry out group selection, while the \code{gMCP} penalty
  carries out bi-level selection.  For bi-level selection, see also the
  \code{\link{gBridge}} function.

  Please note the distinction between \code{grMCP} and \code{gMCP},
  which have unfortunately similar names.  The former involves an MCP
  penalty being applied to an L2-norm of each group.  The latter
  involves a hierarchical penalty which places an outer MCP penalty on a
  sum of inner MCP penalties for each group, as proposed in Breheny &
  Huang, 2009.  Either penalty may be referred to as the "group MCP",
  depending on the publication.  The \code{gMCP} method is referred to
  as "composite MCP" in Huang et al. (2012).

  In keeping with the notation from the original MCP paper, the tuning
  parameter of the MCP penalty is denoted 'gamma'.  Note, however, that
  in Breheny and Huang (2009), \code{gamma} is denoted 'a'.

  The objective function is defined to be
  \deqn{\frac{1}{2n}RSS + penalty}{RSS/(2*n) + penalty}
  for \code{"gaussian"} and
  \deqn{-\frac{1}{n} loglik + penalty}{-loglik/n + penalty}
  for \code{"binomial"}, where the likelihood is from a traditional
  generalized linear model for the log-odds of an event.  For logistic
  regression models, some care is taken to avoid model saturation; the
  algorithm  may exit early in this setting. 
  
  For the bi-level selection methods, a locally approximated coordinate
  descent algorithm is employed.  For the group selection methods, group
  descent algorithms are employed.
  
  The algorithms employed by \code{grpreg} are stable and generally
  converge quite rapidly to values close to the solution.  However,
  especially when p is large compared with n, \code{grpreg} may fail to
  converge at low values of \code{lambda}, where models are
  nonidentifiable or nearly singular.  Often, this is not the region of
  the coefficient path that is most interesting.  The default behavior
  warning the user when convergence criteria are not met may be
  distracting in these cases, and can be modified with \code{warn}
  (convergence can always be checked later by inspecting the value of
  \code{iter}).

  If models are not converging, increasing \code{max.iter} may not be
  the most efficient way to correct this problem.  Consider increasing
  \code{n.lambda} or \code{lambda.min} in addition to increasing
  \code{max.iter}.

  Seemingly unrelated regressions/multitask learning can be carried out
  using \code{grpreg} by passing a matrix to \code{y}.  In this case,
  \code{X} will be used in separate regressions for each column of
  \code{y}, with the coefficients grouped across the responses.  In
  other words, each column of \code{X} will form a group with m
  members, where m is the number of columns of \code{y}.  For multiple
  Gaussian responses, it is recommended to standardize the columns of
  \code{y} prior to fitting, in order to apply the penalization equally
  across columns.
}
\value{
  An object with S3 class \code{"grpreg"} containing:
  \item{beta}{The fitted matrix of coefficients.  The number of rows is
    equal to the number of coefficients, and the number of columns is
    equal to \code{nlambda}.}
  \item{family}{Same as above.}
  \item{group}{Same as above.}
  \item{lambda}{The sequence of \code{lambda} values in the path.}
  \item{alpha}{Same as above.}
  \item{loss}{A vector containing either the residual sum of squares
  (\code{"gaussian"}) or negative log-likelihood (\code{"binomial"}) of
  the fitted model at each value of \code{lambda}.}
  \item{n}{Number of observations.}
  \item{penalty}{Same as above.}
  \item{df}{A vector of length \code{nlambda} containing estimates of
    effective number of model parameters all the points along the
    regularization path.  For details on how this is calculated, see
    Breheny and Huang (2009).}
  \item{iter}{A vector of length \code{nlambda} containing the number
    of iterations until convergence at each value of \code{lambda}.}
}
\references{
  \itemize{
    \item Breheny, P. and Huang, J. (2009) Penalized methods for
    bi-level variable selection.  \emph{Statistics and its interface},
    \strong{2}: 369-380.
    \item Huang J., Breheny, P. and Ma, S. (2012). A selective
    review of group selection in high dimensional
    models. \emph{Statistical Science}, \strong{27}: 481-499.}
  }
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{\link{cv.grpreg}}, as well as \code{\link[=plot.grpreg]{plot}} and \code{\link[=select.grpreg]{select}} methods.}
\examples{
## Linear regression
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X,y,group,penalty="grLasso")
plot(fit)
fit <- grpreg(X,y,group,penalty="grMCP")
plot(fit)
fit <- grpreg(X,y,group,penalty="grSCAD")
plot(fit)
fit <- grpreg(X,y,group,penalty="gMCP")
plot(fit)
select(fit,"AIC")

## Logistic regression
y <- birthwt.grpreg$low
fit <- grpreg(X,y,group,penalty="grLasso", family="binomial")
plot(fit)
fit <- grpreg(X,y,group,penalty="grMCP", family="binomial")
plot(fit)
fit <- grpreg(X,y,group,penalty="grSCAD", family="binomial")
plot(fit)
fit <- grpreg(X,y,group,penalty="gMCP", family="binomial")
plot(fit)
select(fit,"BIC")
}
\keyword{models}

