\name{grpreg}
\alias{grpreg}
\title{Fit a group penalized regression path}
\description{Fit regularization paths for models with grouped penalties
  over a grid of values for the regularization parameter lambda. Fits
  linear and logistic regression models.}
\usage{
grpreg(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD",
"gel", "cMCP"), family=c("gaussian", "binomial", "poisson"),
nlambda=100, lambda, lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05},
log.lambda = TRUE, alpha=1, eps=1e-4, max.iter=10000, dfmax=p,
gmax=length(unique(group)), gamma=ifelse(penalty == "grSCAD", 4, 3),
tau = 1/3, group.multiplier, warn=TRUE, returnX = FALSE, ...)
}

\arguments{
  \item{X}{The design matrix, without an intercept.  \code{grpreg}
    standardizes the data and includes an intercept by default.}
  \item{y}{The response vector, or a matrix in the case of multitask
    learning (see details).}
  \item{group}{A vector describing the grouping of the coefficients.
    For greatest efficiency and least ambiguity (see details), it is
    best if \code{group} is a factor or vector of consecutive integers,
    although unordered groups and character vectors are also allowed.
    If there are coefficients to be included in the model without being
    penalized, assign them to group 0 (or \code{"0"}).}
  \item{penalty}{The penalty to be applied to the model.  For group
    selection, one of \code{grLasso}, \code{grMCP}, or \code{grSCAD}.
    For bi-level selection, one of \code{gel} or \code{cMCP}.  See
    below for details.}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{nlambda}{The number of \code{lambda} values.  Default is 100.}
  \item{lambda}{A user supplied sequence of \code{lambda} values.
    Typically, this is left unspecified, and the function automatically
    computes a grid of lambda values that ranges uniformly on the log
    scale over the relevant range of lambda values.}
  \item{lambda.min}{The smallest value for \code{lambda}, as a fraction
    of \code{lambda.max}.  Default is .0001 if the number of
    observations is larger than the number of covariates and .05
    otherwise.}
  \item{log.lambda}{Whether compute the grid values of lambda on log
    scale (default) or linear scale.}
  \item{alpha}{\code{grpreg} allows for both a group penalty and an L2
    (ridge) penalty; \code{alpha} controls the proportional weight of
    the regularization parameters of these two penalties.  The group
    penalties' regularization parameter is \code{lambda*alpha}, while
    the regularization parameter of the ridge penalty is
    \code{lambda*(1-alpha)}.  Default is 1: no ridge penalty.}
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    RMSD for the change in linear predictors for each coefficient is
    less than \code{eps}.  Default is \code{1e-4}.  See details.}
  \item{max.iter}{Maximum number of iterations (total across entire
    path). Default is 10000.  See details.}
  \item{dfmax}{Limit on the number of parameters allowed to be nonzero.
    If this limit is exceeded, the algorithm will exit early from the
    regularization path.}
  \item{gmax}{Limit on the number of groups allowed to have nonzero
    elements.  If this limit is exceeded, the algorithm will exit early
    from the regularization path.}
  \item{gamma}{Tuning parameter of the group or composite MCP/SCAD
    penalty (see details).  Default is 3 for MCP and 4 for SCAD.}
  \item{tau}{Tuning parameter for the group exponential lasso; defaults
    to 1/3.}
  \item{group.multiplier}{A vector of values representing multiplicative
    factors by which each group's penalty is to be multiplied.
    Often, this is a function (such as the square root) of the number of
    predictors in each group.  The default is to use the square root of
    group size for the group selection methods, and a vector of 1's
    (i.e., no adjustment for group size) for bi-level selection.}
  \item{warn}{Should the function give a warning if it fails to
    converge?  Default is TRUE.  See details.}
  \item{returnX}{Return the standardized design matrix (and associated
    group structure information)?  Default is FALSE.}
  \item{...}{Arguments passed to other functions (such as gBridge).}
}
\details{
  There are two general classes of methods involving grouped penalties:
  those that carry out bi-level selection and those that carry out group
  selection.  Bi-level means carrying out variable selection at the
  group level as well as the level of individual covariates (i.e.,
  selecting important groups as well as important members of those
  groups).  Group selection selects important groups, and not members
  within the group -- i.e., within a group, coefficients will either all
  be zero or all nonzero.  The \code{grLasso}, \code{grMCP}, and
  \code{grSCAD} penalties carry out group selection, while the
  \code{gel} and \code{cMCP} penalties carry out bi-level selection.
  For bi-level selection, see also the \code{\link{gBridge}} function.
  For historical reasons and backwards compatibility, some of these
  penalties have aliases; e.g., \code{gLasso} will do the same thing as
  \code{grLasso}, but users are encouraged to use \code{grLasso}.

  Please note the distinction between \code{grMCP} and \code{cMCP}.  The
  former involves an MCP penalty being applied to an L2-norm of each
  group.  The latter involves a hierarchical penalty which places an
  outer MCP penalty on a sum of inner MCP penalties for each group, as
  proposed in Breheny & Huang, 2009.  Either penalty may be referred to
  as the "group MCP", depending on the publication.  To resolve this
  confusion, Huang et al. (2012) proposed the name "composite MCP" for
  the \code{cMCP} penalty.

  For more information about the penalties and their properties, please
  consult the references below, many of which contain discussion, case
  studies, and simulation studies comparing the methods.  If you use
  \code{grpreg} for an analysis, please cite the appropriate reference.

  In keeping with the notation from the original MCP paper, the tuning
  parameter of the MCP penalty is denoted 'gamma'.  Note, however, that
  in Breheny and Huang (2009), \code{gamma} is denoted 'a'.

  The objective function for \code{grpreg} optimization is defined to be
  \deqn{Q(\beta|X, y) = \frac{1}{n} L(\beta|X, y) +
    P_\lambda(\beta)}{Q(\beta|X, y) = (1/n)*L(\beta|X, y) +
    P(\beta, \lambda),}
  where the loss function L is the deviance (-2 times the log
  likelihood) for the specified outcome distribution
  (gaussian/binomial/poisson).
  \href{http://pbreheny.github.io/ncvreg/articles/web/models.html}{See
    here for more details}.

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

  Although \code{grpreg} allows groups to be unordered and given
  arbitary names, it is recommended that you specify groups as
  consecutive integers.  The first reason is efficiency: if groups are
  out of order, \code{X} must be reordered prior to fitting, then this
  process reversed to return coefficients according to the original
  order of \code{X}.  This is inefficient if \code{X} is very large.
  The second reason is ambiguity with respect to other arguments such as
  \code{group.multiplier}.  With consecutive integers, \code{group=3}
  unambiguously denotes the third element of \code{group.multiplier}.

  Seemingly unrelated regressions/multitask learning can be carried out
  using \code{grpreg} by passing a matrix to \code{y}.  In this case,
  \code{X} will be used in separate regressions for each column of
  \code{y}, with the coefficients grouped across the responses.  In
  other words, each column of \code{X} will form a group with m
  members, where m is the number of columns of \code{y}.  For multiple
  Gaussian responses, it is recommended to standardize the columns of
  \code{y} prior to fitting, in order to apply the penalization equally
  across columns.

  \code{grpreg} requires groups to be non-overlapping.
}
\value{
  An object with S3 class \code{"grpreg"} containing:
  \describe{
    \item{beta}{The fitted matrix of coefficients.  The number of rows
      is equal to the number of coefficients, and the number of columns
      is equal to \code{nlambda}.}
    \item{family}{Same as above.}
    \item{group}{Same as above.}
    \item{lambda}{The sequence of \code{lambda} values in the path.}
    \item{alpha}{Same as above.}
    \item{loss}{A vector containing either the residual sum of squares
      (\code{"gaussian"}) or negative log-likelihood (\code{"binomial"})
      of the fitted model at each value of \code{lambda}.}
    \item{n}{Number of observations.}
    \item{penalty}{Same as above.}
    \item{df}{A vector of length \code{nlambda} containing estimates of
      effective number of model parameters all the points along the
      regularization path.  For details on how this is calculated, see
      Breheny and Huang (2009).}
    \item{iter}{A vector of length \code{nlambda} containing the number
      of iterations until convergence at each value of \code{lambda}.}
    \item{group.multiplier}{A named vector containing the multiplicative
      constant applied to each group's penalty.}
  }
}
\references{
  \itemize{
    \item Breheny, P. and Huang, J. (2009) Penalized methods for
    bi-level variable selection.  \emph{Statistics and its interface},
    \strong{2}: 369-380.
    \url{http://myweb.uiowa.edu/pbreheny/pdf/Breheny2009.pdf}

    \item Huang J., Breheny, P. and Ma, S. (2012). A selective
    review of group selection in high dimensional
    models. \emph{Statistical Science}, \strong{27}: 481-499.
    \url{http://myweb.uiowa.edu/pbreheny/pdf/Huang2012.pdf}

    \item Breheny, P. and Huang, J. (2015) Group descent algorithms for
    nonconvex penalized linear and logistic regression models with grouped
    predictors.  \emph{Statistics and Computing}, \strong{25}: 173-187.
    \url{http://dx.doi.org/10.1007/s11222-013-9424-2}

    \item Breheny, P. (2015) The group exponential lasso for bi-level
    variable selection. \emph{Biometrics}, \strong{71}: 731-740.
    \url{http://dx.doi.org/10.1111/biom.12300}
  }
}

\author{Patrick Breheny}
\seealso{\code{\link{cv.grpreg}}, as well as
  \code{\link[=plot.grpreg]{plot}} and
  \code{\link[=select.grpreg]{select}} methods.}
\examples{
# Birthweight data
data(Birthwt)
X <- Birthwt$X
group <- Birthwt$group

# Linear regression
y <- Birthwt$bwt
fit <- grpreg(X, y, group, penalty="grLasso")
plot(fit)
fit <- grpreg(X, y, group, penalty="grMCP")
plot(fit)
fit <- grpreg(X, y, group, penalty="grSCAD")
plot(fit)
fit <- grpreg(X, y, group, penalty="gel")
plot(fit)
fit <- grpreg(X, y, group, penalty="cMCP")
plot(fit)
select(fit, "AIC")

# Logistic regression
y <- Birthwt$low
fit <- grpreg(X, y, group, penalty="grLasso", family="binomial")
plot(fit)
fit <- grpreg(X, y, group, penalty="grMCP", family="binomial")
plot(fit)
fit <- grpreg(X, y, group, penalty="grSCAD", family="binomial")
plot(fit)
fit <- grpreg(X, y, group, penalty="gel", family="binomial")
plot(fit)
fit <- grpreg(X, y, group, penalty="cMCP", family="binomial")
plot(fit)
select(fit, "BIC")

# Multitask learning (simulated example)
set.seed(1)
n <- 50
p <- 10
k <- 5
X <- matrix(runif(n*p), n, p)
y <- matrix(rnorm(n*k, X[,1] + X[,2]), n, k)
fit <- grpreg(X, y)
# Note that group is set up automatically
fit$group
plot(fit)
}
