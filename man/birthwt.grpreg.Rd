\name{birthwt.grpreg}
\alias{birthwt.grpreg}
\title{
Risk Factors Associated with Low Infant Birth Weight
}
\description{
The \code{birthwt.grpreg} data frame has 189 rows and 18 columns.
The data were collected at Baystate Medical Center, Springfield, Mass
during 1986.  This data frame is a reparameterization of the
\code{birthwt} data frame from the \code{MASS} package.}
\usage{
data(birthwt.grpreg)
}
\format{
  This data frame contains the following columns:
  \itemize{
    \item{\code{low}} Indicator of birth weight less than 2.5kg
    \item{\code{bwt}} Birth weight in kilograms
    \item{\code{age1,age2,age3}} Orthogonal polynomials of
      first, second, and third degree representing mother's age in years
    \item{\code{lwt1,lwt2,lwt3}} Orthogonal polynomials of
      first, second, and third degree representing mother's weight in pounds
      at last menstrual period
    \item{\code{white,black}} Indicator functions for mother's race;
      "other" is reference group
    \item{\code{smoke}} smoking status during pregnancy
    \item{\code{ptl1,ptl2m}} Indicator functions for one or for two or
      more previous premature labors, respectively.  No previous premature
      labors is the reference category.
    \item{\code{ht}} History of hypertension
    \item{\code{ui}} Presence of uterine irritability
    \item{\code{ftv1,ftv2,ftv3m}} Indicator functions for
      one, for two, or for three or more physician visits during the first
      trimester, respectively.  No visits is the reference category.
    }
  }
  
\source{\code{MASS}.  \code{R} package.  Available at \code{http://cran.r-project.org}.}
\references{
  \itemize{
    \item Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied
      Statistics with S.} Fourth edition. Springer.
    \item Hosmer, D.W. and Lemeshow, S. (1989) \emph{Applied Logistic
      Regression.} New York: Wiley}
}
\examples{
data(birthwt.grpreg)
hist(birthwt.grpreg$bwt, xlab="Child's birth weight", main="")
table(birthwt.grpreg$low)
## See examples in ?birthwt (MASS package)
##   for more about the data set
## See examples in grpreg for use of this data set
##   with group penalized regression models
}
\keyword{datasets}
