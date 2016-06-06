\name{birthwt.grpreg}
\alias{birthwt.grpreg}
\title{
Risk Factors Associated with Low Infant Birth Weight
}
\description{
  This version of the data set has been deprecated and will not be
  supported in future versions.  Please use \code{\link{Birthwt}}
  instead.}
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
\seealso{\code{\link{Birthwt}}}  
