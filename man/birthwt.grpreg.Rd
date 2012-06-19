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
\describe{
\item{low}{indicator of birth weight less than 2.5kg}
\item{bwt}{birth weight in kilograms}
\item{age1,age2,age3}{Orthogonal polynomials of
  first, second, and third degree representing mother's age in years}
\item{lwt1,lwt2,lwt3}{Orthogonal polynomials of
  first, second, and third degree representing mother's weight in pounds
  at last menstrual period}
\item{white,black}{Indicator functions for mother's race;
  "other" is reference group}
\item{smoke}{smoking status during pregnancy}
\item{ptl1,ptl2m}{Indicator functions for one or for two or
  more previous premature labors, respectively.  No previous premature
  labors is the reference category.}
\item{ht}{history of hypertension}
\item{ui}{presence of uterine irritability}
\item{ftv1,ftv2,ftv3m}{Indicator functions for
  one, for two, or for three or more physician visits during the first
  trimester, respectively.  No visits is the reference category.}}}
\source{\code{MASS}.  \code{R} package.  Available at \code{http://cran.r-project.org}.}
\references{
  Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied
    Statistics with S.} Fourth edition. Springer. \cr
  Hosmer, D.W. and Lemeshow, S. (1989) \emph{Applied Logistic
    Regression.} New York: Wiley}
\keyword{datasets}
