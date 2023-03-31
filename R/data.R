#' Risk Factors Associated with Low Infant Birth Weight
#' 
#' This version of the data set has been deprecated and will not be supported
#' in future versions.  Please use \code{\link{Birthwt}} instead.
#' 
#' @format This data frame contains the following columns:
#' \itemize{
#' \item\code{low} Indicator of birth weight less than 2.5kg \item\code{bwt}
#' Birth weight in kilograms \item\code{age1,age2,age3} Orthogonal polynomials
#' of first, second, and third degree representing mother's age in years
#' \item\code{lwt1,lwt2,lwt3} Orthogonal polynomials of first, second, and
#' third degree representing mother's weight in pounds at last menstrual period
#' \item\code{white,black} Indicator functions for mother's race; "other" is
#' reference group \item\code{smoke} smoking status during pregnancy
#' \item\code{ptl1,ptl2m} Indicator functions for one or for two or more
#' previous premature labors, respectively.  No previous premature labors is
#' the reference category.  \item\code{ht} History of hypertension
#' \item\code{ui} Presence of uterine irritability \item\code{ftv1,ftv2,ftv3m}
#' Indicator functions for one, for two, or for three or more physician visits
#' during the first trimester, respectively.  No visits is the reference
#' category.
#' }
#' @seealso \code{\link{Birthwt}}
"birthwt.grpreg"


#' Risk Factors Associated with Low Infant Birth Weight
#' 
#' The `Birthwt` data contains 189 observations, 16 predictors, and an
#' outcome, birthweight, available both as a continuous measure and a binary
#' indicator for low birth weight.The data were collected at Baystate Medical
#' Center, Springfield, Mass during 1986. This data frame is a
#' reparameterization of the `birthwt` data frame from the **MASS** package.
#' 
#' @format The \code{Birthwt} object is a list containing four elements (`X`, `bwt`, `low`, and `group`):
#' \describe{
#'   \item{bwt}{Birth weight in kilograms}
#'   \item{low}{Indicator of birth weight less than 2.5kg}
#'   \item{group}{Vector describing how the columns of X are grouped}
#'   \item{X}{A matrix with 189 observations (rows) and 16 predictor variables (columns).}
#' }
#' The matrix `X` contains the following columns:
#' \describe{
#'   \item{age1,age2,age3}{Orthogonal polynomials of first, second, and third degree representing mother's age in years}
#'   \item{lwt1,lwt2,lwt3}{Orthogonal polynomials of first, second, and third degree representing mother's weight in pounds at last menstrual period}
#'   \item{white,black}{Indicator functions for mother's race; "other" is reference group}
#'   \item{smoke}{Smoking status during pregnancy}
#'   \item{ptl1,ptl2m}{Indicator functions for one or for two or more previous premature labors, respectively. No previous premature labors is the reference category.}
#'   \item{ht}{History of hypertension}
#'   \item{ui}{Presence of uterine irritability}
#'   \item{ftv1,ftv2,ftv3m}{Indicator functions for one, for two, or for three or more physician visits during the first trimester, respectively. No visits is the reference category.}
#' }
#' 
#' @seealso [MASS::birthwt], [grpreg()]
#' 
#' @references
#' \itemize{
#'   \item Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics with S.* Fourth edition. Springer.
#'   \item Hosmer, D.W. and Lemeshow, S. (1989) *Applied Logistic Regression.* New York: Wiley
#' }
#' 
#' @source <https://cran.r-project.org/package=MASS>
#' 
#' @examples
#' data(Birthwt)
#' hist(Birthwt$bwt, xlab="Child's birth weight", main="")
#' table(Birthwt$low)
#' ## See examples in ?birthwt (MASS package)
#' ##   for more about the data set
#' ## See examples in ?grpreg for use of this data set
#' ##   with group penalized regression models
"Birthwt"


#' VA lung cancer data set
#' 
#' Data from a randomised trial of two treatment regimens for lung cancer. This
#' is a standard survival analysis data set from the classic textbook by
#' Kalbfleisch and Prentice.
#' 
#' @format A list of two objects: `y` and `X`
#' \describe{
#'   \item{y}{A two column matrix (`Surv` object) containing the follow-up
#'   time (in days) and an indicator variable for whether the patient died
#'   while on the study or not.}
#'   \item{X}{A matrix with 137 observations (rows) and 9 predictor variables
#'   (columns). The remainder of this list describes the columns of `X`}
#'   \item{trt}{Treatment indicator (1=control group, 2=treatment group)}
#'   \item{karno}{Karnofsky performance score (0=bad, 100=good)}
#'   \item{diagtime}{Time from diagnosis to randomization (months)}
#'   \item{age}{Age (years, at baseline)}
#'   \item{prior}{Prior therapy (0=no, 1=yes)}
#'   \item{squamous}{Indicator for whether the cancer type is squamous cell
#'   carcinoma (0=no, 1=yes)}
#'   \item{small}{Indicator for whether the cancer type is small cell lung
#'   cancer (0=no, 1=yes)}
#'   \item{adeno}{Indicator for whether the cancer type is adenocarcinoma
#'   (0=no, 1=yes)}
#'   \item{large}{Indicator for whether the cancer type is large cell carcinoma
#'   (0=no, 1=yes)}
#' }
#' 
#' @seealso `grpsurv()`
#' @references
#' \itemize{
#'   \item Kalbfleisch D and Prentice RL (1980), *The Statistical Analysis of
#'   Failure Time Data*. Wiley, New York.
#' }
#' @source \url{https://cran.r-project.org/package=survival}
"Lung"
