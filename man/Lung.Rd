\name{Lung}
\alias{Lung}
\title{VA lung cancer data set}
\description{Data from a randomised trial of two treatment regimens for
  lung cancer. This is a standard survival analysis data set from the 
  classic textbook by Kalbfleisch and Prentice.
}
\usage{
data(Lung)
}
\format{
  The \code{Lung} object is a list containing three elements:
  \itemize{
    \item{\code{y}:}{A two column matrix (\code{\link[survival]{Surv}}
    object) containing the follow-up time (in days) and an indicator
    variable for whether the patient died while on the study or not.}
    \item{\code{X}:}{Matrix of predictors}
    \item{\code{group}:}{Vector describing how the columns of X are
    grouped}
  }
  The matrix \code{X} contains the following columns:
  \itemize{
    \item{\code{trt}:}{Treatment indicator (1=control group, 2=treatment
    group)}
    \item{\code{karno1}, \code{karno2}, \code{karno3}:}{Orthogonal
    polynomials of first, second, and third degree representing
    Karnofsky performance score (0=bad, 100=good)}
    \item{\code{diagtime1}, \code{diagtime2}:}{Orthogonal
    polynomials of first and second degree representing time from
    diagnosis to randomization (months)}
    \item{\code{age1}, \code{age2}, \code{age3}:}{Orthogonal
    polynomials of first, second, and third degree representing the
    patient's age in years}
    \item{\code{prior}:}{Prior therapy (0=no, 1=yes)}
    \item{\code{squamous}, \code{small}, \code{adeno},
      \code{large}:}{Indicators for the lung cancer cell type.  For
    each subject, exactly one of these columns will be 1 and the rest
    0.}
  }
}
\source{\url{https://cran.r-project.org/package=survival}}
\references{
  \itemize{
    \item Kalbfleisch D and Prentice RL (1980), \emph{The Statistical
  Analysis of Failure Time Data}. Wiley, New York.
  }
}
\seealso{\code{\link[survival]{veteran}}, \code{\link{grpsurv}}}
\examples{
data(Lung)
hist(Lung$y[,1], xlab="Follow-up time", main="")
table(Lung$y[,2])
}
