#' @keywords internal 
#' @aliases grpreg-package NULL
#' @author Patrick Breheny
#' 
#' @references
#' \itemize{
#' \item Yuan M and Lin Y. (2006) Model selection and estimation in regression
#' with grouped variables. *Journal of the Royal Statistical Society Series B*,
#' **68**: 49-67. \doi{10.1111/j.1467-9868.2005.00532.x}
#' 
#' \item Huang J, Ma S, Xie H, and Zhang C. (2009) A group bridge approach for
#' variable selection. *Biometrika*, **96**: 339-355. \doi{10.1093/biomet/asp020}
#' 
#' \item Breheny P and Huang J. (2009) Penalized methods for bi-level variable
#' selection. *Statistics and its interface*, **2**: 369-380.
#' \doi{10.4310/sii.2009.v2.n3.a10}
#' 
#' \item Huang J, Breheny P, and Ma S. (2012). A selective review of group
#' selection in high dimensional models. *Statistical Science*, **27**: 481-499.
#' \doi{10.1214/12-sts392}
#' 
#' \item Breheny P and Huang J. (2015) Group descent algorithms for nonconvex
#' penalized linear and logistic regression models with grouped predictors.
#' *Statistics and Computing*, **25**: 173-187. \doi{10.1007/s11222-013-9424-2}
#' 
#' \item Breheny P. (2015) The group exponential lasso for bi-level variable
#' selection. *Biometrics*, **71**: 731-740. \doi{10.1111/biom.12300}
#' }
#' 
#' @examples
#' \donttest{vignette("getting-started", package="grpreg")}
"_PACKAGE"

#' @useDynLib grpreg, .registration = TRUE
#' @import stats
#' @import graphics
#' @import grDevices
#' @import utils
#' @importFrom Matrix bdiag
NULL
