#' Plots the cross-validation curve from a \code{cv.grpreg} object
#' 
#' Plots the cross-validation curve from a \code{cv.grpreg} object, along with
#' standard error bars.
#' 
#' Error bars representing approximate +/- 1 SE (68\% confidence intervals) are
#' plotted along with the estimates at value of \code{lambda}.  For \code{rsq}
#' and \code{snr}, these confidence intervals are quite crude, especially near
#' zero, and will hopefully be improved upon in later versions of
#' \code{grpreg}.
#' 
#' @param x A \code{cv.grpreg} object.
#' @param log.l Should horizontal axis be on the log scale?  Default is TRUE.
#' @param type What to plot on the vertical axis.  \code{cve} plots the
#' cross-validation error (deviance); \code{rsq} plots an estimate of the
#' fraction of the deviance explained by the model (R-squared); \code{snr}
#' plots an estimate of the signal-to-noise ratio; \code{scale} plots, for
#' \code{family="gaussian"}, an estimate of the scale parameter (standard
#' deviation); \code{pred} plots, for \code{family="binomial"}, the estimated
#' prediction error; \code{all} produces all of the above.
#' @param selected If \code{TRUE} (the default), places an axis on top of the
#' plot denoting the number of groups in the model (i.e., that contain a
#' nonzero regression coefficient) at that value of \code{lambda}.
#' @param vertical.line If \code{TRUE} (the default), draws a vertical line at
#' the value where cross-validaton error is minimized.
#' @param col Controls the color of the dots (CV estimates).
#' @param \dots Other graphical parameters to \code{plot}
#' 
#' @seealso [grpreg()], [cv.grpreg()]
#' 
#' @examples
#' \dontshow{set.seed(1)}
#' # Birthweight data
#' data(Birthwt)
#' X <- Birthwt$X
#' group <- Birthwt$group
#' 
#' # Linear regression
#' y <- Birthwt$bwt
#' cvfit <- cv.grpreg(X, y, group)
#' plot(cvfit)
#' op <- par(mfrow=c(2,2))
#' plot(cvfit, type="all")
#' 
#' ## Logistic regression
#' y <- Birthwt$low
#' cvfit <- cv.grpreg(X, y, group, family="binomial")
#' par(op)
#' plot(cvfit)
#' par(mfrow=c(2,2))
#' plot(cvfit, type="all")
#' @export

plot.cv.grpreg <- function(x, log.l=TRUE, type=c("cve", "rsq", "scale", "snr", "pred", "all"), selected=TRUE, vertical.line=TRUE, col="red", ...) {
  type <- match.arg(type)
  if (type=="all") {
    plot(x, log.l=log.l, type="cve", selected=selected, ...)
    plot(x, log.l=log.l, type="rsq", selected=selected, ...)
    plot(x, log.l=log.l, type="snr", selected=selected, ...)
    if (length(x$fit$family)) {
      if (x$fit$family == "binomial") plot(x, log.l=log.l, type="pred", selected=selected, ...)
      if (x$fit$family == "gaussian") plot(x, log.l=log.l, type="scale", selected=selected, ...)
    }
    return(invisible(NULL))
  }
  l <- x$lambda
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)

  ## Calculate y
  L.cve <- x$cve - x$cvse
  U.cve <- x$cve + x$cvse
  if (type=="cve") {
    y <- x$cve
    L <- L.cve
    U <- U.cve
    ylab <- "Cross-validation error"
  } else if (type=="rsq" | type == "snr") {
    if (length(x$fit$family) && x$fit$family=='gaussian') {
      rsq <- pmin(pmax(1 - x$cve/x$null.dev, 0), 1)
      rsql <- pmin(pmax(1 - U.cve/x$null.dev, 0), 1)
      rsqu <- pmin(pmax(1 - L.cve/x$null.dev, 0), 1)
    } else {
      rsq <- pmin(pmax(1 - exp(x$cve-x$null.dev), 0), 1)
      rsql <- pmin(pmax(1 - exp(U.cve-x$null.dev), 0), 1)
      rsqu <- pmin(pmax(1 - exp(L.cve-x$null.dev), 0), 1)
    }
    if (type == "rsq") {
      y <- rsq
      L <- rsql
      U <- rsqu
      ylab <- ~R^2
    } else if(type=="snr") {
      y <- pmin(rsq/(1-rsq), 1e6)
      L <- pmin(rsql/(1-rsql), 1e6)
      U <- pmin(rsqu/(1-rsqu), 1e6)
      if (max(c(y,L,U)) == 1e6) warning('Signal-to-noise ratio is infinite')
      ylab <- "Signal-to-noise ratio"
    }
  } else if (type=="scale") {
    if (x$fit$family == "binomial") stop("Scale parameter for binomial family fixed at 1", call.=FALSE)
    y <- sqrt(x$cve)
    L <- sqrt(L.cve)
    U <- sqrt(U.cve)
    ylab <- ~hat(sigma)
  } else if (type=="pred") {
    y <- x$pe
    n <- x$fit$n
    CI <- sapply(y, function(x) {binom.test(x*n, n, conf.level=0.68)$conf.int})
    L <- CI[1,]
    U <- CI[2,]
    ylab <- "Prediction error"
  }

  ind <- if (type=="pred") is.finite(l[1:length(x$pe)]) else is.finite(l[1:length(x$cve)])
  ylim <- range(c(L[ind], U[ind]))
  aind <- ((U-L)/diff(ylim) > 1e-3) & ind
  plot.args = list(x=l[ind], y=y[ind], ylim=ylim, xlab=xlab, ylab=ylab, type="n", xlim=rev(range(l[ind])), las=1, bty="n")
  new.args = list(...)
  if (length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  if (vertical.line) abline(v=l[x$min], lty=2, lwd=.5)
  suppressWarnings(arrows(x0=l[aind], x1=l[aind], y0=L[aind], y1=U[aind], code=3, angle=90, col="gray80", length=.05))
  points(l[ind], y[ind], col=col, pch=19, cex=.5)
  if (selected) {
    n.s <- sapply(predict(x$fit, lambda=x$lambda, type="groups"), length)
    axis(3, at=l, labels=n.s, tick=FALSE, line=-0.5)
    mtext("Groups selected", cex=0.8, line=1.5)
  }
}
