#' Plot coefficients from a "grpreg" object
#' 
#' Produces a plot of the coefficient paths for a fitted \code{grpreg} object.
#' 
#' 
#' @param x Fitted \code{"grpreg"} model.
#' @param alpha Controls alpha-blending.  Default is alpha=1.
#' @param legend.loc Where should the legend go?  If left unspecified, no
#' legend is drawn.  See \code{\link[graphics]{legend}} for details.
#' @param label If TRUE, annotates the plot with text labels in the right
#' margin describing which variable/group the corresponding line belongs to.
#' @param log.l Should horizontal axis be on the log scale?  Default is FALSE.
#' @param norm If \code{TRUE}, plot the norm of each group, rather than the
#' individual coefficients.
#' @param \dots Other graphical parameters to \code{plot}, \code{matlines}, or
#' \code{legend}
#' 
#' @seealso [grpreg()]
#' 
#' @examples
#' # Fit model to birthweight data
#' data(Birthwt)
#' X <- Birthwt$X
#' y <- Birthwt$bwt
#' group <- Birthwt$group
#' fit <- grpreg(X, y, group, penalty="grLasso")
#' 
#' # Plot (basic)
#' plot(fit)
#' 
#' # Plot group norms, with labels in right margin
#' plot(fit, norm=TRUE, label=TRUE)
#' 
#' # Plot (miscellaneous options)
#' myColors <- c("black", "red", "green", "blue", "yellow", "purple",
#' "orange", "brown")
#' plot(fit, legend.loc="topleft", col=myColors)
#' labs <- c("Mother's Age", "# Phys. visits", "Hypertension", "Mother's weight",
#'           "# Premature", "Race", "Smoking", "Uterine irritability")
#' plot(fit, legend.loc="topleft", lwd=6, alpha=0.5, legend=labs)
#' plot(fit, norm=TRUE, legend.loc="topleft", lwd=6, alpha=0.5, legend=labs)
#' @export

plot.grpreg <- function(x, alpha=1, legend.loc, label=FALSE, log.l=FALSE, norm=FALSE, ...) {
  if (norm) {
    Y <- predict(x, type="norm")
    if (any(x$group==0)) Y <- Y[-1,]
    nonzero <- which(apply(abs(Y), 1, sum)!=0)
    Y <- Y[nonzero,]
    g <- 1:nrow(Y)
  } else {
    if (length(dim(x$beta))==3) {
      beta <- matrix(x$beta[, -1, , drop=FALSE], ncol=dim(x$beta)[3])
    } else if (inherits(x, "grpsurv")) {
      beta <- x$beta
    } else {
      beta <- x$beta[-1, , drop=FALSE]
    }
    penalized <- which(x$group!=0)
    nonzero <- which(apply(abs(beta), 1, sum)!=0)
    ind <- intersect(penalized, nonzero)
    Y <- beta[ind, , drop=FALSE]
    g <- as.integer(as.factor(x$group[ind]))
  }
  p <- nrow(Y)
  l <- x$lambda
  n.g <- max(g)

  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)

  plot.args <- list(x=l, y=1:length(l), ylim=range(Y), xlab=xlab, ylab="", type="n", xlim=rev(range(l)), las=1, bty="n")
  new.args <- list(...)
  if (length(new.args)) {
    new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
    plot.args[names(new.plot.args)] <- new.plot.args
  }
  do.call("plot", plot.args)
  if (plot.args$ylab=="") {
    ylab <- if (norm) expression("||"*hat(beta)*"||") else expression(hat(beta))
    mtext(ylab, 2, 3.5, las=1, adj=0)
  }
  abline(h=0, lwd=0.5, col="gray")

  cols <- hcl(h=seq(15, 375, len=max(4, n.g+1)), l=60, c=150, alpha=alpha)
  cols <- if (n.g==2) cols[c(1,3)] else cols[1:n.g]
  line.args <- list(col=cols, lwd=1+2*exp(-p/20), lty=1, pch="")
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- l
  line.args$y <- t(Y)
  line.args$col <- line.args$col[g]
  line.args$lty <- rep(line.args$lty, length.out=max(g))
  line.args$lty <- line.args$lty[g]
  do.call("matlines", line.args)

  if(!missing(legend.loc)) {
    legend.args <- list(col=cols, lwd=line.args$lwd, lty=line.args$lty, legend=names(x$group.multiplier))
    if (length(new.args)) {
      new.legend.args <- new.args[names(new.args) %in% names(formals(legend))]
      legend.args[names(new.legend.args)] <- new.legend.args
    }
    legend.args$x <- legend.loc
    do.call("legend", legend.args)
  }
  if (label) {
    ypos <- Y[, ncol(Y)]
    text(-0.001, ypos, names(ypos), xpd=NA, adj=c(0, NA))
  }
}
