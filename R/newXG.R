newXG <- function(X, g, m, ncolY, bilevel) {
  # Coerce X to matrix
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (any(is.na(X))) stop("Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg")

  # Reorder groups, if necessary
  if (length(g) != ncol(X)) stop ("Dimensions of group is not compatible with X")
  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  grp <- reorderGroups(g, m, bilevel)
  g <- grp$g
  m <- grp$m
  if (grp$reorder) X <- X[,grp$ord]

  # Make multiX, if necessary
  if (ncolY > 1) {
    X <- multiX(X, ncolY)
    g <- c(rep(0, ncolY-1), rep(g, each=ncolY))
  }

  # Standardize
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)                        # non-constant columns
  zg <- setdiff(unique(g), unique(g[nz]))  # constant groups
  if (length(zg)) m <- m[-zg]
  if (length(nz) != ncol(X)) {
    XX <- XX[ ,nz, drop=FALSE]
    g <- g[nz]
  }
  if (!bilevel) {
    XX <- orthogonalize(XX, g)
    g <- attr(XX, "group")
  }

  # Return
  return(list(X=XX, g=g, m=m, reorder=grp$reorder, ord.inv=grp$ord.inv, names=xnames,
              center=center[nz], scale=scale[nz], nz=nz))
}
