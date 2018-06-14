newXG <- function(X, g, m, ncolY, bilevel) {
  # Coerce X to matrix
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (any(is.na(X))) stop("Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg")
  if (length(g) != ncol(X)) stop ("Dimensions of group is not compatible with X")
  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)

  # Setup group
  G <- setupG(g, m, bilevel)

  # Reconfigure for multiple outcomes, if necessary
  if (ncolY > 1) {
    X <- multiX(X, ncolY)
    G <- multiG(G, ncolY)
  }

  # Feature-level standardization
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)                # non-constant columns
  if (length(nz) != ncol(X)) {
    XX <- XX[ ,nz, drop=FALSE]
    G <- subsetG(G, nz)
  }

  # Reorder groups, if necessary
  G <- reorderG(G, attr(G, 'm'), bilevel)
  if (attr(G, 'reorder')) XX <- XX[,attr(G, 'ord')]

  # Group-level standardization
  if (!bilevel) {
    XX <- orthogonalize(XX, G)
    g <- attr(XX, "group")
  } else {
    g <- as.integer(G)
  }

  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- if (bilevel) rep(1, max(g)) else sqrt(table(g[g!=0]))
  }

  # Return
  return(list(X=XX, g=g, m=m, reorder=attr(G, 'reorder'), ord.inv=attr(G, 'ord.inv'), names=xnames,
              center=center, scale=scale, nz=nz))
}
