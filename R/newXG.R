newXG <- function(X, g, m, ncolY) {
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (any(is.na(X))) stop("Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg")

  if (length(g) != ncol(X)) stop ("Dimensions of group is not compatible with X")
  grp <- reorderGroups(group, group.multiplier, strtrim(penalty,2)=="gr")

  if (grp$reorder) X <- X[,grp$ord]

  xnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)

  if (missing(m)) {
    m <- if (sqr) sqrt(table(g[g!=0])) else rep(1, J)
  }
  if (length(m)!=max(g)) stop("Length of group.multiplier must equal number of penalized groups")
  if ()

  obj <- XG$new(X=X, g=g, m=m)
  ## Error checking

  # Reorder groups, if necessary

  if (ncolY > 1) {
    X <- multiX(X, ncolY)
    group <- grp$g <- c(rep(0, m-1), rep(grp$g, each=m))
    grp$m <- rep(grp$m, each=m)
  }
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)
  zg <- setdiff(unique(grp$g), unique(grp$g[nz]))
  if (length(zg)) grp$m <- grp$m[-zg]
  XX <- XX[ ,nz, drop=FALSE]
  grp$g <- grp$g[nz]
  if (strtrim(penalty,2)=="gr") {
    XX <- orthogonalize(XX, grp$g)
    grp$g <- attr(XX, "group")
  }
  K <- as.numeric(table(grp$g))


}
