reorderGroups <- function(group, m, bilevel) {
  if (any(order(group) != 1:length(group)) | !is.numeric(group)) {
    reorder <- TRUE
    gf <- factor(group)
    if (any(levels(gf)=="0")) {
      gf <- relevel(gf, "0")
      g <- as.numeric(gf) - 1
    } else {
      g <- as.numeric(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    g <- group
    ord <- ord.inv <- NULL
  }
  J <- max(g)
  if (missing(m)) {
    m <- if (bilevel) rep(1, J) else sqrt(table(g[g!=0]))
  }
  if (length(m)!=max(g)) stop("Length of group.multiplier must equal number of penalized groups")
  if (reorder) {
    names(m) <- setdiff(levels(gf), "0")
  } else {
    names(m) <- paste0("G", unique(g[g!=0]))
  }
  if (storage.mode(m) != "double") storage.mode(m) <- "double"
  list(g=g, m=m, ord=ord, ord.inv=ord.inv, reorder=reorder)
}
