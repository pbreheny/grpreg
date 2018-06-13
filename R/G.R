setupG <- function(group, m, bilevel) {
  gf <- factor(group)
  if (any(levels(gf)=='0')) {
    g <- as.integer(gf) - 1
    lev <- levels(gf)[levels(gf)!='0']
  } else {
    g <- as.integer(gf)
    lev <- levels(gf)
  }
  if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  if (missing(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    #if (all.equal(sort(names(m)), sort(group)))
    TRY <- try(as.integer(group)==g)
    if (class(TRY)=='try-error' || !TRY) stop('Attempting to set group.multiplier is ambiguous if group is not a factor')
    if (length(m) != length(lev)) stop("Length of group.multiplier must equal number of penalized groups")
    if (storage.mode(m) != "double") storage.mode(m) <- "double"
  }
  structure(g, levels=lev, m=m)
}
subsetG <- function(g, nz) {
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  new <- g[nz]
  dropped <- setdiff(g, new)
  if (length(dropped)) {
    lev <- lev[-dropped]
    m <- m[-dropped]
    gf <- factor(new)
    new <- as.integer(gf) - 1*any(levels(gf)=='0')
  }
  structure(new, levels=lev, m=m)
}
reorderG <- function(g, m, bilevel) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(g==0)) {
    g <- as.integer(relevel(factor(g), "0"))-1
  }
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    if (any(levels(gf)=="0")) {
      gf <- relevel(gf, "0")
      g <- as.integer(gf) - 1
    } else {
      g <- as.integer(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(g, levels=lev, m=m, ord=ord, ord.inv=ord.inv, reorder=reorder)
}
