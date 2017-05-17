newY <- function(y, family) {
  if (is.data.frame(y)) y <- as.matrix(y)
  if (is.matrix(y)) {
    d <- dim(y)
    y <- t(y)
  } else {
    d <- c(length(y), 1)
  }

  # Convert fuzzy binomial data
  if (family=="binomial" && typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2) stop("Attemping to use family='binomial' with non-binary data", call.=FALSE)
    if (!identical(names(tab), c("0", "1"))) {
      print(paste0("Logistic regression modeling Pr(y=", names(tab)[2], ")"))
      y <- as.numeric(as.character(y) == names(tab)[2])
      if (d[2] > 1) attr(y, "dim") <- d
    }
  }

  # Convert to double, if necessary
  if (typeof(y) != "double") {
    tryCatch(storage.mode(y) <- "double", warning=function(w) {stop("y must be numeric or able to be coerced to numeric", call.=FALSE)})
  }
  if (any(is.na(y))) stop("Missing data (NA's) detected in outcome y.  You must eliminate missing data (e.g., by removing cases or imputation) before passing y to grpreg")

  # Handle multi
  if (is.matrix(y)) {
    if (ncol(y) > 1) {
      if (is.null(colnames(y))) paste("Y", 1:ncol(y), sep="")
    }
    attributes(y) <- NULL
  }

  if (family=="gaussian") {
    meanY <- mean(y)
    y <- y - meanY
    attr(y, "mean") <- meanY
  }
  attr(y, "m") <- d[2]
  y
}
