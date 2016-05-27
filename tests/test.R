require(grpreg)
test <- function(name) {
  filename <- system.file(paste0("tests/", name, ".R"), package = "grpreg")
  source(filename)
}
check <- function(x, y, ...) {
  if (missing(y)) {
    xname <- gsub("()", "", match.call()[2])
    if (x==TRUE) return(TRUE)
    message <- paste0("Problem in ", .test, "\n", xname, " FALSE")
  }
  checkResult <- all.equal(x, y, ...)
  if (class(checkResult)[1]=="logical") return(TRUE)
  xname <- gsub("()", "", match.call()[2])
  yname <- gsub("()", "", match.call()[3])
  message <- paste0("Problem in ", .test, "\n", xname, " not equal to ", yname, "\n", checkResult)
  stop(message, call.=FALSE)
}

test("standardization-orthogonalization")
