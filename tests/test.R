if (requireNamespace("tinytest", quietly=TRUE)) {
  if (length(unclass(packageVersion("grpreg"))[[1]]) == 4) {
    tinytest::test_package("grpreg", pattern="^[^_].*\\.[rR]$")
  }
}
