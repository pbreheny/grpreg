if (requireNamespace("tinytest", quietly=TRUE)) {
  if (length(unclass(packageVersion("grpreg"))[[1]]) == 4 | Sys.getenv('R_FORCE_TEST') == 'TRUE') {
    tinytest::test_package("grpreg", pattern="^[^_].*\\.[rR]$")
  }
}
