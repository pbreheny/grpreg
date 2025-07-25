<!-- badges: start -->
[![GitHub version](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/pbreheny/grpreg/master/.version.json&style=flat&logo=github)](https://github.com/pbreheny/grpreg)
[![CRAN version](https://img.shields.io/cran/v/grpreg?logo=R)](https://cran.r-project.org/package=grpreg)
[![downloads](https://cranlogs.r-pkg.org/badges/grpreg)](https://cran.r-project.org/package=grpreg)
[![R-CMD-check](https://github.com/pbreheny/grpreg/workflows/R-CMD-check/badge.svg)](https://github.com/pbreheny/grpreg/actions)
[![codecov.io](https://codecov.io/github/pbreheny/grpreg/coverage.svg?branch=master)](https://app.codecov.io/gh/pbreheny/grpreg)
<!-- badges: end -->

# [Regularization Paths for Regression Models with Grouped Covariates](https://pbreheny.github.io/grpreg/)

`grpreg` is an R package for fitting the regularization path of linear regression, GLM, and Cox regression models with grouped penalties.  This includes group selection methods such as group lasso, group MCP, and group SCAD as well as bi-level selection methods such as the group exponential lasso, the composite MCP, and the group bridge.  Utilities for carrying out cross-validation as well as post-fitting visualization, summarization, and prediction are also provided.

### Install

* To install the latest release version from CRAN:

``` r
install.packages("grpreg")
```

* To install the latest development version from GitHub:

``` r
remotes::install_github("pbreheny/grpreg")
```

### Get started

See the ["getting started" vignette](https://pbreheny.github.io/grpreg/articles/grpreg.html)

### Learn more

Follow the links under "Articles" at the [grpreg website](https://pbreheny.github.io/grpreg/)

### References

For more on the mathematical foundations and algorithmic details, see:

* [Breheny, P. and Huang, J. (2009) Penalized methods for bi-level variable selection.  *Statistics and its interface*, **2**: 369-380.](https://myweb.uiowa.edu/pbreheny/pdf/Breheny2009.pdf)
* [Breheny, P. and Huang, J. (2015) Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors. *Statistics and Computing*, **25**: 173-187.](https://dx.doi.org/10.1007/s11222-013-9424-2)
