[![version](http://www.r-pkg.org/badges/version/grpreg)](https://cran.r-project.org/package=grpreg)
[![downloads](http://cranlogs.r-pkg.org/badges/grpreg)](https://cran.r-project.org/package=grpreg)
[![codecov.io](https://codecov.io/github/pbreheny/grpreg/coverage.svg?branch=master)](https://codecov.io/github/pbreheny/grpreg?branch=master)
[![Travis build
status](https://travis-ci.org/pbreheny/grpreg.svg?branch=master)](https://travis-ci.org/pbreheny/grpreg)

# [Regularization Paths for Regression Models with Grouped Covariates](http://pbreheny.github.io/grpreg)

`grpreg` is an R package for fitting the regularization path of linear regression, GLM, and Cox regression models with grouped penalties.  This includes group selection methods such as group lasso, group MCP, and group SCAD as well as bi-level selection methods such as the group exponential lasso, the composite MCP, and the group bridge.  Utilities for carrying out cross-validation as well as post-fitting visualization, summarization, and prediction are also provided.

### Install

* To install the latest release version from CRAN: `install.packages("grpreg")`
* To install the latest development version from GitHub: `devtools::install_github("pbreheny/grpreg")`

### Get started

See the ["getting started" vignette](http://pbreheny.github.io/grpreg/articles/getting-started.html)

### Learn more

Follow the links under "Learn more" at the [grpreg website](http://pbreheny.github.io/grpreg)

### Details of the algorithms used

* [Breheny, P. and Huang, J. (2009) Penalized methods for bi-level variable selection.  *Statistics and its interface*, **2**: 369-380.](http://myweb.uiowa.edu/pbreheny/pdf/Breheny2009.pdf)
* [Breheny, P. and Huang, J. (2015) Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors. *Statistics and Computing*, **25**: 173-187.](http://dx.doi.org/10.1007/s11222-013-9424-2)
