[![GitHub version](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/pbreheny/grpreg/master/.version.json&style=flat&logo=github)](https://github.com/pbreheny/grpreg)
[![CRAN version](https://img.shields.io/cran/v/grpreg?logo=R)](https://cran.r-project.org/package=grpreg)
[![downloads](https://cranlogs.r-pkg.org/badges/grpreg)](https://cran.r-project.org/package=grpreg)
[![Travis build status](https://api.travis-ci.org/pbreheny/grpreg.svg?branch=master)](https://travis-ci.org/pbreheny/grpreg)
[![codecov.io](https://codecov.io/github/pbreheny/grpreg/coverage.svg?branch=master)](https://codecov.io/github/pbreheny/grpreg?branch=master)

# [Regularization Paths for Regression Models with Grouped Covariates](https://pbreheny.github.io/grpreg/)

`grpreg` is an R package for fitting the regularization path of linear regression, GLM, and Cox regression models with grouped penalties.  This includes group selection methods such as group lasso, group MCP, and group SCAD as well as bi-level selection methods such as the group exponential lasso, the composite MCP, and the group bridge.  Utilities for carrying out cross-validation as well as post-fitting visualization, summarization, and prediction are also provided.

### Install

* To install the latest release version from CRAN: `install.packages("grpreg")`
* To install the latest development version from GitHub: `remotes::install_github("pbreheny/grpreg")`

### Get started

See the ["getting started" vignette](https://pbreheny.github.io/grpreg/articles/getting-started.html)

### Learn more

Follow the links under "Learn more" at the [grpreg website](https://pbreheny.github.io/grpreg/)

### Details of the algorithms used

* [Breheny, P. and Huang, J. (2009) Penalized methods for bi-level variable selection.  *Statistics and its interface*, **2**: 369-380.](https://myweb.uiowa.edu/pbreheny/pdf/Breheny2009.pdf)
* [Breheny, P. and Huang, J. (2015) Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors. *Statistics and Computing*, **25**: 173-187.](https://dx.doi.org/10.1007/s11222-013-9424-2)
