![version](http://www.r-pkg.org/badges/version/grpreg)
![downloads](http://cranlogs.r-pkg.org/badges/grpreg)
[![codecov.io](https://codecov.io/github/pbreheny/grpreg/coverage.svg?branch=master)](https://codecov.io/github/pbreheny/grpreg?branch=master)

# Regularization Paths for Regression Models with Grouped Covariates

`grpreg` is an R package for fitting the regularization path of linear regression, GLM, and Cox regression models with grouped penalties.  This includes group selection methods such as group lasso, group MCP, and group SCAD as well as bi-level selection methods such as the group exponential lasso, the composite MCP, and the group bridge.  Utilities for carrying out cross-validation as well as post-fitting visualization, summarization, and prediction are also provided.

## Installation

* To install the latest release version from CRAN: `install.packages("grpreg")`
* To install the latest development version from GitHub: `devtools::install_github("pbreheny/grpreg")`
