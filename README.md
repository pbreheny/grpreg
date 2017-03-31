[![Build Status](https://travis-ci.org/YaohuiZeng/grpreg.svg?branch=master)](https://travis-ci.org/YaohuiZeng/grpreg)
![version](http://www.r-pkg.org/badges/version/grpreg)
![downloads](http://cranlogs.r-pkg.org/badges/grpreg)
[![codecov.io](https://codecov.io/github/pbreheny/grpreg/coverage.svg?branch=master)](https://codecov.io/github/pbreheny/grpreg?branch=master)

`grpreg` fits regularization paths for linear, logistic, or Poisson regression models with grouped penalties, such as the group lasso, group MCP, group SCAD, group exponential lasso, and group bridge. The algorithms are based on the idea of either locally approximated coordinate descent or group descent, depending on the penalty. All of the algorithms (with the exception of group bridge) are stable and fast.

## News:
New efficient group-wise feature screening rules for group lasso are implemented into the package since Version 3.1-1:
* SSR: sequential Strong Rule;
* SEDPP: Sequential Enhanced Dual Polytope Projection;
* SSR-BEDPP: the hybrid safe strong rule that combines SSR and BEDPP (Basic EDPP)

## Installation:
* the latest released version: 
```R
install.packages("grpreg")
```

* the latest version (requires `devtools`): 
```R
install_github("pbreheny/grpreg")
```

## Reference:
* Tibshirani, R., Bien, J., Friedman, J., Hastie, T., Simon, N., Taylor, J., and Tibshirani, R. J. (2012). Strong rules for discarding predictors in lasso-type problems. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 74 (2), 245-266.
* Wang, J., Zhou, J., Liu, J., Wonka, P., and Ye, J. (2014). A safe screening rule for sparse logistic regression. In Advances in Neural Information Processing Systems, pp. 1053-1061.




