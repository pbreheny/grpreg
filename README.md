`grpreg` fits regularization paths for linear, logistic, or Poisson regression models with grouped penalties, such as the group lasso, group MCP, group SCAD, group exponential lasso, and group bridge. The algorithms are based on the idea of either locally approximated coordinate descent or group descent, depending on the penalty. All of the algorithms (with the exception of group bridge) are stable and fast.

To install:

* the latest released version: `install.packages("grpreg")`
* the latest version (requires `devtools`): `install_github("pbreheny/grpreg")`
