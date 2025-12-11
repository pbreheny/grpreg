# Extract residuals from a grpreg or grpsurv fit

Currently, only deviance residuals are supported.

## Usage

``` r
# S3 method for class 'grpreg'
residuals(object, lambda, which = 1:length(object$lambda), drop = TRUE, ...)
```

## Arguments

- object:

  Object of class `grpreg` or `grpsurv`.

- lambda:

  Values of the regularization parameter at which residuals are
  requested (numeric vector). For values of lambda not in the sequence
  of fitted models, linear interpolation is used.

- which:

  Index of the penalty parameter at which residuals are requested
  (default = all indices). If `lambda` is specified, this take
  precedence over `which`.

- drop:

  By default, if a single value of lambda is supplied, a vector of
  residuals is returned (logical; default=`TRUE`). Set `drop=FALSE` if
  you wish to have the function always return a matrix (see
  [`drop()`](https://rdrr.io/r/base/drop.html)).

- ...:

  Not used.

## Examples

``` r
data(Birthwt)
X <- Birthwt$X
y <- Birthwt$bwt
group <- Birthwt$group
fit <- grpreg(X, y, group, returnX=TRUE)
residuals(fit)[1:5, 1:5]
#>          0.2065     0.1882    0.1714     0.1562     0.1423
#> [1,] -0.4215873 -0.3775988 -0.337518 -0.3009980 -0.2677223
#> [2,] -0.3935873 -0.4012375 -0.408208 -0.4145594 -0.4203464
#> [3,] -0.3875873 -0.3952375 -0.402208 -0.4085594 -0.4143464
#> [4,] -0.3505873 -0.3065988 -0.266518 -0.2299980 -0.1967223
#> [5,] -0.3445873 -0.3005988 -0.260518 -0.2239980 -0.1907223
head(residuals(fit, lambda=0.1))
#> [1] -0.17267168 -0.45110465 -0.41799818 -0.08848411 -0.08248411 -0.38010465
```
