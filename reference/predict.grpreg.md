# Model predictions based on a fitted `grpreg` object

Similar to other predict methods, this function returns predictions from
a fitted `"grpreg"` object.

## Usage

``` r
# S3 method for class 'cv.grpreg'
predict(
  object,
  X,
  lambda = object$lambda.min,
  which = object$min,
  type = c("link", "response", "class", "coefficients", "vars", "groups", "nvars",
    "ngroups", "norm"),
  ...
)

# S3 method for class 'cv.grpreg'
coef(object, lambda = object$lambda.min, which = object$min, ...)

# S3 method for class 'grpreg'
predict(
  object,
  X,
  type = c("link", "response", "class", "coefficients", "vars", "groups", "nvars",
    "ngroups", "norm"),
  lambda,
  which = 1:length(object$lambda),
  ...
)

# S3 method for class 'grpreg'
coef(object, lambda, which = 1:length(object$lambda), drop = TRUE, ...)
```

## Arguments

- object:

  Fitted `"grpreg"` or `"cv.grpreg"` model object.

- X:

  Matrix of values at which predictions are to be made. Not used for
  `type="coefficients"`

- lambda:

  Values of the regularization parameter `lambda` at which predictions
  are requested. For values of `lambda` not in the sequence of fitted
  models, linear interpolation is used.

- which:

  Indices of the penalty parameter `lambda` at which predictions are
  required. By default, all indices are returned. If `lambda` is
  specified, this will override `which`.

- type:

  Type of prediction: `"link"` returns the linear predictors;
  `"response"` gives the fitted values; `"class"` returns the binomial
  outcome with the highest probability; `"coefficients"` returns the
  coefficients; `"vars"` returns the indices for the nonzero
  coefficients; `"groups"` returns the indices for the groups with at
  least one nonzero coefficient; `"nvars"` returns the number of nonzero
  coefficients; `"ngroups"` returns the number of groups with at least
  one nonzero coefficient; `"norm"` returns the L2 norm of the
  coefficients in each group.

- ...:

  Not used.

- drop:

  By default, if a single value of `lambda` is supplied, a vector of
  coefficients is returned. Set `drop=FALSE` if you wish to have `coef`
  always return a matrix (see
  [`drop`](https://rdrr.io/r/base/drop.html)).

## Value

The object returned depends on type.

## Details

`coef` and `predict` methods are provided for `"cv.grpreg"` options as a
convenience. They simply call `coef.grpreg` and `predict.grpreg` with
`lambda` set to the value that minimizes the cross-validation error.

## See also

`grpreg`

## Author

Patrick Breheny

## Examples

``` r
# Fit penalized logistic regression model to birthweight data
data(Birthwt)
X <- Birthwt$X
y <- Birthwt$low
group <- Birthwt$group
fit <- grpreg(X, y, group, penalty="grLasso", family="binomial")

# Coef and predict methods
coef(fit, lambda=.001)
#> (Intercept)        age1        age2        age3        lwt1        lwt2 
#>  -1.5594452 -11.3887695 -18.3241083 -13.5534407  -6.7400536  -2.0424956 
#>        lwt3       white       black       smoke        ptl1       ptl2m 
#>  -4.4466007  -0.6955048   0.5156693   0.8198631   1.6532574  -0.2816916 
#>          ht          ui        ftv1        ftv2       ftv3m 
#>   2.0050588   0.7726613  -0.3866688  -0.1593881   0.6674573 
predict(fit, X, type="link", lambda=.07)[1:10]
#>  [1] -0.8193011 -0.8702524 -0.8639390 -0.8129877 -0.8129877 -0.8702524
#>  [7] -0.8702524 -0.8702524 -0.8639390 -0.8639390
predict(fit, X, type="response", lambda=.07)[1:10]
#>  [1] 0.3059120 0.2952018 0.2965170 0.3072542 0.3072542 0.2952018 0.2952018
#>  [8] 0.2952018 0.2965170 0.2965170
predict(fit, X, type="class", lambda=.01)[1:15]
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0
predict(fit, type="vars", lambda=.07)
#> smoke  ptl1 ptl2m    ht    ui 
#>     9    10    11    12    13 
predict(fit, type="groups", lambda=.07)
#> [1] smoke ptl   ht    ui   
#> Levels: age lwt race smoke ptl ht ui ftv
predict(fit, type="norm", lambda=.07)
#>         age         lwt        race       smoke         ptl          ht 
#> 0.000000000 0.000000000 0.000000000 0.006313417 0.487194445 0.032052407 
#>          ui         ftv 
#> 0.050951339 0.000000000 

# Coef and predict methods for cross-validation
cvfit <- cv.grpreg(X, y, group, family="binomial", penalty="grMCP")
coef(cvfit)
#> (Intercept)        age1        age2        age3        lwt1        lwt2 
#> -1.33722466  0.00000000  0.00000000  0.00000000 -7.87460993 -1.94199655 
#>        lwt3       white       black       smoke        ptl1       ptl2m 
#> -3.95042561 -0.80627119  0.41880008  0.82100516  1.45766647 -0.23855135 
#>          ht          ui        ftv1        ftv2       ftv3m 
#>  1.91824312  0.73405273 -0.04308188 -0.01391116  0.04301617 
predict(cvfit, X)[1:10]
#>  [1] -0.6065848 -1.5010329 -0.9905203 -0.3142638 -0.2725193 -1.3597704
#>  [7] -2.1194819 -0.9413244 -1.3750892 -1.1572809
predict(cvfit, X, type="response")[1:10]
#>  [1] 0.3528386 0.1822715 0.2708093 0.4220743 0.4322887 0.2042776 0.1072177
#>  [8] 0.2806329 0.2017989 0.2391617
predict(cvfit, type="groups")
#> [1] lwt   race  smoke ptl   ht    ui    ftv  
#> Levels: age lwt race smoke ptl ht ui ftv
```
