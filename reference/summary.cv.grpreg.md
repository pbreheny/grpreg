# Summarizing inferences based on cross-validation

Summary method for `cv.grpreg` or `cv.grpsurv` objects

## Usage

``` r
# S3 method for class 'cv.grpreg'
summary(object, ...)

# S3 method for class 'summary.cv.grpreg'
print(x, digits, ...)
```

## Arguments

- object:

  A `"cv.grpreg"` object.

- ...:

  Further arguments passed to or from other methods.

- x:

  A `"summary.cv.grpreg"` object.

- digits:

  Number of digits past the decimal point to print out. Can be a vector
  specifying different display digits for each of the five non-integer
  printed values.

## Value

`summary(cvfit)` produces an object with S3 class `"summary.cv.grpreg"`.
The class has its own print method and contains the following list
elements:

- penalty:

  The penalty used by `grpreg`/`grpsurv`.

- model:

  The type of model: `"linear"`, `"logistic"`, `"Poisson"`, `"Cox"`,
  etc.

- n:

  Number of observations

- p:

  Number of regression coefficients (not including the intercept).

- min:

  The index of `lambda` with the smallest cross-validation error.

- lambda:

  The sequence of `lambda` values used by `cv.grpreg`/`cv.grpsurv`.

- cve:

  Cross-validation error (deviance).

- r.squared:

  Proportion of variance explained by the model, as estimated by
  cross-validation.

- snr:

  Signal to noise ratio, as estimated by cross-validation.

- sigma:

  For linear regression models, the scale parameter estimate.

- pe:

  For logistic regression models, the prediction error
  (misclassification error).

## See also

[`grpreg`](https://pbreheny.github.io/grpreg/reference/grpreg.md),
[`cv.grpreg`](https://pbreheny.github.io/grpreg/reference/cv.grpreg.md),
[`cv.grpsurv`](https://pbreheny.github.io/grpreg/reference/cv.grpreg.md),
[`plot.cv.grpreg`](https://pbreheny.github.io/grpreg/reference/plot.cv.grpreg.md)

## Author

Patrick Breheny

## Examples

``` r
# Birthweight data
data(Birthwt)
X <- Birthwt$X
group <- Birthwt$group

# Linear regression
y <- Birthwt$bwt
cvfit <- cv.grpreg(X, y, group)
summary(cvfit)
#> grLasso-penalized linear regression with n=189, p=16
#> At minimum cross-validation error (lambda=0.0167):
#> -------------------------------------------------
#>   Nonzero coefficients: 16
#>   Nonzero groups: 8
#>   Cross-validation error of 0.43
#>   Maximum R-squared: 0.19
#>   Maximum signal-to-noise ratio: 0.23
#>   Scale estimate (sigma) at lambda.min: 0.654

# Logistic regression
y <- Birthwt$low
cvfit <- cv.grpreg(X, y, group, family="binomial")
summary(cvfit)
#> grLasso-penalized logistic regression with n=189, p=16
#> At minimum cross-validation error (lambda=0.0180):
#> -------------------------------------------------
#>   Nonzero coefficients: 16
#>   Nonzero groups: 8
#>   Cross-validation error of 1.17
#>   Maximum R-squared: 0.07
#>   Maximum signal-to-noise ratio: 0.06
#>   Prediction error at lambda.min: 0.312

# Cox regression
data(Lung)
cvfit <- with(Lung, cv.grpsurv(X, y, group))
summary(cvfit)
#> grLasso-penalized Cox regression with n=137, p=14
#> At minimum cross-validation error (lambda=0.1063):
#> -------------------------------------------------
#>   Nonzero coefficients: 7
#>   Nonzero groups: 2
#>   Cross-validation error of 7.60
#>   Maximum R-squared: 0.25
#>   Maximum signal-to-noise ratio: 0.04
```
