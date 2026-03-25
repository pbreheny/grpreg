# Marginal false discovery rates

Estimates the marginal false discovery rate (mFDR) of a group penalized
regression model.

## Usage

``` r
mfdr(fit, X)
```

## Arguments

- fit:

  A `grpreg` or `grpsurv` object.

- X:

  The model matrix corresponding to `fit`. This is not necessary for
  linear regression, but in logistic and Cox regression, the mFDR
  depends on X. It is not necessary to supply `X` if it is already
  contained in `fit`; i.e., if `ncvreg`/`ncvsurv` was run with
  `returnX = TRUE`.

## Value

An object with S3 class `mfdr` inheriting from `data.frame`, containing:

- ef:

  The number of variables selected at each value of `lambda`, averaged
  over the permutation fits.

- s:

  The actual number of selected variables for the non-permuted data.

- mfdr:

  The estimated marginal false discovery rate (`ef/s`).

## Details

The function estimates the marginal false discovery rate (mFDR) for
groups in a group lasso or group MCP penalized regression model. The
estimate tends to be accurate in most settings, but will be somewhat
conservative if predictors are highly correlated.

## See also

[`grpreg()`](https://pbreheny.github.io/grpreg/reference/grpreg.md),
[`grpsurv()`](https://pbreheny.github.io/grpreg/reference/grpsurv.md)

## Examples

``` r
# Birthweight data ---------------------------
data(Birthwt)
x <- Birthwt$X
group <- Birthwt$group

# Linear regression --------------------------
y <- Birthwt$bwt
fit <- grpreg(x, y, group)
obj <- mfdr(fit)
head(obj)
#>                   ef s        mfdr
#> 0.206495 0.000000000 0 0.000000000
#> 0.188151 0.001073303 1 0.001073303
#> 0.171436 0.003251466 1 0.003251466
#> 0.156206 0.008496796 1 0.008496796
#> 0.142329 0.019651242 1 0.019651242
#> 0.129685 0.040878376 2 0.020439188

# Logistic regression ------------------------------
y <- Birthwt$low
fit <- grpreg(x, y, group, penalty="grMCP", family="binomial", returnX = TRUE)
obj <- mfdr(fit)
# If returnX is not TRUE, user must supply X
fit <- grpreg(x, y, group, penalty="grMCP", family="binomial")
obj <- mfdr(fit, x)
head(obj)
#>                   ef s       mfdr
#> 0.0960555 0.00000000 0 0.00000000
#> 0.0875222 0.03048445 1 0.03048445
#> 0.0797470 0.06129080 1 0.06129080
#> 0.0726625 0.11468468 1 0.11468468
#> 0.0662073 0.20446253 3 0.06815418
#> 0.0603256 0.34473624 4 0.08618406

# Cox regression -----------------------------------
data(Lung)
x_lung <- Lung$X
y_lung <- Lung$y
g_lung <- Lung$group
fit <- grpsurv(x_lung, y_lung, g_lung, penalty = "grSCAD")
obj <- mfdr(fit, x_lung)
head(obj)
#>                 ef s        mfdr
#> 0.2632 0.000000000 0 0.000000000
#> 0.2455 0.005998595 1 0.005998595
#> 0.2290 0.012415206 1 0.012415206
#> 0.2135 0.024215529 1 0.024215529
#> 0.1991 0.044357731 2 0.022178865
#> 0.1857 0.077020234 2 0.038510117
```
