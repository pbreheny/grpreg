# Generate nonlinear example data

Mainly intended to demonstrate the use of basis expansion models for
sparse additive modeling; intended for use with
[`expand_spline()`](https://pbreheny.github.io/grpreg/reference/expand_spline.md).

## Usage

``` r
gen_nonlinear_data(n = 100, p = 16, seed)
```

## Arguments

- n:

  Sample size (numeric; default = 100).

- p:

  Number of features (numeric; default = 16).

- seed:

  Set to get different random data sets, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html)

## Examples

``` r
Data <- gen_nonlinear_data()
```
