# grpreg: Regularization Paths for Regression Models with Grouped Covariates

Efficient algorithms for fitting the regularization path of linear
regression, GLM, and Cox regression models with grouped penalties. This
includes group selection methods such as group lasso, group MCP, and
group SCAD as well as bi-level selection methods such as the group
exponential lasso, the composite MCP, and the group bridge. For more
information, see Breheny and Huang (2009)
[doi:10.4310/sii.2009.v2.n3.a10](https://doi.org/10.4310/sii.2009.v2.n3.a10)
, Huang, Breheny, and Ma (2012)
[doi:10.1214/12-sts392](https://doi.org/10.1214/12-sts392) , Breheny and
Huang (2015)
[doi:10.1007/s11222-013-9424-2](https://doi.org/10.1007/s11222-013-9424-2)
, and Breheny (2015)
[doi:10.1111/biom.12300](https://doi.org/10.1111/biom.12300) , or visit
the package homepage <https://pbreheny.github.io/grpreg/>.

## References

- Yuan M and Lin Y. (2006) Model selection and estimation in regression
  with grouped variables. *Journal of the Royal Statistical Society
  Series B*, **68**: 49-67.
  [doi:10.1111/j.1467-9868.2005.00532.x](https://doi.org/10.1111/j.1467-9868.2005.00532.x)

- Huang J, Ma S, Xie H, and Zhang C. (2009) A group bridge approach for
  variable selection. *Biometrika*, **96**: 339-355.
  [doi:10.1093/biomet/asp020](https://doi.org/10.1093/biomet/asp020)

- Breheny P and Huang J. (2009) Penalized methods for bi-level variable
  selection. *Statistics and its interface*, **2**: 369-380.
  [doi:10.4310/sii.2009.v2.n3.a10](https://doi.org/10.4310/sii.2009.v2.n3.a10)

- Huang J, Breheny P, and Ma S. (2012). A selective review of group
  selection in high dimensional models. *Statistical Science*, **27**:
  481-499. [doi:10.1214/12-sts392](https://doi.org/10.1214/12-sts392)

- Breheny P and Huang J. (2015) Group descent algorithms for nonconvex
  penalized linear and logistic regression models with grouped
  predictors. *Statistics and Computing*, **25**: 173-187.
  [doi:10.1007/s11222-013-9424-2](https://doi.org/10.1007/s11222-013-9424-2)

- Breheny P. (2015) The group exponential lasso for bi-level variable
  selection. *Biometrics*, **71**: 731-740.
  [doi:10.1111/biom.12300](https://doi.org/10.1111/biom.12300)

## See also

Useful links:

- <https://pbreheny.github.io/grpreg/>

- <https://github.com/pbreheny/grpreg>

- Report bugs at <https://github.com/pbreheny/grpreg/issues>

## Author

Patrick Breheny

## Examples

``` r
vignette("getting-started", package="grpreg")
#> Warning: vignette ‘getting-started’ not found
```
