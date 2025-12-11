# Penalties

**grpreg** fits models that fall into the penalized likelihood
framework, in which we estimate \boldsymbol{\beta} by minimizing the
objective function

Q(\boldsymbol{\beta}\|\mathbf{X}, \mathbf{y}) =
L(\boldsymbol{\beta}\|\mathbf{X},\mathbf{y}) +
P\_\lambda(\boldsymbol{\beta}),

where L(\boldsymbol{\beta}\|\mathbf{X},\mathbf{y}) is the loss
(deviance) and P\_\lambda(\boldsymbol{\beta}) is the penalty. This
article describes the different penalties available in `grpreg`; see
[models](https://pbreheny.github.io/grpreg/articles/models.md) for more
information on the different loss functions available.

The following notation is used throughout (recall that the design matrix
\mathbf{X} is decomposed into groups \mathbf{X}\_1, \mathbf{X}\_2,
\ldots:

- \boldsymbol{\beta} denotes the entire vector of regression
  coefficients
- \boldsymbol{\beta}\_j denotes the vector of regression coefficients
  corresponding to the jth group
- \beta\_{jk} denotes kth regression coefficient in the jth group
- \lVert\boldsymbol{\beta}\_j\rVert_2 denotes the Euclidean (L_2) norm
  of \boldsymbol{\beta}\_j: \lVert x\rVert_2 = \sqrt{x_1^2 + x_2^2 +
  \ldots}
- \lVert\boldsymbol{\beta}\_j\rVert_1 denotes the L_1 norm of \beta_j:
  \lVert x\rVert_1 = \left\lvert x_1\right\rvert + \left\lvert
  x_2\right\rvert + \ldots

## Group selection

These penalties are sparse at the group level – the coefficients within
a group will either all equal zero or none will equal zero.

If you use any of these penalties, please cite

- Breheny P and Huang J (2015). Group descent algorithms for nonconvex
  penalized linear and logistic regression models with grouped
  predictors. *Statistics and Computing*, **25**: 173-187.
  \[[pdf](https://myweb.uiowa.edu/pbreheny/pdf/group-computing.pdf)\].

The article goes into more mathematical details, discusses issues of
standardization in the group sense, and provides references.

The group lasso was originally proposed in

- Yuan M. and Lin Y. (2006) Model selection and estimation in regression
  with grouped variables. *Journal of the Royal Statistical Society
  Series B*, **68**: 49-67.

### Group lasso

``` r
grpreg(X, y, group, penalty="grLasso")
```

P(\beta) = \lambda\sum_j \lVert\boldsymbol{\beta}\_j\rVert_2

### Group MCP

``` r
grpreg(X, y, group, penalty="grMCP")
```

P(\boldsymbol{\beta}) = \sum_j \textrm{MCP}\_{\lambda,
\gamma}(\lVert\boldsymbol{\beta}\_j\rVert_2)

where \textrm{MCP}\_{\lambda, \gamma}(\cdot) denotes the MCP penalty
with regularization parameter \lambda and tuning parameter \gamma.

### Group SCAD

``` r
grpreg(X, y, group, penalty="grSCAD")
```

P(\boldsymbol{\beta}) = \sum_j \textrm{SCAD}\_{\lambda,
\gamma}(\lVert\boldsymbol{\beta}\_j\rVert_2)

where \textrm{SCAD}\_{\lambda, \gamma}(\cdot) denotes the SCAD penalty
with regularization parameter \lambda and tuning parameter \gamma.

## Bi-level selection

These penalties are sparse at both the group and individual levels. In
some groups, all coefficients will equal zero. However, even if a group
is selected, some of the coefficients within that group may still be
zero.

### Group exponential lasso (GEL)

``` r
grpreg(X, y, group, penalty="gel")
```

P(\beta) = \sum_j f\_{\lambda,
\tau}(\lVert\boldsymbol{\beta}\_j\rVert_1)

where f(\cdot) denotes the exponential penalty with regularization
parameter \lambda and tuning parameter \tau:

f\_{\lambda, \tau}(\theta) =
\frac{\lambda^2}{\tau}\left\\1-\exp\left(-\frac{\tau\theta}{\lambda}\right)\right\\

If you use the GEL penalty, please cite

- Breheny P (2015). The group exponential lasso for bi-level variable
  selection. *Biometrics*, **71**: 731-740.
  \[[pdf](https://onlinelibrary.wiley.com/doi/10.1111/biom.12300/pdf)\].

### Composite MCP

``` r
grpreg(X, y, group, penalty="cMCP")
```

P(\boldsymbol{\beta}) = \sum_j \textrm{MCP}\_{\lambda, \gamma_1} \left(
\sum_k \textrm{MCP}\_{\lambda, \gamma_2}
(\left\lvert\beta\_{jk}\right\rvert) \right)

where \textrm{MCP}\_{\lambda, \gamma}(\cdot) denotes the MCP penalty
with regularization parameter \lambda and tuning parameter \gamma.

If you use the composite MCP penalty, please cite either of the
following papers:

- Breheny P and Huang J (2009). Penalized methods for bi-level variable
  selection. *Statistics and Its Interface*, **2**: 369-380.
  \[[pdf](https://myweb.uiowa.edu/pbreheny/pdf/Breheny2009.pdf)\]
- Huang J, Breheny P and Ma S (2012). A selective review of group
  selection in high-dimensional models. *Statistical Science*, **27**:
  481-499. \[[pdf](https://myweb.uiowa.edu/pbreheny/pdf/Huang2012.pdf)\]

Please note that there is some confusion around the name “group MCP”. In
the first paper above (2009), the composite MCP penalty was referred to
as the “group MCP” penalty; the second paper (2012), in reviewing the
various kinds of group penalties that had been proposed, recommended
changing the name to “composite MCP” to avoid confusion with the “group
MCP” [defined above](#group-mcp).

### Group bridge

``` r
gBridge(X, y, group)
```

P(\boldsymbol{\beta}) = \lambda \sum_j K_j^\gamma
\lVert\boldsymbol{\beta}\_j\rVert_1^\gamma

where K_j denotes the number of elements in group j.

Please note that the group bridge penalty uses a very different
algorithm from the other penalties. Due to the nature of the penalty,
model fitting is slower and less stable for group bridge models. This
is, in fact, the main motivation of the GEL penalty of Section~: to
offer a more tractable alternative to group bridge that has similar
estimation properties but is much better behaved from a numerical
optimization perspective.

If you use the group bridge penalty, please cite either of the following
papers:

- Huang J, Ma S, Xie H and Zhang C (2009). A group bridge approach for
  variable selection. *Biometrika*, **96**: 339-355.
- Breheny P and Huang J (2009). Penalized methods for bi-level variable
  selection. *Statistics and Its Interface*, **2**: 369-380.
  \[[pdf](NA)\]

The first paper proposed the method; the second paper proposed the
algorithm that is used in the `grpreg` package.

## Specifying an additional ridge component

For all of the penalties in the previous section, `grpreg` allows the
specification of an additional ridge (L_2) component to the penalty.
This will set \lambda_1 = \alpha\lambda and \lambda_2=(1-\alpha)\lambda,
with the penalty given by

P(\boldsymbol{\beta}) = P_1(\boldsymbol{\beta}\|\lambda_1) +
\frac{\lambda_2}{2}\lVert\boldsymbol{\beta}\rVert_2^2,

where P_1 is any of the penalties from the earlier sections. So, for
example

``` r
grpreg(X, y, group, penalty="grLasso", alpha=0.75)
```

will fit a model with penalty

P(\beta) = 0.75\lambda\sum_j \lVert\boldsymbol{\beta}\_j\rVert_2 +
\frac{0.25\lambda}{2}\lVert\boldsymbol{\beta}\rVert_2^2.
