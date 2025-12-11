# Package index

## Model fitting

- [`grpreg()`](https://pbreheny.github.io/grpreg/reference/grpreg.md) :
  Fit a group penalized regression path
- [`grpsurv()`](https://pbreheny.github.io/grpreg/reference/grpsurv.md)
  : Fit an group penalized survival model
- [`gBridge()`](https://pbreheny.github.io/grpreg/reference/gBridge.md)
  : Fit a group bridge regression path

## Model selection and cross-validation

- [`cv.grpreg()`](https://pbreheny.github.io/grpreg/reference/cv.grpreg.md)
  [`cv.grpsurv()`](https://pbreheny.github.io/grpreg/reference/cv.grpreg.md)
  : Cross-validation for grpreg/grpsurv

- [`plot(`*`<cv.grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/plot.cv.grpreg.md)
  :

  Plots the cross-validation curve from a `cv.grpreg` object

- [`summary(`*`<cv.grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/summary.cv.grpreg.md)
  [`print(`*`<summary.cv.grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/summary.cv.grpreg.md)
  : Summarizing inferences based on cross-validation

- [`AUC()`](https://pbreheny.github.io/grpreg/reference/AUC.cv.grpsurv.md)
  : Calculates AUC for cv.grpsurv objects

- [`select()`](https://pbreheny.github.io/grpreg/reference/select.md) :
  Select an value of lambda along a grpreg path

## Plotting and extracting model features

- [`logLik(`*`<grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/logLik.grpreg.md)
  [`logLik(`*`<grpsurv>`*`)`](https://pbreheny.github.io/grpreg/reference/logLik.grpreg.md)
  : logLik method for grpreg

- [`predict(`*`<cv.grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/predict.grpreg.md)
  [`coef(`*`<cv.grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/predict.grpreg.md)
  [`predict(`*`<grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/predict.grpreg.md)
  [`coef(`*`<grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/predict.grpreg.md)
  :

  Model predictions based on a fitted `grpreg` object

- [`predict(`*`<grpsurv>`*`)`](https://pbreheny.github.io/grpreg/reference/predict.grpsurv.md)
  : Model predictions for grpsurv objects

- [`plot(`*`<grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/plot.grpreg.md)
  : Plot coefficients from a "grpreg" object

- [`plot(`*`<grpsurv.func>`*`)`](https://pbreheny.github.io/grpreg/reference/plot.grpsurv.func.md)
  : Plot survival curve for grpsurv model

- [`residuals(`*`<grpreg>`*`)`](https://pbreheny.github.io/grpreg/reference/residuals.grpreg.md)
  : Extract residuals from a grpreg or grpsurv fit

## Utilities for additive models

- [`expand_spline()`](https://pbreheny.github.io/grpreg/reference/expand_spline.md)
  : Expand feature matrix using basis splines
- [`gen_nonlinear_data()`](https://pbreheny.github.io/grpreg/reference/gen_nonlinear_data.md)
  : Generate nonlinear example data
- [`plot_spline()`](https://pbreheny.github.io/grpreg/reference/plot_spline.md)
  : Plot spline curve for a fitted additive model

## Data sets

- [`Birthwt`](https://pbreheny.github.io/grpreg/reference/Birthwt.md) :
  Risk Factors Associated with Low Infant Birth Weight
- [`Lung`](https://pbreheny.github.io/grpreg/reference/Lung.md) : VA
  lung cancer data set
