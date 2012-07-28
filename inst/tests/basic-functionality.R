test_that("grpreg() reproduces simple linear regression", {
  n <- 5
  p <- 1
  X <- matrix(rnorm(n*p), ncol=p)
  y <- rnorm(n)
  group <- 1
  reg <- lm(y~X)$coef
  nlam=100
  par(mfrow=c(3,2))
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", nlambda=nlam, lambda.min=0))[,nlam]
  plot(fit, main=fit$penalty)
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  bridge <- coef(fit <- gBridge(X, y, group, nlambda=nlam, lambda.min=0))[,1]
  plot(fit, main=fit$penalty)
  expect_that(bridge, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", nlambda=nlam, lambda.min=0))[,nlam]
  plot(fit, main=fit$penalty)
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", nlambda=nlam, lambda.min=0))[,nlam]
  plot(fit, main=fit$penalty)
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", nlambda=nlam, lambda.min=0))[,nlam]
  plot(fit, main=fit$penalty)
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
})

test_that("grpreg() reproduces linear regression", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  y <- rnorm(n)
  group <- rep(0:3,4:1)
  fit.mle <- lm(y~X)
  reg <- coef(fit.mle)
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", lambda.min=0))[,100]
  if (interactive()) {par(mfrow=c(3,2));plot(fit, main=fit$penalty)}
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  bridge <- coef(fit <- gBridge(X, y, group, lambda.min=0))[,1]
  if (interactive()) {plot(fit, main=fit$penalty)}
  expect_that(bridge, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0))[,100]
  if (interactive()) {plot(fit, main=fit$penalty)}
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", lambda.min=0))[,100]
  if (interactive()) plot(fit, main=fit$penalty)
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", lambda.min=0))[,100]
  if (interactive()) plot(fit, main=fit$penalty)
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
  expect_that(predict(fit, X)[,100], equals(predict(fit.mle), tolerance=.001, check.attributes=FALSE))
})

test_that("grpreg() reproduces simple logistic regression", {
  n <- 20
  p <- 1
  X <- matrix(rnorm(n*p), ncol=p)
  y <- runif(n) > .5
  group <- 1
  reg <- glm(y~X, family="binomial")$coef
  nlam <- 100
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", nlambda=nlam, lambda.min=0, family="binomial", gamma=12))[,nlam]
  if (interactive()) {par(mfrow=c(3,2));plot(fit, main=fit$penalty)}
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  bridge <- coef(fit <- gBridge(X, y, group, lambda.min=0, nlambda=nlam, family="binomial"))[,1]
  if (interactive()) {plot(fit, main=fit$penalty)}
  expect_that(bridge, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", nlambda=nlam, lambda.min=0, family="binomial"))[,nlam]
  if (interactive()) {plot(fit, main=fit$penalty)}
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", nlambda=nlam, lambda.min=0, family="binomial"))[,nlam]
  if (interactive()) {plot(fit, main=fit$penalty)}
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", nlambda=nlam, lambda.min=0, family="binomial"))[,nlam]
  if (interactive()) {plot(fit, main=fit$penalty)}
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
})

test_that("grpreg() reproduces logistic regression", {
  n <- 50
  group <- rep(0:3,1:4)
  ## group <- rep(1:3,rep(3,3))
  p <- length(group)
  X <- matrix(rnorm(n*p),ncol=p)
  y <- runif(n) > .5
  fit.mle <- glm(y~X, family="binomial")
  reg <- coef(fit.mle)
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", family="binomial"))[,100]
  if (interactive()) {par(mfrow=c(3,2));plot(fit, main=fit$penalty)}
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  bridge <- coef(fit <- gBridge(X, y, group, family="binomial"))[,1]
  if (interactive()) plot(fit, main=fit$penalty)
  expect_that(bridge, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", family="binomial"))[,100]
  if (interactive()) plot(fit, main=fit$penalty)
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", family="binomial", gamma=2))[,100]
  if (interactive()) plot(fit, main=fit$penalty)
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", family="binomial", gamma=2))[,100]
  if (interactive()) plot(fit, main=fit$penalty)
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
  expect_that(predict(fit, X)[,100], equals(predict(fit.mle), tolerance=.001, check.attributes=FALSE))
  expect_that(predict(fit, X, type="response")[,100], equals(predict(fit.mle, type="response"), tolerance=.001, check.attributes=FALSE))
})
