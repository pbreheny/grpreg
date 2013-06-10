test_that("grLasso reproduces lasso", {
  require(glmnet)
  n <- 50
  group <- 1:10
  p <- length(group)
  X <- matrix(rnorm(n*p),ncol=p)
  y <- rnorm(n)
  yy <- runif(n) > .5
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso"))
  par(mfrow=c(2,2)); plot(fit, log=TRUE)
  lasso <- as.matrix(coef(fit <- glmnet(X, y, lambda=fit$lambda)))
  plot(fit, "lambda")
  expect_that(grLasso, equals(lasso, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial"))
  plot(fit, log=TRUE)
  lasso <- as.matrix(coef(fit <- glmnet(X, yy, family="binomial", lambda=fit$lambda)))
  plot(fit, "lambda")  
  expect_that(grLasso, equals(lasso, tolerance=.01, check.attributes=FALSE))
})

test_that("grMCP and grSCAD reproduce MCP and SCAD", {
  require(ncvreg)
  n <- 50
  group <- 1:10
  p <- length(group)
  X <- matrix(rnorm(n*p),ncol=p)
  y <- rnorm(n)
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP"))
  par(mfrow=c(2,2)); plot(fit)
  mcp <- coef(fit <- ncvreg(X, y, lambda=fit$lambda, penalty="MCP"))
  plot(fit)
  expect_that(grMCP, equals(mcp, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD"))
  plot(fit)
  scad <- coef(ncvreg(X, y, lambda=fit$lambda, penalty="SCAD"))
  plot(fit)
  expect_that(grSCAD, equals(scad, tolerance=.01, check.attributes=FALSE))
})

test_that("gel reproduces lasso", {
  require(glmnet)
  n <- 100
  group <- rep(1,10)
  p <- length(group)
  X <- matrix(rnorm(n*p),ncol=p)
  y <- rnorm(n)
  yy <- rnorm(n) > 0
  gel <- coef(fit <- grpreg(X, y, group, penalty="gel", tau=0))
  par(mfrow=c(2,2)); plot(fit, log=TRUE)
  lasso <- as.matrix(coef(fit <- glmnet(X, y, lambda=fit$lambda)))
  plot(fit, "lambda")
  expect_that(gel, equals(lasso, tolerance=.01, check.attributes=FALSE))
  gel <- coef(fit <- grpreg(X, yy, group, penalty="gel", family="binomial", tau=0))
  plot(fit, log=TRUE)
  lasso <- as.matrix(coef(fit <- glmnet(X, yy, family="binomial", lambda=fit$lambda)))
  plot(fit, "lambda")  
  expect_that(gel, equals(lasso, tolerance=.01, check.attributes=FALSE))
})
