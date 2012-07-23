## source("../../../.debug-grpreg.R")

test_that("standardize() standardizes correctly", {
  X <- matrix(rnorm(500),ncol=10)
  XX <- standardize(X)
  expect_that(apply(XX,2,mean), equals(rep(0,10)))
  expect_that(apply(XX,2,crossprod), equals(rep(50,10)))
})

test_that("unstandardize() unstandardizes correctly", {
  n <- 50
  p <- 10
  l <- 5
  X <- matrix(rnorm(n*p),ncol=p)
  XX <- standardize(X)
  bb <- matrix(rnorm(l*(p+1)), nrow=p+1)
  b <- unstandardize(bb, attr(XX,"center"), attr(XX, "scale"))
  expect_that(cbind(1,XX) %*% bb, equals(cbind(1,X) %*% b))
})

test_that("orthogonalize() orthogonalizes correctly", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  group <- c(1,1,2,2,2,3,3,3,3,4)
  XX <- orthogonalize(X, group)
  for (j in 1:group[p]) {
    ind <- which(group==j)
    expect_that(crossprod(XX[,ind])/n, equals(diag(length(ind))))
  }
})

test_that("unorthogonalize() unorthogonalizes correctly", {
  n <- 50
  p <- 10
  l <- 5
  group <- c(1,1,2,2,2,3,3,3,3,4)
  X <- matrix(rnorm(n*p), ncol=p)
  XX <- orthogonalize(X, group)
  bb <- matrix(rnorm(l*(p+1)), nrow=p+1)
  b <- unorthogonalize(bb, XX, group)
  expect_that(cbind(1,XX) %*% bb, equals(cbind(1,X) %*% b))
})

test_that("orthogonalize() orthogonalizes correctly w/ 0's present", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  group <- rep(0:3,4:1)
  XX <- orthogonalize(X, group)
  for (j in 1:group[p]) {
    ind <- which(group==j)
    expect_that(crossprod(XX[,ind])/n, equals(diag(length(ind))))
  }
})

test_that("grpreg() reproduces simple linear regression", {
  n <- 5
  p <- 1
  X <- matrix(rnorm(n*p), ncol=p)
  y <- rnorm(n)
  group <- 1
  reg <- lm(y~X)$coef
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", nlambda=10, lambda.min=0))[,10]
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  bridge <- coef(fit <- gBridge(X, y, group, nlambda=10, lambda.min=0))[,1]
  expect_that(bridge, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", nlambda=10, lambda.min=0))[,10]
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", nlambda=10, lambda.min=0))[,10]
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", nlambda=10, lambda.min=0))[,10]
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
})

test_that("grpreg() reproduces linear regression", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  y <- rnorm(n)
  group <- rep(1:4,1:4)
  fit.mle <- lm(y~X)
  reg <- coef(fit.mle)
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", lambda.min=0))[,100]
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  bridge <- coef(fit <- gBridge(X, y, group, lambda.min=0))[,1]
  expect_that(bridge, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0))[,100]
  if (interactive()) {par(mfrow=c(3,1));plot(fit, main=fit$penalty)}
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
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", nlambda=10, lambda.min=0, family="binomial"))[,10]
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  bridge <- coef(fit <- gBridge(X, y, group, lambda.min=0, nlambda=10, family="binomial"))[,1]
  expect_that(bridge, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", nlambda=10, lambda.min=0, family="binomial"))[,10]
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", nlambda=10, lambda.min=0, family="binomial"))[,10]
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", nlambda=10, lambda.min=0, family="binomial"))[,10]
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
})

test_that("grpreg() reproduces logistic regression", {
  n <- 50
  group <- rep(1:3,1:3)
  ## group <- rep(1:3,rep(3,3))
  p <- length(group)
  X <- matrix(rnorm(n*p),ncol=p)
  y <- runif(n) > .5
  fit.mle <- glm(y~X, family="binomial")
  reg <- coef(fit.mle)
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", family="binomial"))[,100]
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  bridge <- coef(fit <- gBridge(X, y, group, family="binomial"))[,1]
  expect_that(bridge, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", family="binomial"))[,100]
  if (interactive()) {par(mfrow=c(3,1));plot(fit, main=fit$penalty)}
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

test_that("grLasso reproduces lasso", {
  require(glmnet)
  n <- 50
  group <- 1:10
  p <- length(group)
  X <- matrix(rnorm(n*p),ncol=p)
  y <- rnorm(n)
  yy <- runif(n) > .5
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso"))
  lasso <- as.matrix(coef(glmnet(X, y, lambda=fit$lambda)))
  expect_that(grLasso, equals(lasso, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial"))
  lasso <- as.matrix(coef(glmnet(X, yy, family="binomial", lambda=fit$lambda)))
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
  mcp <- coef(ncvreg(X, y, lambda=fit$lambda, penalty="MCP"))
  expect_that(grMCP, equals(mcp, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD"))
  scad <- coef(ncvreg(X, y, lambda=fit$lambda, penalty="SCAD"))
  expect_that(grSCAD, equals(scad, tolerance=.01, check.attributes=FALSE))
})

test_that("logLik() is correct", {
  n <- 50
  group <- rep(1:4,1:4)
  p <- length(group)
  X <- matrix(rnorm(n*p),ncol=p)
  y <- rnorm(n)
  yy <- runif(n) > .5
  fit.mle <- lm(y~X)
  fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0)
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, tol=.001))
  fit <- grpreg(X, y, group, penalty="gMCP", lambda.min=0)
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, tol=.001))
  fit.mle <- glm(yy~X, family="binomial")
  fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="binomial")
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, , tol=.001))
  fit <- grpreg(X, yy, group, penalty="gMCP", lambda.min=0, family="binomial")
  expect_that(logLik(fit)[100], equals(logLik(fit.mle)[1], check.attributes=FALSE, tol=.001))  
  expect_that(AIC(fit)[100], equals(AIC(fit.mle), check.attributes=FALSE, tol=.001))
})

if (interactive()) {
  test_that("cv.grpreg() seems to work", {
    n <- 50
    group <- rep(1:4,1:4)
    p <- length(group)
    X <- matrix(rnorm(n*p),ncol=p)
    b <- rnorm(p)
    b[abs(b) < 1] <- 0
    y <- rnorm(n, mean=X%*%b)
    yy <- runif(n) > .5
    cvfit <- cv.grpreg(X, y, group)
    par(mfrow=c(2,1))
    plot(cvfit)
    cvfit <- cv.grpreg(X, yy, group, family="binomial")
    plot(cvfit)
  })  
}

