source("../../../.debug-grpreg.R")

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
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", nlambda=10, max.iter=10, lambda.min=0))[,10]
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", nlambda=10, max.iter=10, lambda.min=0))[,10]
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", nlambda=10, max.iter=10, lambda.min=0))[,10]
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", nlambda=10, max.iter=10, lambda.min=0))[,10]
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
})

test_that("grpreg() reproduces linear regression", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  y <- rnorm(n)
  group <- rep(1:4,1:4)
  reg <- lm(y~X)$coef
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", lambda.min=0))[,100]
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0))[,100]
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", lambda.min=0))[,100]
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", lambda.min=0))[,100]
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
})

test_that("grpreg() reproduces simple logistic regression", {
  n <- 20
  p <- 1
  X <- matrix(rnorm(n*p), ncol=p)
  y <- runif(n) > .5
  group <- 1
  reg <- glm(y~X, family="binomial")$coef
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", nlambda=10, max.iter=10, lambda.min=0, family="binomial"))[,10]
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", nlambda=10, max.iter=10, lambda.min=0, family="binomial"))[,10]
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", nlambda=10, max.iter=10, lambda.min=0, family="binomial"))[,10]
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", nlambda=10, max.iter=10, lambda.min=0, family="binomial"))[,10]
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
})

test_that("grpreg() reproduces logistic regression", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  group <- rep(1:4,1:4)
  y <- runif(n) > .5
  reg <- glm(y~X, family="binomial")$coef
  gMCP <- coef(fit <- grpreg(X, y, group, penalty="gMCP", family="binomial"))[,100]
  expect_that(gMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", family="binomial"))[,100]
  expect_that(grLasso, equals(reg, tolerance=.01, check.attributes=FALSE))
  grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", family="binomial"))[,100]
  expect_that(grMCP, equals(reg, tolerance=.01, check.attributes=FALSE))
  grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", , family="binomial"))[,100]
  expect_that(grSCAD, equals(reg, tolerance=.01, check.attributes=FALSE))
})
