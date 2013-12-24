## Tests concerning seemingly unrelated regressions/multitask learning
test_that("multitask learning works", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  b <- c(-3, 3, -1, 1, rep(0, 6))
  Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
             5+rnorm(n, mean=X%*%b, sd=1),
             10+rnorm(n, mean=X%*%b, sd=1))
  colnames(X) <- LETTERS[1:10]
  
  par(mfcol=c(3,2))
  fit <- grpreg(X, Y, penalty="grLasso"); plot(fit)
  fit <- grpreg(X, Y, penalty="cMCP"); plot(fit)
  fit <- gBridge(X, Y); plot(fit)
  
  n <- 200
  X <- matrix(rnorm(n*p),ncol=p)
  Y <- cbind(rnorm(n, mean=X%*%b/10, sd=1),
             rnorm(n, mean=X%*%b/10, sd=1),
             rnorm(n, mean=X%*%b/10, sd=1)) > 0
  fit <- grpreg(X, Y, penalty="grLasso", family="binomial", lambda.min=0.1); plot(fit)
  fit <- grpreg(X, Y, penalty="cMCP", family="binomial", lambda.min=0.2); plot(fit)
  fit <- gBridge(X, Y, family="binomial", lambda.min=0.2); plot(fit)
})

test_that("coef/predict works for multitask learning", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  b <- c(-3, 3, -1, 1, rep(0, 6))
  Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
             5+rnorm(n, mean=X%*%b, sd=1),
             10+rnorm(n, mean=X%*%b, sd=1))
  colnames(X) <- LETTERS[1:10]
  
  par(mfcol=c(3,2))
  fit <- grpreg(X, Y)
  coef(fit, which=1:2)
  coef(fit, lambda=1)
  predict(fit, lambda=1, type="vars")
  predict(fit, which=30, type="vars")
  predict(fit, lambda=1, type="groups")
  predict(fit, which=c(30,60), type="groups")
  predict(fit, lambda=1, type="norm")
  predict(fit, which=c(30,60), type="norm")
  predict(fit, X, lambda=1)
  predict(fit, X, which=c(30,60))
  
  n <- 200
  X <- matrix(rnorm(n*p),ncol=p)
  Y <- cbind(rnorm(n, mean=X%*%b/3, sd=1),
             rnorm(n, mean=X%*%b/3, sd=1),
             rnorm(n, mean=X%*%b/3, sd=1)) > 0
  fit <- grpreg(X, Y, family="binomial")
  predict(fit, lambda=0.1, type="vars")
  predict(fit, lambda=0.1, type="groups")
  predict(fit, lambda=0.1, type="norm")
  predict(fit, X, lambda=0.1)
  predict(fit, X, lambda=0.1, type="response")
  predict(fit, X, lambda=0.1, type="class")
})

test_that("multitask learning reproduces linear regression", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  Y <- matrix(rnorm(n*3),ncol=3)
  fit.mle <- lm(Y~X)
  reg <- coef(fit.mle)
  
  cMCP <- coef(fit <- grpreg(X, Y, penalty="cMCP", lambda.min=0), which=100)
  expect_that(t(cMCP), equals(reg, tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=100)
  expect_that(p, equals(predict(fit.mle), tolerance=.01, check.attributes=FALSE))
  
  bridge <- coef(fit <- gBridge(X, Y, lambda.min=0), which=1)
  expect_that(t(bridge), equals(reg, tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=1)
  expect_that(p, equals(predict(fit.mle), tolerance=.01, check.attributes=FALSE))
  
  grLasso <- coef(fit <- grpreg(X, Y, penalty="grLasso", lambda.min=0), which=100)
  expect_that(t(grLasso), equals(reg, tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=100)
  expect_that(p, equals(predict(fit.mle), tolerance=.01, check.attributes=FALSE))
  
  grMCP <- coef(fit <- grpreg(X, Y, penalty="grMCP", lambda.min=0), which=100)
  expect_that(t(grMCP), equals(reg, tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=100)
  expect_that(p, equals(predict(fit.mle), tolerance=.01, check.attributes=FALSE))
  
  grSCAD <- coef(fit <- grpreg(X, Y, penalty="grSCAD", lambda.min=0), which=100)
  expect_that(t(grSCAD), equals(reg, tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=100)
  expect_that(p, equals(predict(fit.mle), tolerance=.01, check.attributes=FALSE))
})

test_that("multitask learning reproduces logistic regression", {
  n <- 200
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  Y <- matrix(rnorm(n*3),ncol=3)>0
  fit.mle <- glm(Y[,3]~X, family=binomial)
  mle <- coef(fit.mle)
  
  beta <- coef(fit <- grpreg(X, Y, lambda.min=0, family="binomial"), which=100)[3,]
  expect_that(beta, equals(mle, tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=100)[,3]
  expect_that(p, equals(predict(fit.mle), tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=100, type="response")[,3]
  expect_that(p, equals(predict(fit.mle, type="response"), tolerance=.01, check.attributes=FALSE))
  
  bridge <- coef(fit <- gBridge(X, Y, family="binomial", lambda.min=0), which=1)[3,]
  expect_that(bridge, equals(mle, tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=1)[,3]
  expect_that(p, equals(predict(fit.mle), tolerance=.01, check.attributes=FALSE))
  p <- predict(fit, X, which=1, type="response")[,3]
  expect_that(p, equals(predict(fit.mle, type="response"), tolerance=.01, check.attributes=FALSE))
})

test_that("cross-validation for multitask learning works", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  colnames(X) <- LETTERS[1:10]
  b <- c(-3, 3, -1, 1, rep(0, 6))
  Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
             5+rnorm(n, mean=X%*%b, sd=1),
             10+rnorm(n, mean=X%*%b, sd=1))
  
  par(mfcol=c(2,2))
  cvfit <- cv.grpreg(X, Y); plot(cvfit, type="all")
  cvfit <- cv.grpreg(X, Y, penalty="cMCP"); plot(cvfit, type="all")
  
  b <- rep(0,10)
  Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
             5+rnorm(n, mean=X%*%b, sd=1),
             10+rnorm(n, mean=X%*%b, sd=1))
  cvfit <- cv.grpreg(X, Y); plot(cvfit, type="all")
  cvfit <- cv.grpreg(X, Y, penalty="cMCP"); plot(cvfit, type="all")
  
  n <- 200
  X <- matrix(rnorm(n*p),ncol=p)
  b <- c(-3, 3, -1, 1, rep(0, 6))
  Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
             rnorm(n, mean=X%*%b, sd=1),
             rnorm(n, mean=X%*%b, sd=1)) > 0
  cvfit <- cv.grpreg(X, Y, family="binomial"); plot(cvfit, type="all")
  cvfit <- cv.grpreg(X, Y, penalty="cMCP", family="binomial"); plot(cvfit, type="all")

  b <- rep(0,10)
  Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
             rnorm(n, mean=X%*%b, sd=1),
             rnorm(n, mean=X%*%b, sd=1)) > 0
  cvfit <- cv.grpreg(X, Y, family="binomial"); plot(cvfit, type="all")
  cvfit <- cv.grpreg(X, Y, penalty="cMCP", family="binomial"); plot(cvfit, type="all")
})
