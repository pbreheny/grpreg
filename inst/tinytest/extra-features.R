if (interactive()) library(tinytest)

# logLik is correct
n <- 50
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- runif(n) > .5
fit.mle <- lm(y~X)
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0)
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
expect_equivalent(apply(grpreg:::deviance_grpreg(y, predict(fit, X=X, type='response'), family='gaussian'), 2, sum), fit$deviance)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0)
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
expect_equivalent(apply(grpreg:::deviance_grpreg(y, predict(fit, X=X, type='response'), family='gaussian'), 2, sum), fit$deviance, tol=0.0001)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol=.001)

# Binomial
fit.mle <- glm(yy~X, family="binomial")
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="binomial")
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
expect_equivalent(apply(grpreg:::deviance_grpreg(yy, predict(fit, X=X, type='response'), family='binomial'), 2, sum), fit$deviance, tol=0.0001)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit <- grpreg(X, yy, group, penalty="gel", lambda.min=0, family="binomial")
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
expect_equivalent(apply(grpreg:::deviance_grpreg(yy, predict(fit, X=X, type='response'), family='binomial'), 2, sum), fit$deviance, tol=0.0001)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit.mle <- glm(yy~X, family="poisson")
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="poisson")
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
expect_equivalent(apply(grpreg:::deviance_grpreg(yy, predict(fit, X=X, type='response'), family='poisson'), 2, sum), fit$deviance, tol=0.0001)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit <- grpreg(X, yy, group, penalty="gel", lambda.min=0, family="poisson")
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
expect_equivalent(apply(grpreg:::deviance_grpreg(yy, predict(fit, X=X, type='response'), family='poisson'), 2, sum), fit$deviance, tol=0.0001)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol=.001)

# linear predictors are correct
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0)
expect_equivalent(fit$linear.predictor, predict(fit, X))
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0)
expect_equivalent(fit$linear.predictor, predict(fit, X))
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="binomial")
expect_equivalent(fit$linear.predictor, predict(fit, X))
fit <- grpreg(X, yy, group, penalty="gel", lambda.min=0, family="binomial")
expect_equivalent(fit$linear.predictor, predict(fit, X))
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="poisson")
expect_equivalent(fit$linear.predictor, predict(fit, X))
fit <- grpreg(X, yy, group, penalty="gel", lambda.min=0, family="poisson")
expect_equivalent(fit$linear.predictor, predict(fit, X))

# residuals are correct
fit.mle <- lm(y ~ X)
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-12)
expect_equivalent(residuals(fit, lambda=0), residuals(fit.mle))
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, eps=1e-12)
expect_equivalent(residuals(fit, lambda=0), residuals(fit.mle))
fit.mle <- glm(yy ~ X, family="binomial")
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="binomial", eps=1e-12, max.iter=1e6)
expect_equivalent(residuals(fit, lambda=0), residuals(fit.mle), tolerance=1e-4)
fit <- grpreg(X, yy, group, penalty="gel", lambda.min=0, family="binomial", eps=1e-12, max.iter=1e6)
expect_equivalent(residuals(fit, lambda=0), residuals(fit.mle), tolerance=1e-4)
fit.mle <- glm(yy ~ X, family="poisson")
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="poisson", eps=1e-12, max.iter=1e6)
expect_equivalent(residuals(fit, lambda=0), residuals(fit.mle), tolerance=1e-4)
fit <- grpreg(X, yy, group, penalty="gel", lambda.min=0, family="poisson", eps=1e-12, max.iter=1e6)
expect_equivalent(residuals(fit, lambda=0), residuals(fit.mle), tolerance=1e-4)

# grpreg handles user-specified lambda
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- y > 0
fit1 <- grpreg(X, y, group, penalty="grLasso")
fit2 <- grpreg(X, y, group, penalty="grLasso", lambda=fit1$lambda)
expect_equivalent(fit1$beta, fit2$beta)
fit1 <- grpreg(X, y, group, penalty="gel")
fit2 <- grpreg(X, y, group, penalty="gel", lambda=fit1$lambda)
expect_equivalent(fit1$beta, fit2$beta)
fit1 <- grpreg(X, yy, group, penalty="grLasso", family="binomial")
fit2 <- grpreg(X, yy, group, penalty="grLasso", family="binomial", lambda=fit1$lambda)
expect_equivalent(fit1$beta, fit2$beta)
fit1 <- grpreg(X, yy, group, penalty="gel", family="binomial")
fit2 <- grpreg(X, yy, group, penalty="gel", family="binomial", lambda=fit1$lambda)
expect_equivalent(fit1$beta, fit2$beta)

# grpreg named groups
n <- 50
group1 <- rep(0:3,4:1)
group2 <- rep(c("0", "A", "B", "C"), 4:1)
p <- length(group1)
X <- matrix(rnorm(n*p), ncol=p)
X[, group1==2] <- 0
y <- rnorm(n)
yy <- y > 0
fit1 <- grpreg(X, y, group1, penalty="grLasso")
fit2 <- grpreg(X, y, group2, penalty="grLasso")
expect_equivalent(coef(fit1), coef(fit2), tol=0.001)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")

# group.multiplier works
n <- 50
p <- 11
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
group <- rep(0:3, c(1, 2, 3, 5))
gm <- 1:3
plot(fit <- grpreg(X, y, group, penalty="cMCP", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- gBridge(X, y, group, lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grMCP", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grSCAD", lambda.min=0, group.multiplier=gm), main=fit$penalty)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso", group.multiplier=gm)

# group.multiplier works to assign unpenalized groups
## n <- 50
## p <- 11
## X <- matrix(rnorm(n*p),ncol=p)
## y <- rnorm(n)
## group <- rep(0:3, c(1, 2, 3, 5))
## gm <- 0:2
## plot(fit <- grpreg(X, y, group, penalty="cMCP", lambda.min=0, group.multiplier=gm), main=fit$penalty)
## fit <- grpreg(X, y, group, penalty="cMCP", lambda.min=0, group.multiplier=gm)


# dfmax works
n <- 100
group <- rep(1:10, rep(3,10))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- runif(n) > .5
dfmax <- 21
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
expect_true(max(head(nv, length(nv)-1)) <= dfmax)
expect_true(max(nv) > 3)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
expect_true(max(head(nv, length(nv)-1)) <= dfmax)
expect_true(max(nv) > 3)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
expect_true(max(head(nv, length(nv)-1)) <= dfmax)
expect_true(max(nv) > 3)
fit <- grpreg(X, yy, group, penalty="gel", family="binomial", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
expect_true(max(head(nv, length(nv)-1)) <= dfmax)
expect_true(max(nv) > 3)

# gmax works
gmax <- 7
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
expect_true(max(head(ng, length(ng)-1)) <= gmax)
expect_true(max(ng) > 2)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
expect_true(max(head(ng, length(ng)-1)) <= gmax)
expect_true(max(ng) > 2)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
expect_true(max(head(ng, length(ng)-1)) <= gmax)
expect_true(max(ng) > 2)
fit <- grpreg(X, yy, group, penalty="gel", family="binomial", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
expect_true(max(head(ng, length(ng)-1)) <= gmax)
expect_true(max(ng) > 2)
