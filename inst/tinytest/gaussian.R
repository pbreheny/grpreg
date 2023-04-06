# Simple linear regression ------------------------------------------------

n <- 5
p <- 1
X <- matrix(rnorm(n*p), ncol=p)
y <- rnorm(n)
group <- 1
reg <- lm(y~X)$coef
nlam=100
par(mfcol=c(3,2))
gel <- coef(fit <- grpreg(X, y, group, penalty="gel", nlambda=nlam, lambda.min=0))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(gel, reg, tolerance=1e-7)
cMCP <- coef(fit <- grpreg(X, y, group, penalty="cMCP", nlambda=nlam, lambda.min=0))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(cMCP, reg, tolerance=1e-7)
bridge <- coef(fit <- gBridge(X, y, group, nlambda=nlam, lambda.min=0))[,1]
plot(fit, main=fit$penalty)
expect_equivalent(bridge, reg, tolerance=1e-7)
grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", nlambda=nlam, lambda.min=0))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(grLasso, reg, tolerance=1e-7)
grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", nlambda=nlam, lambda.min=0))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(grMCP, reg, tolerance=1e-7)
grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", nlambda=nlam, lambda.min=0))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(grSCAD, reg, tolerance=1e-7)


# Linear regression -------------------------------------------------------

n <- 50
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
group <- rep(0:3,4:1)
fit.mle <- lm(y~X)
reg <- coef(fit.mle)
nlam=100
par(mfcol=c(3,2))
gel <- coef(fit <- grpreg(X, y, group, penalty="gel", nlambda=nlam, lambda.min=0, eps=1e-10))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(gel, reg, tolerance=1e-7)
cMCP <- coef(fit <- grpreg(X, y, group, penalty="cMCP", lambda.min=0, eps=1e-10))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(cMCP, reg, tolerance=1e-7)
bridge <- coef(fit <- gBridge(X, y, group, lambda.min=0, eps=1e-10))[,1]
plot(fit, main=fit$penalty)
expect_equivalent(bridge, reg, tolerance=1e-7)
grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-10))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grLasso, reg, tolerance=1e-7)
grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", lambda.min=0, eps=1e-10))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grMCP, reg, tolerance=1e-7)
grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", lambda.min=0, eps=1e-10))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grSCAD, reg, tolerance=1e-7)
expect_equivalent(predict(fit, X)[,100], predict(fit.mle), tolerance=1e-7)
