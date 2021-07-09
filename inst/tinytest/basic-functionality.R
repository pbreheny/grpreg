# grpreg() reproduces simple linear regression
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

# grpreg reproduces linear regression
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

# grpreg reproduces simple logistic regression
n <- 20
p <- 1
X <- matrix(rnorm(n*p), ncol=p)
y <- runif(n) > .5
group <- 1
reg <- glm(y~X, family="binomial")$coef
nlam <- 100
par(mfcol=c(3,2))
gel <- coef(fit <- grpreg(X, y, group, penalty="gel", nlambda=nlam, lambda.min=0, family="binomial", eps=1e-10))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(gel, reg, tolerance=1e-7)
cMCP <- coef(fit <- grpreg(X, y, group, penalty="cMCP", nlambda=nlam, lambda.min=0, family="binomial", gamma=9, eps=1e-10))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(cMCP, reg, tolerance=1e-7)
bridge <- coef(fit <- gBridge(X, y, group, lambda.min=0, nlambda=nlam, family="binomial", eps=1e-10))[,1]
plot(fit, main=fit$penalty)
expect_equivalent(bridge, reg, tolerance=1e-7)
grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", nlambda=nlam, lambda.min=0, family="binomial", eps=1e-10))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(grLasso, reg, tolerance=1e-7)
grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", nlambda=nlam, lambda.min=0, family="binomial", eps=1e-10))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(grMCP, reg, tolerance=1e-7)
grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", nlambda=nlam, lambda.min=0, family="binomial", eps=1e-10))[,nlam]
plot(fit, main=fit$penalty)
expect_equivalent(grSCAD, reg, tolerance=1e-7)

# grpreg() reproduces logistic regression
n <- 100
group <- rep(0:3,1:4)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- runif(n) > .5
fit.mle <- glm(y~X, family="binomial")
reg <- coef(fit.mle)
par(mfcol=c(3,2))
gel <- coef(fit <- grpreg(X, y, group, penalty="gel", family="binomial", eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(gel, reg, tol=1e-6)
cMCP <- coef(fit <- grpreg(X, y, group, penalty="cMCP", family="binomial", eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(cMCP, reg, tol=1e-6)
bridge <- coef(fit <- gBridge(X, y, group, family="binomial", eps=1e-10, lambda.min=0))[,1]
plot(fit, main=fit$penalty)
expect_equivalent(bridge, reg, tol=1e-6)
grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", family="binomial", eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grLasso, reg, tol=1e-6)
grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", family="binomial", gamma=2, eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grMCP, reg, tol=1e-6)
grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", family="binomial", gamma=2.1, eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grSCAD, reg, tol=1e-6)
expect_equivalent(predict(fit, X)[,100], predict(fit.mle), tol=1e-6)
expect_equivalent(predict(fit, X, type="response")[,100], predict(fit.mle, type="response"), tol=1e-6)

# grpreg() reproduces poisson regression
n <- 50
group <- rep(0:3,1:4)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- sample(1:5, n, replace=TRUE)
fit.mle <- glm(y~X, family="poisson")
reg <- coef(fit.mle)
par(mfcol=c(3,2))
gel <- coef(fit <- grpreg(X, y, group, penalty="gel", family="poisson", eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(gel, reg, tolerance=1e-7)
cMCP <- coef(fit <- grpreg(X, y, group, penalty="cMCP", family="poisson", eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(cMCP, reg, tolerance=1e-7)
bridge <- coef(fit <- gBridge(X, y, group, family="poisson", eps=1e-10, lambda.min=0))[,1]
plot(fit, main=fit$penalty)
expect_equivalent(bridge, reg, tolerance=1e-7)
grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso", family="poisson", eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grLasso, reg, tolerance=1e-7)
grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", family="poisson", gamma=2, eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grMCP, reg, tolerance=1e-7)
grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", family="poisson", gamma=2.1, eps=1e-10, lambda.min=0))[,100]
plot(fit, main=fit$penalty)
expect_equivalent(grSCAD, reg, tolerance=1e-7)
expect_equivalent(predict(fit, X)[,100], predict(fit.mle), tolerance=1e-7)
expect_equivalent(predict(fit, X, type="response")[,100], predict(fit.mle, type="response"), tolerance=1e-7)
