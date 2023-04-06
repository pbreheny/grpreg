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
