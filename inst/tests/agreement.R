suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(ncvreg))

# gel reproduces lasso
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
expect_equivalent(gel, lasso, tolerance=.01)
gel <- coef(fit <- grpreg(X, yy, group, penalty="gel", family="binomial", tau=0))
plot(fit, log=TRUE)
lasso <- as.matrix(coef(fit <- glmnet(X, yy, family="binomial", lambda=fit$lambda)))
plot(fit, "lambda")
expect_equivalent(gel, lasso, tolerance=.01)
gel <- coef(fit <- grpreg(X, yy, group, penalty="gel", family="poisson", tau=0))
plot(fit, log=TRUE)
lasso <- as.matrix(coef(fit <- glmnet(X, yy, family="poisson", lambda=fit$lambda)))
plot(fit, "lambda")
expect_equivalent(gel, lasso, tolerance=.01)

# grLasso reproduces lasso
n <- 50
group <- 1:10
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- runif(n) > .5

grLasso <- coef(fit <- grpreg(X, y, group, penalty="grLasso"))
par(mfrow=c(3,2)); plot(fit, log=TRUE)
fit1 <- fit
lasso <- as.matrix(coef(fit <- glmnet(X, y, lambda=fit$lambda)))
plot(fit, "lambda")
fit2 <- fit
expect_equivalent(grLasso, lasso, tolerance=.01)

grLasso <- coef(fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial"))
plot(fit, log=TRUE)
lasso <- as.matrix(coef(fit <- glmnet(X, yy, family="binomial", lambda=fit$lambda)))
plot(fit, "lambda")
expect_equivalent(grLasso, lasso, tolerance=.01)
grLasso <- coef(fit <- grpreg(X, yy, group, penalty="grLasso", family="poisson"))
plot(fit, log=TRUE)
lasso <- as.matrix(coef(fit <- glmnet(X, yy, family="poisson", lambda=fit$lambda)))
plot(fit, "lambda")
expect_equivalent(grLasso, lasso, tolerance=.01)


# grMCP and grSCAD reproduce MCP and SCAD lasso
n <- 50
group <- 1:10
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
grMCP <- coef(fit <- grpreg(X, y, group, penalty="grMCP", gamma=3))
par(mfrow=c(2,2)); plot(fit)
mcp <- coef(fit <- ncvreg(X, y, lambda=fit$lambda, penalty="MCP", gamma=3))
plot(fit)
expect_equivalent(grMCP, mcp, tolerance=.01)
grSCAD <- coef(fit <- grpreg(X, y, group, penalty="grSCAD", gamma=4))
plot(fit)
scad <- coef(fit <- ncvreg(X, y, lambda=fit$lambda, penalty="SCAD", gamma=4))
plot(fit)
expect_equivalent(grSCAD, scad, tolerance=.01)
