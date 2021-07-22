if (basename(getwd()) == 'tinytest') source('_median-survfit.R') else {library(tinytest); source ('inst/tinytest/_median-survfit.R')}
suppressPackageStartupMessages(library(survival))

# Test that grpsurv works when x has 1 column
y <- Surv(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*1), 50, 1)
g <- 1
fit <- grpsurv(X, y, g, lambda.min=0)

# Test that grpsurv works
y <- Surv(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*6), 50, 6)
g <- rep(1:3, each=2)
fit <- grpsurv(X, y, g, lambda.min=0)

# $ grpsurv equals MLE when lam=0
fit.mle <- coxph(y~X)
expect_equivalent(coef(fit)[,100], coef(fit.mle), tol=0.01)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="grSCAD")
expect_equivalent(coef(fit)[,100], coef(fit.mle), tol=0.01)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="grMCP")
expect_equivalent(coef(fit)[,100], coef(fit.mle), tol=0.01)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="gel")
expect_equivalent(coef(fit)[,100], coef(fit.mle), tol=0.01)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="cMCP")
expect_equivalent(coef(fit)[,100], coef(fit.mle), tol=0.01)

# grpsurv returns correct logLik
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol=0.1)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol=0.1)

# grpsurv returns correct linear predictors
fit <- grpsurv(X, y, g, lambda.min=0, penalty="grLasso", eps=1e-12)
expect_equivalent(fit$linear.predictors[, 100], fit.mle$linear.predictors[order(y)])
fit <- grpsurv(X, y, g, lambda.min=0, penalty="gel", eps=1e-12)
expect_equivalent(fit$linear.predictors[, 100], fit.mle$linear.predictors[order(y)])

# residuals are correct (slightly different at final observation because baseline hazard estimated differently)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="grLasso", eps=1e-12)
r1 <- residuals(fit, lambda=0)
r2 <- residuals(fit.mle, type='deviance')
plot(r1, r2)
abline(0,1)
expect_equivalent(residuals(fit, lambda=0), residuals(fit.mle, type='deviance'), tolerance=0.1)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="gel", eps=1e-12)
expect_equivalent(residuals(fit, lambda=0), residuals(fit.mle, type='deviance'), tolerance=0.1)

# predict works for grpsurv
head(predict(fit, X, 'vars'))
head(predict(fit, X, 'nvars'))
head(predict(fit, X, 'groups'))
head(predict(fit, X, 'ngroups'))
head(predict(fit, X, 'link'))
head(predict(fit, X, 'response'))
head(predict(fit, X, 'coef'))
head(predict(fit, X, 'median'))
S <- predict(fit, X, 'survival', lambda=0.1)
plot(S)
S <- predict(fit, X[1,], 'survival', lambda=0.1)
plot(S)
H <- predict(fit, X, 'hazard', lambda=0.1)
plot(H)
H <- predict(fit, X[1,], 'hazard', lambda=0.1)
plot(H)

# Survival curves vs survival package
DF <- as.data.frame(X)
fit.mle <- coxph(y~., DF)
S <- predict(fit, X[1,], 'survival', lambda=0)
plot(S)
lines(survfit(fit.mle, DF[1,], conf.int=FALSE))
H <- predict(fit, X[1,], 'hazard', lambda=0)
plot(H)
lines(survfit(fit.mle, DF[1,]), conf.int=FALSE, cumhaz=TRUE)

# Cumulative hazard
h1 <- H(fit$time)
h2 <- survfit(fit.mle, DF[1,])$cumhaz
plot(h1, h2)
abline(0, 1)

# penalty factor works for grpsurv
fit <- grpsurv(X, y, g, group.multiplier=c(1:3))

# coersion works for grpsurv
fit <- grpsurv(as.data.frame(X), y)

# loss works for grpsurv
eta <- predict(fit, X, 'link', lambda=0.1)
grpreg:::loss.grpsurv(y, eta)
grpreg:::loss.grpsurv(y, eta, total=FALSE)

# cross-validation works for grpsurv
cvfit <- cv.grpsurv(X, y, g, lambda.min=0)
plot(cvfit)
summary(cvfit)
