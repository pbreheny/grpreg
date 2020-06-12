source('_median-survfit.R')

# Test that grpsurv works
y <- survival::Surv(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*6), 50, 6)
g <- rep(1:3, each=2)
fit <- grpsurv(X, y, g, lambda.min=0)

# $ grpsurv equals MLE when lam=0
fit.mle <- survival::coxph(y~X)
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

# predict works for grpsurv
head(predict(fit, X, 'vars'))
head(predict(fit, X, 'nvars'))
head(predict(fit, X, 'groups'))
head(predict(fit, X, 'ngroups'))
head(predict(fit, X, 'link'))
head(predict(fit, X, 'response'))
head(predict(fit, X, 'coef'))
S <- predict(fit, X, 'survival', lambda=0.1)
plot(S)
S <- predict(fit, X[1,], 'survival', lambda=0.1)
plot(S)
head(predict(fit, X, 'median'))

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
