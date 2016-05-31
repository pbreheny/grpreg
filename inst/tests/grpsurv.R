set.seed(1)

.test = "grpsurv works"
y <- survival::Surv(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*6), 50, 6)
g <- rep(1:3, each=2)
fit <- grpsurv(X, y, g, lambda.min=0)

.test = "grpsurv equals MLE when lam=0"
fit.mle <- survival::coxph(y~X)
check(coef(fit)[,100], coef(fit.mle), tol=0.01)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="grSCAD")
check(coef(fit)[,100], coef(fit.mle), tol=0.01)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="grMCP")
check(coef(fit)[,100], coef(fit.mle), tol=0.01)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="gel")
check(coef(fit)[,100], coef(fit.mle), tol=0.01)
fit <- grpsurv(X, y, g, lambda.min=0, penalty="cMCP")
check(coef(fit)[,100], coef(fit.mle), tol=0.01)

.test = "grpsurv returns correct logLik"
check(logLik(fit)[100], logLik(fit.mle)[1], tol=0.1)
check(AIC(fit)[100], AIC(fit.mle), tol=0.1)

.test = "predict works for grpsurv"
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

.test = "penalty factor works for grpsurv"
fit <- grpsurv(X, y, g, group.multiplier=c(1:3))

.test = "coersion works for grpsurv"
fit <- grpsurv(as.data.frame(X), y)

.test = "loss works for grpsurv"
eta <- predict(fit, X, 'link', lambda=0.1)
grpreg:::loss.grpsurv(y, eta)

.test = "cross-validation works for grpsurv"
cvfit <- cv.grpsurv(X, y, g, lambda.min=0)
plot(cvfit)
summary(cvfit)
