require(survival)
source('~/dev/.grpreg.setup.R')
#set.seed(1)
equal <- function(x, y, tol=0.001) {all.equal(x, y, tol=tol, check.attributes=FALSE)}

# Works
y <- Surv(rexp(50), sample(rep(0:1, c(10,40))))
X <- matrix(rnorm(50*6), 50, 6)
g <- rep(1:3, each=2)
fit <- grpsurv(X, y, g, lambda.min=0)

# Equals MLE when lam=0
fit.mle <- coxph(y~X)
stopifnot(equal(coef(fit)[,100], coef(fit.mle), tol=0.01))

# logLik
stopifnot(equal(logLik(fit)[100], logLik(fit.mle)[1]))
stopifnot(equal(AIC(fit)[100], AIC(fit.mle)))

# Other penalties
fit <- grpsurv(X, y, g, lambda.min=0, penalty="grSCAD")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))
fit <- grpsurv(X, y, g, lambda.min=0, penalty="grMCP")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))
fit <- grpsurv(X, y, g, lambda.min=0, penalty="gel")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))
fit <- grpsurv(X, y, g, lambda.min=0, penalty="cMCP")
stopifnot(equal(coef(fit)[,100], coef(fit.mle)))

# Predict
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')
p <- predict(fit, X, 'groups')
p <- predict(fit, X, 'ngroups')
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
S <- predict(fit, X, 'survival', lambda=0.1)
plot(S)
S <- predict(fit, X[1,], 'survival', lambda=0.1)
plot(S)
p <- predict(fit, X, 'median')

# Penalty factor
fit <- grpsurv(X, y, g)
plot(fit, norm=TRUE, legend="topleft", bty="n")
fit <- grpsurv(X, y, g, group.multiplier=c(1:3))
plot(fit, norm=TRUE, legend="topleft", bty="n")

# Coersion
fit <- grpsurv(as.data.frame(X), y)

# Loss
eta <- predict(fit, X, 'link', lambda=0.1)
loss.grpsurv(y, eta)

# Cross-validation
cvfit <- cv.grpsurv(X, y, g, lambda.min=0)
plot(cvfit)
summary(cvfit)

#############################################################
.test = "cox survival predictions agree with Kaplan-Meier" ##
#############################################################
n <- 50
p <- 5
X <- matrix(rnorm(n*p), ncol=p)
b <- c(2, -2, rep(0, p-2))
g <- c(1,3,3,2,2)
y <- Surv(rexp(n, exp(X%*%b)), rbinom(n, 1, 0.75))

fit <- grpsurv(X, y, g)
fit <- grpsurv(X, y, g, penalty="gel")
S <- predict(fit, X[1,], which=1, type='survival')
km <- survfit(y~1)
par(op)
plot(km, conf.int=FALSE, mark.time=FALSE, xlim=c(0,10), lwd=10, col="gray")
lines(fit$time, S(fit$time), type="s", col="slateblue", lwd=2)
median(km)
predict(fit, X[1,], which=1, type='median')
