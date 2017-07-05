.test = "logLik is correct"
n <- 50
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- runif(n) > .5
fit.mle <- lm(y~X)
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0)
check(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
check(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0)
check(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
check(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit.mle <- glm(yy~X, family="binomial")
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="binomial")
check(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
check(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit <- grpreg(X, yy, group, penalty="gMCP", lambda.min=0, family="binomial")
check(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
check(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit.mle <- glm(yy~X, family="poisson")
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="poisson")
check(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
check(AIC(fit)[100], AIC(fit.mle), tol=.001)
fit <- grpreg(X, yy, group, penalty="gMCP", lambda.min=0, family="poisson")
check(logLik(fit)[100], logLik(fit.mle)[1], tol=.001)
check(AIC(fit)[100], AIC(fit.mle), tol=.001)

.test = "grpreg handles user-specified lambda"
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- y > 0
fit1 <- grpreg(X, y, group, penalty="grLasso")
fit2 <- grpreg(X, y, group, penalty="grLasso", lambda=fit1$lambda)
check(fit1$beta, fit2$beta)
fit1 <- grpreg(X, y, group, penalty="gel")
fit2 <- grpreg(X, y, group, penalty="gel", lambda=fit1$lambda)
check(fit1$beta, fit2$beta)
fit1 <- grpreg(X, yy, group, penalty="grLasso", family="binomial")
fit2 <- grpreg(X, yy, group, penalty="grLasso", family="binomial", lambda=fit1$lambda)
check(fit1$beta, fit2$beta)
fit1 <- grpreg(X, yy, group, penalty="gel", family="binomial")
fit2 <- grpreg(X, yy, group, penalty="gel", family="binomial", lambda=fit1$lambda)
check(fit1$beta, fit2$beta)

.test = "grpreg handles constant columns"
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
yy <- y > 0
par(mfrow=c(3,3))
fit <- grpreg(X, y, group, penalty="grLasso"); plot(fit)
fit <- grpreg(X, y, group, penalty="cMCP"); plot(fit)
fit <- gBridge(X, y, group); plot(fit)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial"); plot(fit)
fit <- grpreg(X, yy, group, penalty="cMCP", family="binomial"); plot(fit); fit$beta[,100]
fit <- gBridge(X, yy, group, family="binomial"); plot(fit); fit$beta[,100]
fit <- grpreg(X, yy, group, penalty="grLasso", family="poisson"); plot(fit)
fit <- grpreg(X, yy, group, penalty="cMCP", family="poisson"); plot(fit); fit$beta[,100]
fit <- gBridge(X, yy, group, family="poisson"); plot(fit); fit$beta[,100]
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")
cvfit <- cv.grpreg(X, y, group, penalty="gel")

.test = "grpreg handles groups of non-full rank"
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,7] <- X[,6]
y <- rnorm(n)
yy <- y > 0
par(mfrow=c(2,3))
fit <- grpreg(X, y, group, penalty="grLasso"); plot(fit)
fit <- grpreg(X, y, group, penalty="cMCP"); plot(fit)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial"); plot(fit)
fit <- grpreg(X, yy, group, penalty="cMCP", family="binomial"); plot(fit)
fit <- grpreg(X, yy, group, penalty="grLasso", family="poisson"); plot(fit)
fit <- grpreg(X, yy, group, penalty="cMCP", family="poisson"); plot(fit)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")
cvfit <- cv.grpreg(X, y, group, penalty="gel")

.test = "grpreg out-of-order groups"
n <- 50
group <- rep(0:3,4:1)
ind <- sample(1:length(group))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
yy <- y > 0
fit1 <- grpreg(X, y, group, penalty="grLasso")
fit2 <- grpreg(X[,ind], y, group[ind], penalty="grLasso")
b1 <- coef(fit1)[-1,][ind,]
b2 <- coef(fit2)[-1,]
check(b1, b2, tol=0.01)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")

.test = "grpreg named groups"
n <- 50
group1 <- rep(0:3,4:1)
group2 <- rep(c("0", "A", "B", "C"), 4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
yy <- y > 0
fit1 <- grpreg(X, y, group1, penalty="grLasso")
fit2 <- grpreg(X, y, group2, penalty="grLasso")
check(coef(fit1), coef(fit2), tol=0.001)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")

.test = "grpreg out-of-order groups with constant columns"
n <- 50
group <- rep(c(1,3,0,2),5:2)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
mle <- coef(lm(y~X))
mle[!is.finite(mle)] <- 0
grl <- coef(grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7), lambda=0)
check(mle, grl, tol=0.01)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")

.test = "group.multiplier works"
n <- 50
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
group <- rep(0:3,1:4)
gm <- 1:3
plot(fit <- grpreg(X, y, group, penalty="cMCP", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- gBridge(X, y, group, lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grMCP", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grSCAD", lambda.min=0, group.multiplier=gm), main=fit$penalty)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso", group.multiplier=gm)

.test = "dfmax works"
n <- 100
group <- rep(1:10, rep(3,10))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- runif(n) > .5
dfmax <- 21
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
check(max(head(nv, length(nv)-1)) <= dfmax)
check(max(nv) > 3)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
check(max(head(nv, length(nv)-1)) <= dfmax)
check(max(nv) > 3)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
check(max(head(nv, length(nv)-1)) <= dfmax)
check(max(nv) > 3)
fit <- grpreg(X, yy, group, penalty="gel", family="binomial", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
check(max(head(nv, length(nv)-1)) <= dfmax)
check(max(nv) > 3)

.test = "gmax works"
gmax <- 7
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
check(max(head(ng, length(ng)-1)) <= gmax)
check(max(ng) > 2)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
check(max(head(ng, length(ng)-1)) <= gmax)
check(max(ng) > 2)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
check(max(head(ng, length(ng)-1)) <= gmax)
check(max(ng) > 2)
fit <- grpreg(X, yy, group, penalty="gel", family="binomial", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
check(max(head(ng, length(ng)-1)) <= gmax)
check(max(ng) > 2)

