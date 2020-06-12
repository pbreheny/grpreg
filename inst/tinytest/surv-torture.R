suppressPackageStartupMessages(library(survival))

# constant columns
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,5] <- 0
y <- Surv(rexp(n), rep(0:1, c(10, n-10)))
b.mle <- coef(coxph(y~X)); b.mle[is.na(b.mle)] <- 0
b <- coef(fit <- grpsurv(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b.mle, tol=0.0001)
b <- coef(fit <- grpsurv(X, y, group, penalty="gel", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b.mle, tol=0.0001)
cvfit <- cv.grpsurv(X, y, group, penalty="grLasso")
cvfit <- cv.grpsurv(X, y, group, penalty="gel")

# constant groups
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- Surv(rexp(n), rep(0:1, c(10, n-10)))
b.mle <- coef(coxph(y~X)); b.mle[is.na(b.mle)] <- 0
b <- coef(fit <- grpsurv(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b.mle, tol=0.0001)
b <- coef(fit <- grpsurv(X, y, group, penalty="gel", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b.mle, tol=0.0001)
cvfit <- cv.grpsurv(X, y, group, penalty="grLasso")
cvfit <- cv.grpsurv(X, y, group, penalty="gel")

# constant groups w/ multiplier specified
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- Surv(rexp(n), rep(0:1, c(10, n-10)))
b.mle <- coef(coxph(y~X)); b.mle[is.na(b.mle)] <- 0
b <- coef(fit <- grpsurv(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7, group.multiplier=1:3), 0)
expect_equivalent(b, b.mle, tol=0.0001)
cvfit <- cv.grpsurv(X, y, group, penalty="grLasso")
cvfit <- cv.grpsurv(X, y, group, penalty="gel")

# grpsurv handles groups of non-full rank
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,7] <- X[,6]
y <- Surv(rexp(n), rep(0:1, c(10, n-10)))
b0 <- coef(coxph(y~X)); b0[6:7] <- b0[6]/2
b <- coef(fit <- grpsurv(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b0, tol=0.0001)
cvfit <- cv.grpsurv(X, y, group, penalty="grLasso")
cvfit <- cv.grpsurv(X, y, group, penalty="gel")

# out-of-order groups #1
n <- 50
group <- rep(1:2, each=2)
perm <- sample(1:length(group))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- Surv(rexp(n), rep(0:1, c(10, n-10)))
fit1 <- grpsurv(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7)
fit2 <- grpsurv(X[,perm], y, group[perm], penalty="grLasso", lambda.min=0, eps=1e-7)
b1 <- coef(fit1, 0)[perm]
b2 <- coef(fit2, 0)
expect_equivalent(b1, b2, tol=0.01)
cvfit <- cv.grpsurv(X[,perm], y, group[perm], penalty="grLasso")

# out-of-order groups #2
n <- 50
group <- rep(0:3,4:1)
perm <- sample(1:length(group))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- Surv(rexp(n), rep(0:1, c(10, n-10)))
fit1 <- grpsurv(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7)
fit2 <- grpsurv(X[,perm], y, group[perm], penalty="grLasso", lambda.min=0, eps=1e-7)
b1 <- coef(fit1, 0)[perm]
b2 <- coef(fit2, 0)
expect_equivalent(b1, b2, tol=0.01)
cvfit <- cv.grpsurv(X[,perm], y, group[perm], penalty="grLasso")

# groups order + rank + constant col + constant grp
n <- 50
group <- rep(0:4, c(2, 2:5))
perm <- sample(1:length(group))
p <- length(group)
X <- matrix(rnorm(n*p), ncol=p)
X[,7] <- X[,6]    # Group 2 not full rank
X[,group==3] <- 0 # Group 3 constant
X[,15] <- 0       # Group 4 contains a zero column
y <- Surv(rexp(n), rep(0:1, c(10, n-10)))
fit1 <- grpsurv(X, y, group, penalty="grLasso", lambda.min=0)
fit2 <- grpsurv(X[,perm], y, group[perm], penalty="grLasso", lambda.min=0)
b1 <- coef(fit1, 0)
b2 <- coef(fit2, 0)
expect_equivalent(b1[perm], b2, tol=0.01) # Checking perm ordering
nz <- which(apply(X, 2, sd)!=0)
fit3 <- grpsurv(X[,nz], y, group[nz], penalty="grLasso", lambda.min=0)
b3 <- coef(fit3, 0)
expect_equivalent(b1[nz], b3, tol=0.01)  # Checking dropped group/var
expect_equivalent(coef(fit1)["V6",], coef(fit1)["V7",], tol=0.01)  # Checking rank handled correctly
cvfit <- cv.grpsurv(X[,perm], y, group[perm], penalty="grLasso", lambda.min=0)
plot(cvfit)
summary(cvfit)
plot(cvfit$fit)
