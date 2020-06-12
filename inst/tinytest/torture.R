# constant columns
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,5] <- 0
y <- rnorm(n)
b.lm <- coef(lm(y~X)); b.lm[is.na(b.lm)] <- 0
b <- coef(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b.lm, tol=0.0001)
b <- coef(fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b.lm, tol=0.0001)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")
cvfit <- cv.grpreg(X, y, group, penalty="gel")

# constant groups
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
b.lm <- coef(lm(y~X)); b.lm[is.na(b.lm)] <- 0
b <- coef(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b.lm, tol=0.0001)
b <- coef(fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b.lm, tol=0.0001)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")
cvfit <- cv.grpreg(X, y, group, penalty="gel")

# constant groups w/ multiplier specified
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
b.lm <- coef(lm(y~X)); b.lm[is.na(b.lm)] <- 0
b <- coef(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7, group.multiplier=1:3), 0)
expect_equivalent(b, b.lm, tol=0.0001)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")
cvfit <- cv.grpreg(X, y, group, penalty="gel")

# grpreg handles groups of non-full rank
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,7] <- X[,6]
y <- rnorm(n)
b0 <- coef(lm(y~X)); b0[7:8] <- b0[7]/2
b <- coef(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7), 0)
expect_equivalent(b, b0, tol=0.0001)
cvfit <- cv.grpreg(X, y, group, penalty="grLasso")
cvfit <- cv.grpreg(X, y, group, penalty="gel")

# out-of-order groups #1
n <- 50
group <- rep(1:2, each=2)
perm <- sample(1:length(group))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
fit1 <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7)
fit2 <- grpreg(X[,perm], y, group[perm], penalty="grLasso", lambda.min=0, eps=1e-7)
b1 <- coef(fit1, 0)[-1][perm]
b2 <- coef(fit2, 0)[-1]
expect_equivalent(b1, b2, tol=0.01)
cvfit <- cv.grpreg(X[,perm], y, group[perm], penalty="grLasso")

# out-of-order groups #2
n <- 50
group <- rep(0:3,4:1)
perm <- sample(1:length(group))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
fit1 <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, eps=1e-7)
fit2 <- grpreg(X[,perm], y, group[perm], penalty="grLasso", lambda.min=0, eps=1e-7)
b1 <- coef(fit1, 0)[-1][perm]
b2 <- coef(fit2, 0)[-1]
expect_equivalent(b1, b2, tol=0.01)
cvfit <- cv.grpreg(X[,perm], y, group[perm], penalty="grLasso")

# groups order + rank + constant col + constant grp
n <- 50
group <- rep(0:4, c(2, 2:5))
perm <- sample(1:length(group))
p <- length(group)
X <- matrix(rnorm(n*p), ncol=p)
X[,7] <- X[,6]    # Group 2 not full rank
X[,group==3] <- 0 # Group 3 constant
X[,15] <- 0       # Group 4 contains a zero column
y <- rnorm(n)
fit1 <- grpreg(X, y, group, penalty="grLasso", lambda.min=0)
fit2 <- grpreg(X[,perm], y, group[perm], penalty="grLasso", lambda.min=0)
b1 <- coef(fit1, 0)[-1]
b2 <- coef(fit2, 0)[-1]
expect_equivalent(b1[perm], b2, tol=0.01) # Checking perm ordering
nz <- which(apply(X, 2, sd)!=0)
fit3 <- grpreg(X[,nz], y, group[nz], penalty="grLasso", lambda.min=0)
b3 <- coef(fit3, 0)[-1]
expect_equivalent(b1[nz], b3, tol=0.01)  # Checking dropped group/var
expect_equivalent(coef(fit1)["V6",], coef(fit1)["V7",], tol=0.01)  # Checking rank handled correctly
cvfit <- cv.grpreg(X[,perm], y, group[perm], penalty="grLasso", lambda.min=0)
plot(cvfit)
summary(cvfit)
plot(cvfit$fit)
