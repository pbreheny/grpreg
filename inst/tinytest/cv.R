if (interactive()) library(tinytest)

# Gaussian
n <- 50
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
cvfit <- cv.grpreg(X, y, group, penalty='grLasso')
cvfit <- cv.grpreg(X, y, group, penalty='gel')
cvfit <- cv.grpreg(X, y, group, penalty='grLasso', fold=1:50)
cvfit <- cv.grpreg(X, y, group, penalty='gel', fold=1:50)
plot(cvfit, type='all')

# Binomial
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- runif(n) > 0.5
cvfit <- cv.grpreg(X, y, group, family='binomial', penalty='grLasso')
cvfit <- cv.grpreg(X, y, group, family='binomial', penalty='gel')
cvfit <- cv.grpreg(X, y, group, family='binomial', penalty='grLasso', fold=1:50)
cvfit <- cv.grpreg(X, y, group, family='binomial', penalty='gel', fold=1:50)
plot(cvfit, type='all')

# Poisson
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- sample(1:n)
cvfit <- cv.grpreg(X, y, group, family='poisson', penalty='grLasso')
cvfit <- cv.grpreg(X, y, group, family='poisson', penalty='gel')
cvfit <- cv.grpreg(X, y, group, family='poisson', penalty='grLasso', fold=1:50)
cvfit <- cv.grpreg(X, y, group, family='poisson', penalty='gel', fold=1:50)
plot(cvfit, type='all')
summary(cvfit)

# Multitask learning
n <- 50
p <- 10
m <- 4
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*m), ncol=m)
cvfit <- cv.grpreg(X, Y)
cvfit <- cv.grpreg(X, Y, nfolds=50)
Y <- matrix(rnorm(n*m), ncol=m) > 0
cvfit <- cv.grpreg(X, Y, family='binomial')
cvfit <- cv.grpreg(X, Y, family='binomial', nfolds=50)
plot(cvfit, type='all')

# p > n
n <- 75
p <- 200
X <- matrix(rnorm(n*p), n, p)
mu <- exp(apply(X[,1:10], 1, sum)*0.5)
y <- rpois(n, mu)
g <- rep(LETTERS[1:20], each=10)
cvfit <- cv.grpreg(X, y, group=g)
plot(cvfit, type='all')
cvfit <- cv.grpreg(X, y>0, group=g, family='binomial')
plot(cvfit, type='all')
cvfit <- cv.grpreg(X, y, group=g, family='poisson')
plot(cvfit, type='all')

# summary
set.seed(4)
n <- 75
p <- 200
X <- matrix(rnorm(n*p), n, p)
y <- rpois(n, 1)
g <- rep(LETTERS[1:20], each=10)
cvfit <- cv.grpreg(X, y, group=g)
s <- summary(cvfit)
expect_equivalent(s$ngroups[1], 0)

# predict
expect_equivalent(
  predict(cvfit, X[1:3,]),
  predict(cvfit$fit, X[1:3,], lambda=cvfit$lambda.min))


# R squared ---------------------------------------------------------------

# Gaussian
n <- 5000
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n, mean=X[,6])
cvfit <- cv.grpreg(X, y, group, penalty='grLasso', lambda.min=0)
summary(cvfit, lambda=0)$r.squared
summary(lm(y ~ ., as.data.frame(X)))$r.squared


# Return Y ----------------------------------------------------------------

n <- 50
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
cvfit <- cv.grpreg(X, y, group, returnY=TRUE)
expect_equivalent(
  cvfit$cve,
  apply((cvfit$Y - y)^2, 2, mean))
