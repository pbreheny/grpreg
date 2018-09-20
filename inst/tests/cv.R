.test = "Cross-validation: gaussian"
n <- 50
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
cvfit <- cv.grpreg(X, y, group, penalty='grLasso')
cvfit <- cv.grpreg(X, y, group, penalty='gel')
cvfit <- cv.grpreg(X, y, group, penalty='grLasso', nfolds=50)
cvfit <- cv.grpreg(X, y, group, penalty='gel', nfolds=50)

.test = "Cross-validation: binomial"
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- runif(n) > 0.5
cvfit <- cv.grpreg(X, y, group, family='binomial', penalty='grLasso')
cvfit <- cv.grpreg(X, y, group, family='binomial', penalty='gel')
cvfit <- cv.grpreg(X, y, group, family='binomial', penalty='grLasso', nfolds=50)
cvfit <- cv.grpreg(X, y, group, family='binomial', penalty='gel', nfolds=50)

.test = "Cross-validation: poisson"
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- sample(1:n)
cvfit <- cv.grpreg(X, y, group, family='poisson', penalty='grLasso')
cvfit <- cv.grpreg(X, y, group, family='poisson', penalty='gel')
cvfit <- cv.grpreg(X, y, group, family='poisson', penalty='grLasso', nfolds=50)
cvfit <- cv.grpreg(X, y, group, family='poisson', penalty='gel', nfolds=50)

.test = "Cross-validation: multitask learning"
n <- 50
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
cvfit <- cv.grpreg(X, y, group, returnY=TRUE)
cve <- apply(cvfit$Y - y, 2, crossprod)/n
check(cve, cvfit$cve, tol= .001)
y <- rnorm(n) > 0
cvfit <- cv.grpreg(X, y, group, family='binomial', returnY=TRUE, lambda.min=0.5)
pe <- apply((cvfit$Y>0.5)!=y, 2, mean)
check(pe, cvfit$pe, tol= .001)

.test = "Cross-validation: p > n"
n <- 75
p <- 200
X <- matrix(rnorm(n*p), n, p)
mu <- exp(apply(X[,1:10], 1, sum))
y <- rpois(n, mu)
g <- rep(LETTERS[1:20], each=10)
cvfit <- cv.grpreg(X, y, group=g)
cvfit <- cv.grpreg(X, y>0, group=g, family='binomial')
cvfit <- cv.grpreg(X, y, group=g, family='poisson')
