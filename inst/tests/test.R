test_that("standardize() standardizes correctly", {
  X <- matrix(rnorm(500),ncol=10)
  XX <- standardize(X)
  expect_that(apply(XX,2,mean), equals(rep(0,10)))
  expect_that(apply(XX,2,crossprod), equals(rep(50,10)))
})

test_that("orthogonalize() orthogonalizes correctly", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  group <- c(1,1,2,2,2,3,3,3,3,4)
  XX <- orthogonalize(X, group)
  for (j in 1:group[p]) {
    ind <- which(group==j)
    expect_that(crossprod(XX[,ind])/n, equals(diag(length(ind))))
  }
})

X <- matrix(rnorm(500),ncol=10)
b <- rnorm(10)
y <- rnorm(X%*%b)
group=c(1,1,1,2,2,2,3,3,3,3)
tol <- .01

coef <- lm(y~X)$coef
beta <- grpreg(X,y,group=group,lambda=0,penalty="gLasso",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")
beta <- grpreg(X,y,group=group,lambda=0,penalty="gBridge",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")
beta <- grpreg(X,y,group=group,lambda=0,penalty="gMCP",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("gaussian check failed")

X <- matrix(rnorm(1000),ncol=10)
b <- rnorm(10)
y <- 1*(rnorm(X%*%b)>0)

coef <- glm(y~X,family="binomial")$coef
beta <- grpreg(X,y,group=group,family="binomial",lambda=0,penalty="gLasso",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("binomial check failed")
beta <- grpreg(X,y,group=group,family="binomial",lambda=0,penalty="gBridge",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("binomial check failed")
beta <- grpreg(X,y,group=group,family="binomial",lambda=0,penalty="gMCP",eps=.0001)$beta
if (max(abs(coef - beta)) > tol) stop("binomial check failed")

