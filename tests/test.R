library(grpreg)
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

