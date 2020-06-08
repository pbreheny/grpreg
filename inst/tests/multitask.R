n <- 100
p <- 5

# multitask learning works
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*3), ncol=3)
colnames(X) <- LETTERS[1:p]
fit <- grpreg(X, Y, penalty="grLasso")
expect_equal(dim(fit$beta), c(3, p+1, 100))
fit <- grpreg(X, Y, penalty="cMCP")
expect_equal(dim(fit$beta), c(3, p+1, 100))
fit <- gBridge(X, Y)
expect_equal(dim(fit$beta), c(3, p+1, 100))

# multitask learning works (logistic regression)
fit <- grpreg(X, Y>0, family="binomial", penalty="grLasso", lambda.min=0.4)
expect_equal(dim(fit$beta), c(3, p+1, 100))
fit <- grpreg(X, Y>0, family="binomial", penalty="cMCP", lambda.min=0.4)
expect_equal(dim(fit$beta), c(3, p+1, 100))
fit <- gBridge(X, Y>0, family="binomial", lambda.min=0.4)
expect_equal(dim(fit$beta), c(3, p+1, 100))

# coef/predict work for multitask learning
fit <- grpreg(X, Y)
l <- fit$lambda[20]
coef(fit, which=1:2)
coef(fit, lambda=l)
predict(fit, lambda=l, type="nvars")
predict(fit, which=c(30,60), type="nvars")
predict(fit, lambda=l, type="ngroups")
predict(fit, which=c(30,60), type="ngroups")
predict(fit, lambda=l, type="groups")
predict(fit, which=c(30,60), type="groups")
predict(fit, lambda=l, type="norm")
predict(fit, which=c(30,60), type="norm")
head(predict(fit, X, lambda=l))
predict(fit, X, which=c(30,60))[1:10,,]

# coef/predict work for multitask learning (logistic regression)
fit <- grpreg(X, Y>0, family="binomial")
l <- fit$lambda[20]
predict(fit, lambda=l, type="nvars")
predict(fit, lambda=l, type="ngroups")
predict(fit, lambda=l, type="groups")
predict(fit, lambda=l, type="norm")
head(predict(fit, X, lambda=l))
head(predict(fit, X, lambda=l, type="response"))
head(predict(fit, X, lambda=l, type="class"))

# cross-validation for multitask learning works
cvfit <- cv.grpreg(X, Y)
plot(cvfit)
summary(cvfit)
cvfit <- cv.grpreg(X, Y, penalty="cMCP")
plot(cvfit)
summary(cvfit)

# cross-validation for multitask learning works (logistic regression)
cvfit <- cv.grpreg(X, Y>0, family="binomial")
plot(cvfit)
summary(cvfit)
cvfit <- cv.grpreg(X, Y>0, family="binomial", penalty="cMCP")
plot(cvfit)
summary(cvfit)

# multitask learning reproduces linear regression
fit.mle <- lm(Y~X)
reg <- coef(fit.mle)
cMCP <- coef(fit <- grpreg(X, Y, penalty="cMCP", lambda.min=0), which=100)
expect_equal(t(cMCP), reg, tolerance=0.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)
expect_equal(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)
bridge <- coef(fit <- gBridge(X, Y, lambda.min=0), which=1)
expect_equal(t(bridge), reg, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=1)
expect_equal(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)
grLasso <- coef(fit <- grpreg(X, Y, penalty="grLasso", lambda.min=0), which=100)
expect_equal(t(grLasso), reg, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)
expect_equal(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)
grMCP <- coef(fit <- grpreg(X, Y, penalty="grMCP", lambda.min=0), which=100)
expect_equal(t(grMCP), reg, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)
expect_equal(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)
grSCAD <- coef(fit <- grpreg(X, Y, penalty="grSCAD", lambda.min=0), which=100)
expect_equal(t(grSCAD), reg, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)
expect_equal(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)

# multitask learning reproduces logistic regression
n <- 50
p <- 2
X <- matrix(rnorm(n*p),ncol=p)
Y <- matrix(rnorm(n*3),ncol=3)>0
fit.mle <- glm(Y[,3]~X, family=binomial)
mle <- coef(fit.mle)
beta <- coef(fit <- grpreg(X, Y, lambda.min=0, family="binomial"), which=100)[3,]
expect_equal(beta, mle, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)[,3]
expect_equal(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100, type="response")[,3]
expect_equal(p, predict(fit.mle, type="response"), tolerance=.01, check.attributes=FALSE)
bridge <- coef(fit <- gBridge(X, Y, family="binomial", lambda.min=0), which=1)[3,]
expect_equal(bridge, mle, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=1)[,3]
expect_equal(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=1, type="response")[,3]
expect_equal(p, predict(fit.mle, type="response"), tolerance=.01, check.attributes=FALSE)
