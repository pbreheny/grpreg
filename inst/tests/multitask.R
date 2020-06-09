n <- 100
p <- 5

# multitask learning works
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*3), ncol=3)
colnames(X) <- LETTERS[1:p]
fit <- grpreg(X, Y, penalty="grLasso")
expect_equivalent(dim(fit$beta), c(3, p+1, 100))
fit <- grpreg(X, Y, penalty="cMCP")
expect_equivalent(dim(fit$beta), c(3, p+1, 100))
fit <- gBridge(X, Y)
expect_equivalent(dim(fit$beta), c(3, p+1, 100))

# multitask learning works (logistic regression)
fit <- grpreg(X, Y>0, family="binomial", penalty="grLasso", lambda.min=0.4)
expect_equivalent(dim(fit$beta), c(3, p+1, 100))
fit <- grpreg(X, Y>0, family="binomial", penalty="cMCP", lambda.min=0.4)
expect_equivalent(dim(fit$beta), c(3, p+1, 100))
fit <- gBridge(X, Y>0, family="binomial", lambda.min=0.4)
expect_equivalent(dim(fit$beta), c(3, p+1, 100))

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
expect_equivalent(t(cMCP), reg, tolerance=0.01)
p <- predict(fit, X, which=100)
expect_equivalent(p, predict(fit.mle), tolerance=.01)
bridge <- coef(fit <- gBridge(X, Y, lambda.min=0), which=1)
expect_equivalent(t(bridge), reg, tolerance=.01)
p <- predict(fit, X, which=1)
expect_equivalent(p, predict(fit.mle), tolerance=.01)
grLasso <- coef(fit <- grpreg(X, Y, penalty="grLasso", lambda.min=0), which=100)
expect_equivalent(t(grLasso), reg, tolerance=.01)
p <- predict(fit, X, which=100)
expect_equivalent(p, predict(fit.mle), tolerance=.01)
grMCP <- coef(fit <- grpreg(X, Y, penalty="grMCP", lambda.min=0), which=100)
expect_equivalent(t(grMCP), reg, tolerance=.01)
p <- predict(fit, X, which=100)
expect_equivalent(p, predict(fit.mle), tolerance=.01)
grSCAD <- coef(fit <- grpreg(X, Y, penalty="grSCAD", lambda.min=0), which=100)
expect_equivalent(t(grSCAD), reg, tolerance=.01)
p <- predict(fit, X, which=100)
expect_equivalent(p, predict(fit.mle), tolerance=.01)

# multitask learning reproduces logistic regression
n <- 50
p <- 2
X <- matrix(rnorm(n*p),ncol=p)
Y <- matrix(rnorm(n*3),ncol=3)>0
fit.mle <- glm(Y[,3]~X, family=binomial)
mle <- coef(fit.mle)
beta <- coef(fit <- grpreg(X, Y, lambda.min=0, family="binomial"), which=100)[3,]
expect_equivalent(beta, mle, tolerance=.01)
p <- predict(fit, X, which=100)[,3]
expect_equivalent(p, predict(fit.mle), tolerance=.01)
p <- predict(fit, X, which=100, type="response")[,3]
expect_equivalent(p, predict(fit.mle, type="response"), tolerance=.01)
bridge <- coef(fit <- gBridge(X, Y, family="binomial", lambda.min=0), which=1)[3,]
expect_equivalent(bridge, mle, tolerance=.01)
p <- predict(fit, X, which=1)[,3]
expect_equivalent(p, predict(fit.mle), tolerance=.01)
p <- predict(fit, X, which=1, type="response")[,3]
expect_equivalent(p, predict(fit.mle, type="response"), tolerance=.01)
