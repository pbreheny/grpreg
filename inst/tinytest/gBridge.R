# gBridge reproduces linear regression
n <- 50
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
group <- rep(0:3,4:1)
fit.mle <- lm(y~X)
reg <- coef(fit.mle)
fit <- gBridge(X, y, group, lambda.min=0, eps=1e-10)
expect_equivalent(coef(fit, which=1), coef(fit.mle))
expect_silent(plot(fit))

# CV
cvfit <- cv.grpreg(X, y, group, lambda.min=0, gBridge=TRUE)
plot(cvfit)
summary(cvfit)

# Should probably have more tests for logistic regression, etc.
