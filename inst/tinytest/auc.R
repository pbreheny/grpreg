if (interactive()) library(tinytest)
library(survival, quietly=TRUE)

# The quantities here are not exactly the same thing (one is in-sample, the
# other is out-of-sample), but they should be close
n <- 500
X <- matrix(rnorm(n*10), n, 10)
y <- Surv(rexp(n, rate=exp(X[,1])), rbinom(n, 1, prob=0.8))
cvfit <- cv.grpsurv(X, y, lambda.min=0, returnY=TRUE)
a <- AUC(cvfit)
a[length(a)]
fit <- coxph(y~X)
b <- concordancefit(y, -predict(fit))$concordance
expect_equivalent(a[length(a)], b, tol=0.02)
