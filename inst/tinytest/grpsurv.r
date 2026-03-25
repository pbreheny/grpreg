if (interactive()) library(tinytest)
suppressPackageStartupMessages(library(survival))

# Test that grpsurv works when x has 1 column
y <- Surv(rexp(50), sample(rep(0:1, c(10, 40))))
x <- matrix(rnorm(50 * 1), 50, 1)
g <- 1
fit <- grpsurv(x, y, g, lambda.min = 0)

# Test that grpsurv works
y <- Surv(rexp(50), sample(rep(0:1, c(10, 40))))
x <- matrix(rnorm(50 * 6), 50, 6)
g <- rep(1:3, each = 2)
fit <- grpsurv(x, y, g, lambda.min = 0)

# $ grpsurv equals MLE when lam=0
fit_mle <- coxph(y ~ x)
expect_equivalent(coef(fit)[, 100], coef(fit_mle), tol = 0.01)
fit <- grpsurv(x, y, g, lambda.min = 0, penalty = "grSCAD")
expect_equivalent(coef(fit)[, 100], coef(fit_mle), tol = 0.01)
fit <- grpsurv(x, y, g, lambda.min = 0, penalty = "grMCP")
expect_equivalent(coef(fit)[, 100], coef(fit_mle), tol = 0.01)
fit <- grpsurv(x, y, g, lambda.min = 0, penalty = "gel")
expect_equivalent(coef(fit)[, 100], coef(fit_mle), tol = 0.01)
fit <- grpsurv(x, y, g, lambda.min = 0, penalty = "cMCP")
expect_equivalent(coef(fit)[, 100], coef(fit_mle), tol = 0.01)

# grpsurv returns correct logLik
expect_equivalent(logLik(fit)[100], logLik(fit_mle)[1], tol = 0.1)
expect_equivalent(AIC(fit)[100], AIC(fit_mle), tol = 0.1)

# grpsurv returns correct linear predictors
fit <- grpsurv(x, y, g, lambda.min = 0, penalty = "grLasso", eps = 1e-12)
expect_equivalent(fit$linear.predictors[, 100], fit_mle$linear.predictors[order(y)])
fit <- grpsurv(x, y, g, lambda.min = 0, penalty = "gel", eps = 1e-12)
expect_equivalent(fit$linear.predictors[, 100], fit_mle$linear.predictors[order(y)])

# residuals are correct (slightly different at final observation because
# baseline hazard estimated differently)
fit <- grpsurv(x, y, g, lambda.min = 0, penalty = "grLasso", eps = 1e-12)
r1 <- residuals(fit, lambda = 0)
r2 <- residuals(fit_mle, type = "deviance")
plot(r1, r2)
abline(0, 1)
expect_equivalent(
  residuals(fit, lambda = 0),
  residuals(fit_mle, type = "deviance"),
  tolerance = 0.1
)
fit <- grpsurv(x, y, g, lambda.min = 0, penalty = "gel", eps = 1e-12)
expect_equivalent(
  residuals(fit, lambda = 0),
  residuals(fit_mle, type = "deviance"),
  tolerance = 0.1
)

# predict works for grpsurv
head(predict(fit, x, "vars"))
head(predict(fit, x, "nvars"))
head(predict(fit, x, "groups"))
head(predict(fit, x, "ngroups"))
head(predict(fit, x, "link"))
head(predict(fit, x, "response"))
head(predict(fit, x, "coef"))
head(predict(fit, x, "median"))
s <- predict(fit, x, "survival", lambda = 0.1)
plot(s)
s <- predict(fit, x[1, ], "survival", lambda = 0.1)
plot(s)
h <- predict(fit, x, "hazard", lambda = 0.1)
plot(h)
h <- predict(fit, x[1, ], "hazard", lambda = 0.1)
plot(h)

# Survival curves vs survival package
df <- as.data.frame(x)
fit_mle <- coxph(y ~ ., df)
s <- predict(fit, x[1, ], "survival", lambda = 0)
plot(s)
lines(survfit(fit_mle, df[1, ], se.fit = FALSE))
h <- predict(fit, x[1, ], "hazard", lambda = 0)
plot(h)
lines(survfit(fit_mle, df[1, ]), se.fit = FALSE, cumhaz = TRUE)

# Cumulative hazard
h1 <- h(fit$time)
h2 <- survfit(fit_mle, df[1, ])$cumhaz
plot(h1, h2)
abline(0, 1)

# penalty factor works for grpsurv
fit <- grpsurv(x, y, g, group.multiplier = 1:3)

# coersion works for grpsurv
fit <- grpsurv(as.data.frame(x), y)

# loss works for grpsurv
eta <- predict(fit, x, "link", lambda = 0.1)
grpreg:::deviance_grpsurv(y, eta)
grpreg:::deviance_grpsurv(y, eta, total = FALSE)

# cross-validation works for grpsurv
cvfit <- cv.grpsurv(x, y, g, lambda.min = 0)
plot(cvfit)
summary(cvfit)
