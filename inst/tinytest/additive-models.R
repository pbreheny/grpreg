if (interactive()) library(tinytest)


# Linear regression -------------------------------------------------------

Data <- gen_nonlinear_data()
y <- Data$y

# Basic setup
Xb <- expand_spline(Data$X, df=4, type='bs')
expect_equal(ncol(Xb$X), ncol(Data$X)*4)
expect_equal(length(Xb$group), ncol(Data$X)*4)
expect_false(any(is.na(Xb$X)))
Xn <- expand_spline(Data$X, df=4, type='ns')
expect_equal(ncol(Xn$X), ncol(Data$X)*4)
expect_equal(length(Xn$group), ncol(Data$X)*4)
expect_false(any(is.na(Xn$X)))

# Fit
fit <- grpreg(Xn, y, penalty='grLasso', eps=1e-12)
fit2 <- grpreg(Xn$X, y, rep(1:ncol(Data$X), each=4), penalty='grLasso', eps=1e-12)
expect_equivalent(coef(fit), coef(fit2))
fit <- grpreg(Xb, y, penalty='grLasso', eps=1e-12)
fit2 <- grpreg(Xb$X, y, rep(1:ncol(Data$X), each=4), penalty='grLasso', eps=1e-12)
expect_equivalent(coef(fit), coef(fit2))

# Predict (bs)
fit <- grpreg(Xb, y, penalty='grLasso', eps=1e-12)
XX <- gen_nonlinear_data(seed=2, n=20)$X
expect_warning(P <- predict(fit, XX, type='link'))
expect_equal(dim(P), c(20, length(fit$lambda)))
P <- predict(fit, Data$X, type='link')
L <- apply(grpreg:::deviance_grpreg(y, P, 'gaussian'), 2, sum)
expect_equivalent(L, fit$deviance)
PP <- predict(fit, Data$X, type='response')
expect_equivalent(PP, P)
COEF <- predict(fit, type='coefficients')
expect_equivalent(COEF, coef(fit))
expect_inherits(predict(fit, type='vars', lambda=0.01), 'integer')
expect_inherits(predict(fit, type='groups', lambda=0.01), 'factor')
expect_inherits(predict(fit, type='nvars', lambda=0.01), 'integer')
expect_inherits(predict(fit, type='ngroups', lambda=0.01), 'integer')
N <- predict(fit, type='norm')
expect_true(all(dim(N) == c(ncol(Data$X), length(fit$lambda))))
expect_true(typeof(N) == 'double')

# Predict (ns)
fit <- grpreg(Xn, y, penalty='grLasso', eps=1e-12)
XX <- gen_nonlinear_data(seed=2, n=20)$X
P <- predict(fit, XX, type='link')
expect_equal(dim(P), c(20, length(fit$lambda)))
P <- predict(fit, Data$X, type='link')
L <- apply(grpreg:::deviance_grpreg(y, P, 'gaussian'), 2, sum)
expect_equivalent(L, fit$deviance)
PP <- predict(fit, Data$X, type='response')
expect_equivalent(PP, P)
COEF <- predict(fit, type='coefficients')
expect_equivalent(COEF, coef(fit))
expect_inherits(predict(fit, type='vars', lambda=0.01), 'integer')
expect_inherits(predict(fit, type='groups', lambda=0.01), 'factor')
expect_inherits(predict(fit, type='nvars', lambda=0.01), 'integer')
expect_inherits(predict(fit, type='ngroups', lambda=0.01), 'integer')
N <- predict(fit, type='norm')
expect_true(all(dim(N) == c(ncol(Data$X), length(fit$lambda))))
expect_true(typeof(N) == 'double')

# Plot (bs)
fit <- grpreg(Xb, y, penalty='grLasso', eps=1e-12)
plot_spline(fit, 'V02', lambda=0.01)
plot_spline(fit, 'V02', which=50)
plot_spline(fit, 'V02', which=80, partial=TRUE, type='conditional')
plot_spline(fit, 'V02', which=80, partial=TRUE, type='contrast')

# Plot (ns)
fit <- grpreg(Xn, y, penalty='grLasso', eps=1e-12)
plot_spline(fit, 'V02', lambda=0.01)
plot_spline(fit, 'V02', which=50)
plot_spline(fit, 'V02', which=80, partial=TRUE, type='conditional')
plot_spline(fit, 'V02', which=80, partial=TRUE, type='contrast')

# Cross-validation
cvfit <- cv.grpreg(Xn, y)
expect_silent(plot_spline(cvfit, 'V02'))
expect_warning(plot_spline(cvfit, 'V15', which=3))
plot_spline(cvfit, 'V02', partial=TRUE, type='conditional')
plot_spline(cvfit, 'V02', partial=TRUE, type='contrast')


# Logistic regression -----------------------------------------------------

Data <- gen_nonlinear_data(n=500)
y <- Data$y > quantile(Data$y, 0.25)

# Basic setup
Xb <- expand_spline(Data$X, df=4, type='bs')
Xn <- expand_spline(Data$X, df=4, type='ns')

# Fit
fit <- grpreg(Xn, y, penalty='grLasso', family='binomial', eps=1e-12)
fit2 <- grpreg(Xn$X, y, rep(1:ncol(Data$X), each=4), penalty='grLasso', family='binomial', eps=1e-12)
expect_equivalent(coef(fit), coef(fit2))
fit <- grpreg(Xb, y, penalty='grLasso', family='binomial', eps=1e-12)
fit2 <- grpreg(Xb$X, y, rep(1:ncol(Data$X), each=4), penalty='grLasso', family='binomial', eps=1e-12)
expect_equivalent(coef(fit), coef(fit2))

# Predict (bs)
fit <- grpreg(Xb, y, penalty='grLasso', family='binomial', eps=1e-12)
XX <- gen_nonlinear_data(seed=2, n=20)$X
XX[1,1] <- 1
expect_warning(P <- predict(fit, XX, type='link'))
expect_equal(dim(P), c(20, length(fit$lambda)))
P <- predict(fit, Data$X, type='response')
L <- apply(grpreg:::deviance_grpreg(y, P, 'binomial'), 2, sum)
expect_equivalent(L, fit$deviance, tolerance=1e-5)
COEF <- predict(fit, type='coefficients')
expect_equivalent(COEF, coef(fit))
expect_inherits(predict(fit, type='vars', lambda=0.01), 'integer')
expect_inherits(predict(fit, type='groups', lambda=0.01), 'factor')
expect_inherits(predict(fit, type='nvars', lambda=0.01), 'integer')
expect_inherits(predict(fit, type='ngroups', lambda=0.01), 'integer')
N <- predict(fit, type='norm')
expect_true(all(dim(N) == c(ncol(Data$X), length(fit$lambda))))
expect_true(typeof(N) == 'double')

# Plot (bs)
fit <- grpreg(Xb, y, penalty='grLasso', family='binomial', eps=1e-12)
plot_spline(fit, 'V02', lambda=0.01)
plot_spline(fit, 'V02', which=20)
plot_spline(fit, 'V02', which=20, partial=TRUE, type='conditional')
plot_spline(fit, 'V02', which=20, partial=TRUE, type='contrast')

# Cross-validation
cvfit <- cv.grpreg(Xn, y, family='binomial')
expect_silent(plot_spline(cvfit, 'V02'))
expect_warning(plot_spline(cvfit, 'V10', which=2))
plot_spline(cvfit, 'V02', partial=TRUE, type='conditional')
plot_spline(cvfit, 'V02', partial=TRUE, type='contrast')


# Cox regression -----------------------------------------------------

Data <- gen_nonlinear_data(n=500)
y <- cbind(Data$y, rbinom(500, 1, 0.75))

# Basic setup
Xb <- expand_spline(Data$X, df=4, type='bs')
Xn <- expand_spline(Data$X, df=4, type='ns')

# Fit
fit <- grpsurv(Xn, y, penalty='grLasso', eps=1e-12)
fit2 <- grpsurv(Xn$X, y, rep(1:ncol(Data$X), each=4), penalty='grLasso', eps=1e-12)
expect_equivalent(coef(fit), coef(fit2))
fit <- grpsurv(Xb, y, penalty='grLasso', eps=1e-12)
fit2 <- grpsurv(Xb$X, y, rep(1:ncol(Data$X), each=4), penalty='grLasso', eps=1e-12)
expect_equivalent(coef(fit), coef(fit2))

# Predict (bs)
fit <- grpsurv(Xb, y, penalty='grLasso', eps=1e-12)
XX <- gen_nonlinear_data(seed=2, n=20)$X
XX[1,1] <- 1
expect_warning(P <- predict(fit, XX, type='link'))
expect_equal(dim(P), c(20, length(fit$lambda)))
P <- predict(fit, Data$X, type='link')
L <- apply(grpreg:::deviance_grpsurv(y, P), 2, sum)
expect_equivalent(L, fit$deviance, tolerance=1e-5)
COEF <- predict(fit, type='coefficients')
expect_equivalent(COEF, coef(fit))
expect_inherits(predict(fit, type='vars', lambda=0.05), 'integer')
expect_inherits(predict(fit, type='groups', lambda=0.05), 'factor')
expect_inherits(predict(fit, type='nvars', lambda=0.05), 'integer')
expect_inherits(predict(fit, type='ngroups', lambda=0.05), 'integer')
N <- predict(fit, type='norm')
expect_true(all(dim(N) == c(ncol(Data$X), length(fit$lambda))))
expect_true(typeof(N) == 'double')

# Plot (bs)
fit <- grpsurv(Xb, y, penalty='grLasso', eps=1e-12)
plot_spline(fit, 'V02', lambda=0.05)
plot_spline(fit, 'V02', which=20)
plot_spline(fit, 'V02', which=20, partial=TRUE)
plot_spline(fit, 'V06', which=20, partial=TRUE)

# Cross-validation
cvfit <- cv.grpsurv(Xn, y)
expect_silent(plot_spline(cvfit, 'V02'))
expect_warning(plot_spline(cvfit, 'V11', which=2))
plot_spline(cvfit, 'V02', partial=TRUE)
