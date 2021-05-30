Data <- gen_nonlinear_data()
y <- Data$y

# Basic setup
X <- expand_spline(Data$X, df=4)
expect_equal(ncol(X$X), ncol(Data$X)*4)
expect_equal(length(X$group), ncol(Data$X)*4)
expect_false(any(is.na(X$X)))

# Fit
fit <- grpreg(X, y, penalty='grLasso', eps=1e-12)
fit2 <- grpreg(X$X, y, rep(1:ncol(Data$X), each=4), penalty='grLasso', eps=1e-12)
expect_equivalent(coef(fit), coef(fit2))  # Passes

# Predict
XX <- gen_nonlinear_data(seed=2, n=20)$X
expect_warning(predict(fit, XX, type='link'))
P <- predict(fit, Data$X, type='link')
L <- apply(grpreg:::loss.grpreg(y, P, 'gaussian'), 2, sum)
expect_equivalent(L, fit$loss)
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

# Plot
plot_spline(fit, 'V02', lambda=0.01)
plot_spline(fit, 'V02', which=50)
plot_spline(fit, 'V02', which=80, scatter=TRUE)

# Cross-validation
B <- expand_spline(Data$X, type='ns')
cvfit <- cv.grpreg(B, y)

# Logistic

