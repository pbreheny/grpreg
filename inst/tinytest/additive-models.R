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

# Plot
plot_spline(fit, 'V02', lambda=0.01)
plot_spline(fit, 'V02', which=50)
# plot_spline(fit, 'V02', which=50, scatter=TRUE)

# Cross-validation
B <- expand_spline(Data$X, type='ns')
cvfit <- cv.grpreg(B, y)
