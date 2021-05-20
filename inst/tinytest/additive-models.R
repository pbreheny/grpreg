# Generate test data
n <- 100
J <- 16
X.raw <- matrix(runif(n*J), nrow=n, ncol=J)
f <- vector("list",6)
f[[1]] <- function(x){2*(exp(-10*x)-exp(-10))/(1-exp(-10)) - 1}
f[[2]] <- function(x){-2*(exp(-10*x)-exp(-10))/(1-exp(-10)) + 1}
f[[3]] <- function(x){2*x-1}
f[[4]] <- function(x){-2*x+1}
f[[5]] <- function(x){8*(x-0.5)^2 - 1}
f[[6]] <- function(x){-8*(x-0.5)^2 + 1}
eta <- matrix(NA, nrow=n, ncol=6)
for (j in 1:6) eta[,j] <- f[[j]](X.raw[,j])
mu <- apply(eta,1,sum)
y <- rnorm(n, mean=mu)

# Basic setup
X <- expand_spline(X.raw, df=4)
expect_equal(ncol(X$X), ncol(X.raw)*4)
expect_equal(length(X$group), ncol(X.raw)*4)
expect_false(any(is.na(X$X)))

# Fit
fit <- grpreg(X, y, penalty='grLasso', eps=1e-12)
fit2 <- grpreg(X$X, y, rep(1:ncol(X.raw), each=4), penalty='grLasso', eps=1e-12)
expect_equivalent(coef(fit), coef(fit2))  # Passes

# Plot
plot_spline(fit, 'V2', lambda=0.01)
plot_spline(fit, 'V2', which=50)

# Cross-validation
B <- expand_spline(X.raw, type='ns')
cvfit <- cv.grpreg(B, y)
