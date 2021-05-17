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
X <- grpmat(X.raw, df=4)
expect_equal(ncol(X$x), ncol(X.raw)*4)
expect_equal(length(X$groups), ncol(X.raw)*4)
expect_equal(any(is.na(X)), FALSE) #expect_false

