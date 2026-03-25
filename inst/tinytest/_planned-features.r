# This are tests that do not currently pass, but I would like to improve the
# package to the point where they do pass

# group.multiplier works to assign unpenalized groups
n <- 50
p <- 11
x <- matrix(rnorm(n * p),ncol = p)
y <- rnorm(n)
group <- rep(0:3, c(1, 2, 3, 5))
gm <- 0:2
plot(fit <- grpreg(x, y, group, penalty = "cMCP", lambda.min = 0, group.multiplier = gm), main = fit$penalty)
fit <- grpreg(X, y, group, penalty = "cMCP", lambda.min = 0, group.multiplier = gm)

