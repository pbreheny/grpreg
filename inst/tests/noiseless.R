X <- rbind(diag(8), -diag(8)) * sqrt(8)
beta <- c(-1,-1,-1,-1,2,2,2,2)
y <- as.numeric(X %*% beta)
group <- c(1,1,1,1,2,2,2,2)

fit <- grpreg(X, y, group, lambda=c(2, 1.999999, 1, 0))
expect_equal(predict(fit, type='nvars', which=2), 4)
expect_equivalent(coef(fit, lambda=0), c(0, beta))
expect_equivalent(coef(fit, lambda=1), c(0, 0, 0, 0, 0, 1, 1, 1, 1))
