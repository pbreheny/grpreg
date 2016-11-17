rm(list = ls())
require(grpreg)
data(Birthwt)
X <- Birthwt$X
group <- Birthwt$group

## Linear regression
y <- Birthwt$bwt
fit <- grpreg(X, y, group, penalty="grLasso")
plot(fit)

fit.ssr <- grpreg(X, y, group, penalty = "grLasso", screen = "SSR")
plot(fit.ssr)

all.equal(fit$loss, fit.ssr$loss)
all.equal(fit$beta, fit.ssr$beta)
