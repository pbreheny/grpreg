if (interactive()) library(tinytest)

data(Birthwt)
X <- Birthwt$X
y <- Birthwt$bwt
group <- Birthwt$group
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0)
lmfit <- lm(y~X)

# BIC
expect_silent(sb <- select(fit))
expect_equal(sb$IC[100], BIC(lmfit), tolerance=1e-4)

# AIC
expect_warning(sa <- select(fit, criterion='AIC', df.method='active'))
expect_silent(sa <- select(fit, criterion='AIC'))
expect_equivalent(sa$IC[100], AIC(lmfit), tolerance=1e-4)

# GCV
expect_silent(s <- select(fit, criterion='GCV', smooth=TRUE))

# AICc
expect_silent(s <- select(fit, criterion='AICc'))
expect_true(all(s$IC >= sa$IC))

# EBIC
expect_silent(s <- select(fit, criterion='EBIC'))
expect_true(all(s$IC >= sb$IC))

