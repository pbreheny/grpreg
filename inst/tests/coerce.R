.test = "coersion of X, y"
data(Birthwt)
X <- data.frame(Birthwt$X)
y <- factor(Birthwt$low, labels=c("No", "Yes"))
fit <- grpreg(X, y, family="binomial")

.test = "coersion of group"
y <- Birthwt$low
g1 <- Birthwt$group
g2 <- as.numeric(factor(g1))
g3 <- as.numeric(factor(g1, levels=unique(g1)))
fit1 <- grpreg(X, y, group=g1, family="binomial")
fit2 <- grpreg(X, y, group=g2, family="binomial")
fit3 <- grpreg(X, y, group=g3, family="binomial")
check(coef(fit1, which=50), coef(fit2, which=50), tol=0.001)
check(coef(fit2, which=50), coef(fit3, which=50), tol=0.001)
check(coef(fit1, which=50), coef(fit3, which=50), tol=0.001)
