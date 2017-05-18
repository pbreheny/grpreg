.test = "coersion of y"
char_y <- matrix(LETTERS[1:2], 4, 3)
logi_y <- matrix(rep(c(TRUE, FALSE), each=6), 4, 3)
int_y <- matrix(1:12, 4, 3)
num_y <- matrix(1.0*1:12, 4, 3)
tools::assertError(grpreg:::newY(char_y, 'gaussian'))
tools::assertError(grpreg:::newY(num_y, 'binomial'))
tools::assertError(grpreg:::newY(int_y, 'binomial'))
grpreg:::newY(logi_y, 'binomial')
grpreg:::newY(char_y, 'binomial')
grpreg:::newY(num_y, 'gaussian')
grpreg:::newY(int_y, 'gaussian')
grpreg:::newY(char_y[,1], 'binomial')

.test = "coersion of X, y"
data(Birthwt)
X <- data.frame(Birthwt$X)
y <- factor(Birthwt$low, labels=c("No", "Yes"))
fit <- grpreg(X, y, family="binomial")

.test = "coersion of group"
y <- Birthwt$low
g1 <- Birthwt$group
g2 <- as.numeric(factor(g1))
g3 <- as.numeric(factor(g1, levels=sort(levels(g1)))) ## Tests reordering
fit1 <- grpreg(X, y, group=g1, family="binomial")
fit2 <- grpreg(X, y, group=g2, family="binomial")
fit3 <- grpreg(X, y, group=g3, family="binomial")
check(coef(fit1, which=50), coef(fit2, which=50), tol=0.001)
check(coef(fit2, which=50), coef(fit3, which=50), tol=0.001)
check(coef(fit1, which=50), coef(fit3, which=50), tol=0.001)
