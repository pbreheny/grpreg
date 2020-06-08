n <- 20
p <- 5
l <- 5
group <- c(1,1,2,2,2)
zgroup <- rep(0:1, 2:3)

# standardize() standardizes correctly
X <- matrix(rnorm(n*p),ncol=p)
XX <- .Call("standardize", X)[[1]]
expect_equal(apply(XX,2,mean), rep(0,5))
expect_equal(apply(XX,2,crossprod), rep(20,5))

# unstandardize() unstandardizes correctly
X <- matrix(rnorm(n*p),ncol=p)
std <- .Call("standardize", X)
XX <- std[[1]]
center <- std[[2]]
scale <- std[[3]]
bb <- matrix(rnorm(l*(p+1)), nrow=p+1)
b <- grpreg:::unstandardize(bb, list(center=center, scale=scale, nz=which(scale>1e-10)))
expect_equal(cbind(1,XX) %*% bb, cbind(1,X) %*% b)

# orthogonalize() orthogonalizes correctly
X <- matrix(rnorm(n*p),ncol=p)
XX <- grpreg:::orthogonalize(X, group)
for (j in 1:group[p]) {
  ind <- which(group==j)
  expect_equal(crossprod(XX[,ind])/n, diag(length(ind)))
}

# unorthogonalize() unorthogonalizes correctly
X <- matrix(rnorm(n*p), ncol=p)
XX <- grpreg:::orthogonalize(X, group)
bb <- matrix(rnorm(l*(p+1)), nrow=p+1)
b <- grpreg:::unorthogonalize(bb, XX, attr(XX, "group"))
expect_equal(cbind(1,XX) %*% bb, cbind(1,X) %*% b)

# unorthogonalize() unorthogonalizes correctly w/o intercept
X <- matrix(rnorm(n*p), ncol=p)
XX <- grpreg:::orthogonalize(X, zgroup)
bb <- matrix(rnorm(l*(p)), nrow=p)
b <- grpreg:::unorthogonalize(bb, XX, attr(XX, "group"), intercept=FALSE)
expect_equal(XX %*% bb, X %*% b)

# orthogonalize() orthogonalizes correctly w/o full rank
X <- matrix(rnorm(n*p),ncol=p)
X[,5] <- X[,4]
XX <- grpreg:::orthogonalize(X, group)
for (j in 1:group[p]) {
  ind <- which(attr(XX, "group")==j)
  expect_equal(crossprod(XX[,ind])/n, diag(length(ind)))
}
y <- rnorm(nrow(X))
fit <- grpreg(X, y, group=LETTERS[group])

# orthogonalize() orthogonalizes correctly w/ 0's present and non-full-rank
X <- matrix(rnorm(n*p),ncol=p)
X[,4] <- 0
X[,5] <- X[,3]
grp <- grpreg:::newXG(X, group, rep(1, max(group)), 1, FALSE)
XX <- grp$X
g <- grp$g
for (j in 1:max(g)) {
  ind <- which(g==j)
  expect_equal(crossprod(XX[,ind])/n, diag(length(ind)))
}
y <- rnorm(nrow(X))
fit <- grpreg(X, y, group=LETTERS[group])
