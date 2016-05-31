set.seed(1)
n <- 20
p <- 5
l <- 5
group <- c(1,1,2,2,2)
zgroup <- rep(0:1, 2:3)

.test = "standardize() standardizes correctly"
X <- matrix(rnorm(n*p),ncol=p)
XX <- .Call("standardize", X)[[1]]
check(apply(XX,2,mean), rep(0,5))
check(apply(XX,2,crossprod), rep(20,5))

.test = "unstandardize() unstandardizes correctly"
X <- matrix(rnorm(n*p),ncol=p)
std <- .Call("standardize", X)
XX <- std[[1]]
center <- std[[2]]
scale <- std[[3]]
bb <- matrix(rnorm(l*(p+1)), nrow=p+1)
b <- grpreg:::unstandardize(bb, center, scale)
check(cbind(1,XX) %*% bb, cbind(1,X) %*% b)

.test = "orthogonalize() orthogonalizes correctly"
X <- matrix(rnorm(n*p),ncol=p)
XX <- grpreg:::orthogonalize(X, group)
for (j in 1:group[p]) {
  ind <- which(group==j)
  check(crossprod(XX[,ind])/n, diag(length(ind)))
}

.test = "unorthogonalize() unorthogonalizes correctly"
X <- matrix(rnorm(n*p), ncol=p)
XX <- grpreg:::orthogonalize(X, group)
bb <- matrix(rnorm(l*(p+1)), nrow=p+1)
b <- grpreg:::unorthogonalize(bb, XX, attr(XX, "group"))
check(cbind(1,XX) %*% bb, cbind(1,X) %*% b)

.test = "unorthogonalize() unorthogonalizes correctly w/o intercept"
X <- matrix(rnorm(n*p), ncol=p)
XX <- grpreg:::orthogonalize(X, zgroup)
bb <- matrix(rnorm(l*(p)), nrow=p)
b <- grpreg:::unorthogonalize(bb, XX, attr(XX, "group"), intercept=FALSE)
check(XX %*% bb, X %*% b)

.test = "orthogonalize() orthogonalizes correctly w/ 0's present"
X <- matrix(rnorm(n*p),ncol=p)
XX <- grpreg:::orthogonalize(X, zgroup)
for (j in 1:zgroup[p]) {
  ind <- which(zgroup==j)
  check(crossprod(XX[,ind])/n, diag(length(ind)))
}

.test = "orthogonalize() orthogonalizes correctly w/o full rank"
X <- matrix(rnorm(n*p),ncol=p)
X[,5] <- X[,4]
XX <- grpreg:::orthogonalize(X, group)
for (j in 1:group[p]) {
  ind <- which(attr(XX, "group")==j)
  check(crossprod(XX[,ind])/n, diag(length(ind)))
}
