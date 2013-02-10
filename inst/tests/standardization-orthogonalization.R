test_that("standardize() standardizes correctly", {
  X <- matrix(rnorm(500),ncol=10)
  XX <- standardize(X)
  expect_that(apply(XX,2,mean), equals(rep(0,10)))
  expect_that(apply(XX,2,crossprod), equals(rep(50,10)))
})

test_that("unstandardize() unstandardizes correctly", {
  n <- 50
  p <- 10
  l <- 5
  X <- matrix(rnorm(n*p),ncol=p)
  XX <- standardize(X)
  bb <- matrix(rnorm(l*(p+1)), nrow=p+1)
  b <- unstandardize(bb, attr(XX,"center"), attr(XX, "scale"))
  expect_that(cbind(1,XX) %*% bb, equals(cbind(1,X) %*% b))
})

test_that("orthogonalize() orthogonalizes correctly", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  group <- c(1,1,2,2,2,3,3,3,3,4)
  XX <- orthogonalize(X, group)
  for (j in 1:group[p]) {
    ind <- which(group==j)
    expect_that(crossprod(XX[,ind])/n, equals(diag(length(ind))))
  }
})

test_that("unorthogonalize() unorthogonalizes correctly", {
  n <- 50
  p <- 10
  l <- 5
  group <- c(1,1,2,2,2,3,3,3,3,4)
  X <- matrix(rnorm(n*p), ncol=p)
  XX <- orthogonalize(X, group)
  bb <- matrix(rnorm(l*(p+1)), nrow=p+1)
  b <- unorthogonalize(bb, XX, group, attr(XX, "group"))
  expect_that(cbind(1,XX) %*% bb, equals(cbind(1,X) %*% b))
})

test_that("orthogonalize() orthogonalizes correctly w/ 0's present", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  group <- rep(0:3,4:1)
  XX <- orthogonalize(X, group)
  for (j in 1:group[p]) {
    ind <- which(group==j)
    expect_that(crossprod(XX[,ind])/n, equals(diag(length(ind))))
  }
})

test_that("orthogonalize() orthogonalizes correctly w/o full rank", {
  n <- 50
  p <- 10
  X <- matrix(rnorm(n*p),ncol=p)
  group <- rep(0:3,4:1)
  X[,7] <- X[,6]
  XX <- orthogonalize(X, group)
  for (j in 1:group[p]) {
    ind <- which(attr(XX, "group")==j)
    expect_that(crossprod(XX[,ind])/n, equals(diag(length(ind))))
  }
})
