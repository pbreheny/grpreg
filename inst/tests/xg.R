## source("~/dev/.grpreg.setup.R")
## g1 <- rep(1:3, 1:3)
## g2 <- c(2, 1, 3, 3, 3, 2)
## m <- 1:3

## # Ordinary X
## X <- matrix(rnorm(10*6), 10, 6)
## newXG(X, g1, m, 1, TRUE)
## newXG(X, g2, m, 1, TRUE)

## # X with constant column
## X <- matrix(rnorm(10*6), 10, 6)
## X[,6] <- 0
## newXG(X, g1, m, 1, TRUE)
## newXG(X, g2, m, 1, TRUE)
