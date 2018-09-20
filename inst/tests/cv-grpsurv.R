.test = "Cross-validation: cox"
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- cbind(rexp(n), rep(0:1, c(10, n-10)))
cvfit <- cv.grpsurv(X, y, group, penalty='grLasso')
cvfit <- cv.grpsurv(X, y, group, penalty='gel')
cvfit <- cv.grpsurv(X, y, group, penalty='grLasso', nfolds=50)
cvfit <- cv.grpsurv(X, y, group, penalty='gel', nfolds=50)
cvfit <- cv.grpsurv(X, y, group, penalty='grLasso', se='bootstrap')
cvfit <- cv.grpsurv(X, y, group, penalty='gel', se='bootstrap')

plot(cvfit)
cvfit <- cv.grpsurv(X, y, group, penalty='grLasso')
