cv.grpreg <- function(X, y, group, penalty=c("gMCP","gLasso","gBridge"), family=c("gaussian","binomial"), alpha=1, lambda.min=ifelse(n>p,.001,.05), lambda.max, nlambda=100, lambda, nfolds=10, seed, trace=FALSE, gamma, group.multiplier=rep(1,J),...)
  {
    ## Check for errors
    family <- match.arg(family)
    penalty <- match.arg(penalty)
    if (alpha > 1 | alpha < 0) stop("alpha must be in [0,1]")
    if (length(group)!=ncol(X)) stop("group does not match X")
    if (is.null(colnames(X))) colnames(X) <- paste("V",1:ncol(X),sep="")
    J <- max(group)
    K <- as.numeric(table(group))
    if (!(identical(as.integer(sort(unique(group))),as.integer(1:J)) | identical(as.integer(sort(unique(group))),as.integer(0:J)))) stop("Groups must be consecutively numbered 1,2,3,...")
    if (length(group.multiplier)!=J) stop("Length of group.multiplier must equal number of penalized groups")
    if (missing(gamma))
      {
        if (penalty=="gBridge") gamma=0.5
        else if (family=="gaussian") gamma=3
        else gamma=30
      }
    if (!missing(seed)) set.seed(seed)

    ## Set up XX, yy, lambda
    n <- nrow(X)
    meanx <- apply(X,2,mean)
    normx <- sqrt(apply((t(X)-meanx)^2,1,sum))/sqrt(n)
    if (any(normx < 0.0001)) stop("X contains columns which are numerically constant.  If you are attempting to specify an intercept, please remove this column; an intercept is included automatically.")
    XX <- cbind(1,scale(X,meanx,normx))
    group <- c(0,group)
    colnames(XX)[1] <- "(Intercept)"
    p <- ncol(XX)
    if (missing(lambda))
      {
        lambda <- setupLambda(XX,y,group,family,penalty,lambda.max,lambda.min,nlambda,gamma,group.multiplier)
        user.lambda <- FALSE
      }
    else
      {
        nlambda <- length(lambda)
        user.lambda <- TRUE
      }
    rm(XX)

    error <- mce <- array(NA,dim=c(nfolds,length(lambda)))
    
    if (family=="gaussian")
      {
        cv.ind <- ceiling(sample(1:n)/n*nfolds)
      }
    else if (family=="binomial")
      {
        ind1 <- which(y==1)
        ind0 <- which(y==0)
        n1 <- length(ind1)
        n0 <- length(ind0)
        cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
        cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
        cv.ind <- numeric(n)
        cv.ind[y==1] <- cv.ind1
        cv.ind[y==0] <- cv.ind0
      }

    for (i in 1:nfolds)
      {
        if (trace) cat("Starting CV fold #",i,sep="","\n")
        X1 <- X[cv.ind!=i,]
        y1 <- y[cv.ind!=i]
        X2 <- X[cv.ind==i,]
        y2 <- y[cv.ind==i]

        fit.i <- grpreg(X1,y1,penalty=penalty,family=family,alpha=alpha,lambda=lambda,warn=FALSE,gamma=gamma,...)
        yhat <- predict(fit.i,X2,type="response")
        error[i,1:ncol(yhat)] <- loss.grpreg(y2,yhat,family)
        print(fit.i$df)
        mce[i,1:ncol(yhat)] <- apply(y2!=(yhat>.5),2,sum)
      }

    ## Eliminate saturated lambda values, if any
    ind <- which(apply(is.finite(error),2,all))
    E <- error[,ind]
    lambda <- lambda[ind]
    mce <- mce[,ind]
    
    val <- list(E=E,
                cve=apply(E,2,mean),
                lambda=lambda,
                mce=apply(mce,2,mean))
    val$min <- which.min(val$cve)
    class(val) <- "cv.ncvreg"
    return(val)
  }

