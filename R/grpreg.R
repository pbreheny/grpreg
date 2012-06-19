grpreg <- function(X, y, group=1:ncol(X), family=c("gaussian","binomial"), penalty=c("geMCP","geLasso","gMCP","gBridge","gLasso"), nlambda=100, lambda, lambda.min=ifelse(n>p,.001,.05), lambda.max, alpha=1, tau=1/3, eps=.005, max.iter=1000, dfmax=p, delta=1e-8, gamma, group.multiplier=rep(1,J), verbose=FALSE, warn.conv=TRUE)
  {
    ## Check for errors
    family <- match.arg(family)
    penalty <- match.arg(penalty)
    if (alpha > 1 | alpha < 0) stop("alpha must be in [0,1]")
    if (length(group)!=ncol(X)) stop("group does not match X")
    if (is.null(colnames(X))) colnames(X) <- paste("V",1:ncol(X),sep="")
    if (delta <= 0) stop("Delta must be a positive number")
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

    ## Scale X
    n <- nrow(X)
    meanx <- apply(X,2,mean)
    normx <- sqrt(apply((t(X)-meanx)^2,1,sum))/sqrt(n)
    if (any(normx < 0.0001)) stop("X contains columns which are numerically constant.  If you are attempting to specify an intercept, please remove this column; an intercept is included automatically.")
    XX <- cbind(1,scale(X,meanx,normx))
    group <- c(0,group)
    colnames(XX)[1] <- "(Intercept)"
    p <- ncol(XX)

    ## Setup lambda
    if (missing(lambda))
      {
        lambda1 <- setupLambda(XX,y,group,family,penalty,lambda.max,lambda.min,nlambda,gamma,group.multiplier)
        lambda <- lambda1/alpha
      }
    l <- length(lambda)

    ## Fit
    path <- .C("gpPathFit",double(p*l),integer(l),double(l),as.double(XX),as.double(y),as.integer(group),family,as.integer(n),as.integer(p),as.integer(J),as.integer(K),penalty,as.double(lambda*alpha),as.integer(l),as.double(eps),as.integer(max.iter),as.integer(verbose),as.double(delta),as.double(gamma),as.double(tau),as.double((1-alpha)*lambda),as.double(group.multiplier),as.integer(dfmax),as.integer(warn.conv))

    ## Eliminate saturated lambda values, if any
    b <- matrix(path[[1]],nrow=p,byrow=T)
    ind <- !is.na(b[p,])
    b <- b[,ind,drop=FALSE]
    iter <- path[[2]][ind]
    lambda <- lambda[ind]
    beta <- unstandardize(b,meanx,normx)
    df <- path[[3]][ind]
    iter <- path[[2]][ind]

    ## Names
    if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
    else varnames <- colnames(X)
    varnames <- c("(Intercept)",varnames)
    dimnames(beta) <- list(varnames,round(lambda,digits=4))

    val <- list(beta=beta,
                family=family,
                group=group,
                lambda=lambda,
                alpha=alpha,
                loss = calcL(cbind(1,X),y,beta,family),
                n = length(y),
                penalty=penalty,
                df=df,
                iter=iter)
    class(val) <- "grpreg"
    return(val)
  }
