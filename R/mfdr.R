mfdr <- function(fit, X) {

  # Initial checks
  if (!inherits(fit, 'grpreg')) stop('"fit" must be an grpreg object', call.=FALSE)
  if (!(fit$penalty == "grLasso" | fit$penalty == "grMCP")) stop('"mFDR" is only avaiable for "grLasso" and "grMCP" penalties', call.=FALSE)
  if (inherits(fit, "grpsurv") || fit$family == "binomial") {
    if (!("XG" %in% names(fit))) {
      stop("The argument 'returnX=TRUE' is needed to calculate mFDR for GLM/Cox models", call.=FALSE)
      if (missing(X)) {
        stop("X must be supplied in addition to the model fit in order to calculate mFDR for GLM/Cox models", call.=FALSE)
      }
    }
  }
  
  ## Setup
  S0 <- sum(fit$group.multiplier==0)
  S <- predict(fit, type="ngroups") - S0
  fit$gl <- as.vector(table(fit$group))
  
  # Call C functions
  if (inherits(fit, "grpsurv")) {
    #stop('"mFDR" has not yet been implemented for survival models', call.=FALSE)
    fit$XX <- fit$XG$X
    EF <- .Call("mfdr_cox", fit)
  } else {
    if (fit$family == "binomial") {
       fit$P <- matrix(nrow = length(fit$lambda), ncol = fit$n)
       for(i in 1:length(fit$lambda)){
         fit$P[i,] <- predict(fit, X = X, lambda = fit$lambda[i], type = "response")  
       }
       fit$XX <- fit$XG$X
      EF <- .Call("mfdr_binomial", fit)
    } else if (fit$family == "gaussian") {
      EF <- .Call("mfdr_gaussian", fit)
    }
  }
  
  # Calculate rate, return
  EF <- pmin(EF - S0, S)
  mFDR <- EF/S
  mFDR[S==0] <- 0
  df <- data.frame(EF=EF, S=S, mFDR=mFDR)
  rownames(df) <- lamNames(fit$lambda)
  structure(df, class=c("mfdr", "data.frame"))
    
}