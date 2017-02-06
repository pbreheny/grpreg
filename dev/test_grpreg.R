# rm(list = ls())
# require(grpreg)
# require(microbenchmark)
# 
# data(Birthwt)
# X <- Birthwt$X
# group <- Birthwt$group
# 
# ## Linear regression
# y <- Birthwt$bwt
# fit <- grpreg(X, y, group, penalty="grLasso", screen = "None")
# plot(fit, main = "No screen")
# 
# fit.ssr <- grpreg(X, y, group, penalty = "grLasso", screen = "SSR")
# plot(fit.ssr, main = "SSR")
# 
# fit.sedpp <- grpreg(X, y, group, peanlty = "grLasso", screen = "SEDPP")
# plot(fit.sedpp, main = "SEDPP")
# 
# fit.ssr.bedpp <- grpreg(X, y, group, penalty = 'grLasso', screen = 'SSR-BEDPP')
# plot(fit.ssr.bedpp, main = "SSR-BEDPP")
# 
# 
# fit.nac <- grpreg(X, y, group, penalty="grLasso", screen = "No-Active")
# plot(fit, main = "NAC")
# 
# fit.ssr <- grpreg(X, y, group, penalty = "grLasso", screen = "SSR")
# plot(fit.ssr, main = "SSR")
# 
# fit.sedpp <- grpreg(X, y, group, peanlty = "grLasso", screen = "SEDPP")
# plot(fit.sedpp, main = "SEDPP")
# 
# fit.ssr.bedpp <- grpreg(X, y, group, penalty = 'grLasso', screen = 'SSR-BEDPP')
# plot(fit.ssr.bedpp, main = "SSR-BEDPP")
# 
# 
# all.equal(fit$loss, fit.ssr$loss)
# all.equal(fit$beta, fit.ssr$beta)
# all.equal(fit$loss, fit.sedpp$loss)
# all.equal(fit$beta, fit.ssr.bedpp$beta)
# all.equal(fit$loss, fit.ssr.bedpp$loss)
# all.equal(fit$beta, fit.ssr.bedpp$beta)

# simulated data
rm(list = ls())
require(grpreg)
require(microbenchmark)
setwd("~/GitHub/grpreg_experiment/0116_debug/")

n <- 200
grp <- 500
grp.size <- 10
true.grp <- 10
p <- grp * grp.size
eps <- 0.001
lambda.min <- 0.1
rep <- 2
date <- Sys.Date()
methods <- c("None", "No-Active", "SSR-No-Active", "SSR", "SEDPP-No-Active",
             "SEDPP", "SSR-BEDPP-No-Active", "SSR-BEDPP")
methods.name <- c("NS-AC", "NS-NAC", "SSR-NAC", "SSR-AC", "SEDPP-NAC", "SEDPP-AC", "SSR-BEDPP-NAC", "SSR-BEDPP-AC")

time.model <- matrix(NA, ncol = length(methods), nrow = rep)
time.all <- matrix(NA, ncol = length(methods), nrow = rep)
colnames(time.model) <- methods.name
colnames(time.all) <- methods.name

for (i in 1:rep) {
  print(i)
  set.seed(1234+i-1)
  X <- matrix(rnorm(n*p), ncol = p)
  beta.true <- runif(true.grp, -1, 1)
  beta <- c(rep(beta.true, each = grp.size), rep(0, p - true.grp * grp.size))
  y <- X %*% beta + 0.1 * rnorm(n)
  group <- rep(1:grp, each = grp.size)
  
  pdf(file = paste0(date, "_", i, "_path.pdf"))
  t <- system.time(
    fit1 <- grpreg(X, y, group = group, penalty="grLasso", screen = methods[1],
                  lambda.min = lambda.min, eps = eps, log.lambda = FALSE)
  )
  time.model[i, 1] <- fit1$time
  time.all[i, 1] <- as.numeric(t['elapsed'])
  plot(fit1, main = methods.name[1])

  sink(file = paste0(date, "_", i, "_output.txt"), append = FALSE)
  
  for (j in 2:length(methods)) {
    t <- system.time(
      fit <- grpreg(X, y, group = group, penalty="grLasso", screen = methods[j],
                    lambda.min = lambda.min, eps = eps, log.lambda = FALSE)
    )
    time.model[i, j] <- fit$time
    time.all[i, j] <- as.numeric(t['elapsed'])
    plot(fit, main = methods.name[j])
    
    cat("\n\n--------------------------------------\n\n")
    cat("Method: ", methods[j], "\n")
    cat("--------------------------------------\n\n")
    
    print(all.equal(fit1$loss, fit$loss))
    print(all.equal(fit1$beta, fit$beta))
    
    if (methods[j] %in% c("SSR", "SEDPP", "SSR-BEDPP", "SSR-No-Active", "SEDPP-No-Active", "SSR-BEDPP-No-Active")) {
      cat("# of rejections:", methods[j], "\n")
      print(fit$rejections)
    }

    if (methods[j] %in% c("SSR-BEDPP", "SSR-BEDPP-No-Active")) {
      cat("# of safe rejections:", methods[j], "\n")
      print(fit$safe_rejections)
    }
    
  }
  sink()
  dev.off()
  
  save(X, y, group, beta.true, beta, file = paste0(date, "_", i, "_path.RData"))
  rm(X, y, group, beta.true, beta)
  gc()
  
  
}
print(time.model)
print(colMeans(time.model))

print(time.all)
print(colMeans(time.all))





