rm(list = ls())

library(tidyverse)
library(optimParallel)
library(tictoc)
library(foreach)
library(doParallel)
library(doRNG)
library(xtable)

set.seed(1729)

####### Useful functions #######

# state-dependent distribution
P.M2 <- function(y, sigma.1, sigma.2, beta.1, beta.2) {
  p <- rep(0, times=m)
  for (i in 1:m) {
    cstar.1 <- c.mid[i]
    cstar.2 <- (sigma.2/sigma.1)*cstar.1
    p[i] <- dpois(y[1], w*beta.1*exp(cstar.1))*dpois(y[2], w*beta.2*exp(cstar.2))
  }
  return(diag(p))
}

# convert signed Infs to max-precision values (from 'copula' package)
asFinite <- function(x) {
  if (any(nifi <- !is.finite(x)))
    x[nifi] <- sign(x[nifi]) * .Machine$double.xmax
  x
}

# to transform unconstrained parameters to their natural parameter space
partrans.M2 <- function(parvec.t) {
  c(tanh(parvec.t[1]),
    exp(parvec.t[2]),
    exp(parvec.t[3]),
    exp(parvec.t[4]),
    exp(parvec.t[5])
  )
}




pars.hat <- c(0.98675254, 0.08977194, 0.13132740, 0.14747268, 0.06342463)
phi.1.hat <- pars.hat[1]
sigma.1.hat <- pars.hat[2]
sigma.2.hat <- pars.hat[3]
beta.1.hat <- pars.hat[4]
beta.2.hat <- pars.hat[5]

w <- 50
TT <- 1937
m <- 40
cmin <- -1.5
cmax <- 2
c <- seq(from=cmin, to=cmax, length.out=(m+1))
c.mid <- c[-length(c)] + diff(c)/2

llk.M2 <- function(phi.1.t, lsigma.1, lsigma.2, lbeta.1, lbeta.2, Y)
{
  phi.1 <- tanh(phi.1.t)
  sigma.1 <- exp(lsigma.1)
  sigma.2 <- exp(lsigma.2)
  beta.1 <- exp(lbeta.1)
  beta.2 <- exp(lbeta.2)
  
  delta <- rep(0, m)
  for (j in 1:m) {
    delta[j] <- pnorm(c[j+1], mean=0, sd=sqrt(sigma.1^2/(1-phi.1^2))) - pnorm(c[j], mean=0, sd=sqrt(sigma.1^2/(1-phi.1^2)))
  }
  delta <- delta/sum(delta)
  
  Gam <- matrix(0L, nrow=m, ncol=m)
  for (i in 1:m) {
    for (j in 1:m) {
      Gam[i,j] <- pnorm(c[j+1], mean=phi.1*c.mid[i], sd=sigma.1) - pnorm(c[j], mean=phi.1*c.mid[i], sd=sigma.1)
    }
  }
  Gam <- Gam/rowSums(Gam)
  
  # forward algorithm
  ph <- delta
  lscale <- 0
  for (t in 1:TT) {
    v <- ph %*% Gam %*% P.M2(Y[t,], sigma.1, sigma.2, beta.1, beta.2)
    u <- sum(v)
    lscale <- lscale + log(u)
    ph <- v/u
  }
  
  return(asFinite(lscale))
}




####### Parametric bootstrap CIs and bias estimates #######

gen_data_and_est <- function(phi1, sig1, sig2, bet1, bet2) {
  # sample new underlying MC from the true model
  X1.r <- 0*1:TT
  X1.r[1] <- rnorm(n=1, mean=0, sd=sig1/sqrt(1-phi1^2))
  for (t in 2:TT) {
    X1.r[t] <- phi1*X1.r[t-1] + rnorm(n=1, sd=sig1)
  }
  X2.r <- (sig2/sig1)*X1.r
  
  # sample observed data from the true model
  Y.r <- matrix(0L, nrow=TT, ncol=2)
  for (t in 1:TT) {
    Y.r[t,] <- c(rpois(n=1, lambda=w*bet1*exp(X1.r[t])), rpois(n=1, lambda=w*bet2*exp(X2.r[t])))
  }
  
  pars.r <- optim(par=c(2.5, -2.4, -2.0, -1.9, -2.7),
                  fn = function(p) llk.M2(phi.1.t=p[1],
                                          lsigma.1=p[2],
                                          lsigma.2=p[3],
                                          lbeta.1=p[4],
                                          lbeta.2=p[5],
                                          Y=Y.r),
                  method="L-BFGS-B",
                  control=list(fnscale=-1, trace=5),
                  hessian=FALSE
  )
  
  return(pars.r)
}

B <- 100 # number of bootstrap resamples

cl <- makePSOCKcluster(100) # set the number of processor cores
registerDoParallel(cl)
registerDoRNG(1729, once=FALSE)

tic()
BS <- foreach(b = 1:B, .errorhandling = "pass", .options.RNG = 1729) %dopar% {
  
  MLE.b <- gen_data_and_est(phi.1.hat, sigma.1.hat, sigma.2.hat, beta.1.hat, beta.2.hat)

  sink("log.txt", append = TRUE)
  cat(paste("Finished bootstrap sample", b, "of", B, "\n"))
  sink()
  
  MLE.b
  
}
toc()

stopImplicitCluster()
registerDoSEQ()
toc()

gc()

save.image("Model2_2_bootstrapCI_79.Rdata")



bias.BS <- colMeans(t(sapply(BS, function(x) partrans.M2(x$par)))) - pars.hat # bootstrap estimate of bias
sd.BS <- sqrt(diag(cov(t(sapply(BS, function(x) partrans.M2(x$par)))))) # bootstrap estimate of SEs

pars.hat.cor <- pars.hat - bias.BS # bias-corrected MLE
CI.L.cor <- pars.hat.cor - 1.96*sd.BS
CI.U.cor <- pars.hat.cor + 1.96*sd.BS

tab_nums.M2 <- data.frame(Estimate = pars.hat.cor,
                          Std_Err = sd.BS,
                          CI.L = CI.L.cor,
                          CI.U = CI.U.cor)

colnames(tab_nums.M2) <- c("Estimate",
                           "Standard Error",
                           "CI (Lower)",
                           "CI (Upper)")

rownames(tab_nums.M2) <- c("$\\phi_{1}$",
                           "$\\sigma_{1}$",
                           "$\\sigma_{2}$",
                           "$\\beta_{1}$",
                           "$\\beta_{2}$")

print(xtable(tab_nums.M2, type = "latex", digits = 6, display = rep("f", times = 5)), booktabs = TRUE, sanitize.text.function = function(x) x, math.style.exponents = F)
