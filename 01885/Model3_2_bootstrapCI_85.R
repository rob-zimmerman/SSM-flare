rm(list = ls())

library(tidyverse)
library(optimParallel)
library(tictoc)
library(mvtnorm)
library(pbivnorm)
library(foreach)
library(doParallel)
library(doRNG)
library(xtable)
source("HMM_Functions.R")

set.seed(1729)

####### Useful functions #######

# state-dependent distribution
P.M3 <- function(y, beta.1, beta.2) {
  p <- dpois(y[1], w*beta.1*exp(0.5*(a[pairs[,2]] + a[pairs[,2] + 1]))) * dpois(y[2], w*beta.2*exp(0.5*(b[pairs[,1]] + b[pairs[,1] + 1])))
  return(diag(p))
}

# modified "elegant" pair and unpair functions
pair <- function(i,j) {
  if (i != max(i,j)) {
    return(j^2 - 2*j + i + 1)
  } else {
    return(i^2 + j - i)
  }
}

unpair <- function(j) {
  if (j - floor(sqrt(j-1))^2 - 1 < floor(sqrt(j-1))) {
    return( c(j - floor(sqrt(j-1))^2, floor(sqrt(j-1)) + 1) )
  } else {
    return( c(floor(sqrt(j-1)) + 1, j - floor(sqrt(j-1))^2 - floor(sqrt(j-1))) )
  }
}
unpair <- Vectorize(unpair)

# to transform unconstrained parameters to their natural parameter space
partrans.M3 <- function(CI.t) {
  c(tanh(CI.t[1]),
    tanh(CI.t[2]),
    exp(CI.t[3]),
    exp(CI.t[4]),
    exp(CI.t[5]),
    exp(CI.t[6]),
    tanh(CI.t[7])
  )
}


####### Other useful vectors

w <- 50
TT <- 2027
m1 <- 40
m2 <- 40
m <- m1*m2
a <- seq(from=-1.25, to=2.65, length.out=(m1+1))
b <- seq(from=-1.75, to=3.6, length.out=(m2+1))
a.mid <- a[-length(a)] + diff(a)/2
b.mid <- b[-length(b)] + diff(b)/2

pairran <- rep(0, times=m) # maps 1:m to pair(1:m1, 1:m2)
k <- 0
for (i in 1:m1) {
  for (j in 1:m2) {
    k <- k+1
    pairran[k] <- pair(i,j)
  }
}

unpairran <- rep(0, times=max(pairran)) # maps pair(1:m1, 1:m2) to 1:m
for (i in 1:length(pairran)) {
  unpairran[i] <- max(which(pairran == i),0)
}

ktopair <- function(k) {unpair(pairran[k])}

pairtok <- function(i,j) {unpairran[pair(i,j)]}

pairs <- expand.grid(j1 = 1:m1, j2 = 1:m2)
pairinv <- apply(pairs, MARGIN=1, FUN=function(x) pairtok(x[1],x[2]))




pars.hat <- c(0.97694799, 0.97455621, 0.09867089, 0.15681079, 0.18535641, 0.05872134, 0.99999990)
phi.1.hat <- pars.hat[1]
phi.2.hat <- pars.hat[2]
sigma.1.hat <- pars.hat[3]
sigma.2.hat <- pars.hat[4]
beta.1.hat <- pars.hat[5]
beta.2.hat <- pars.hat[6]
rho.hat <- pars.hat[7]

llk.M3 <- function(phi.1.t, phi.2.t, lsigma.1, lsigma.2, lbeta.1, lbeta.2, rho.t, Y)
{
  phi.1 <- tanh(phi.1.t)
  phi.2 <- tanh(phi.2.t)
  sigma.1 <- exp(lsigma.1)
  sigma.2 <- exp(lsigma.2)
  beta.1 <- exp(lbeta.1)
  beta.2 <- exp(lbeta.2)
  rho <- tanh(rho.t)
  
  Phi <- matrix(c(phi.1, 0,
                  0, phi.2), nrow=2, ncol=2, byrow=T)
  
  S <- matrix(c(sigma.1^2, rho*sigma.1*sigma.2,
                rho*sigma.1*sigma.2, sigma.2^2),
              nrow=2, ncol=2, byrow = T)
  
  L <- matrix(solve(diag(4) - Phi%x%Phi) %*% c(S), ncol=2, nrow=2)
  
  delta <- rep(0, times=m)
  for (j1 in 1:m1) {
    for (j2 in 1:m2) {
      J <- pairtok(j1, j2)
      delta[J] <- mvtnorm::pmvnorm(lower=c(a[j1], b[j2]), upper=c(a[j1+1], b[j2+1]), sigma=L)
    }
  }
  delta <- pmax(delta, 0)
  
  Gam <- matrix(0L, ncol=m1*m2, nrow=m1*m2)
  
  for (j in 1:(m1*m2)) {
    J <- pairtok(pairs[j,1], pairs[j,2])
    
    bstar <- 0.5*c(a[pairs[j,1]] + a[pairs[j,1]+1], b[pairs[j,2]] + b[pairs[j,2]+1]) # midpoint of box
    muJ <- as.vector(Phi%*%bstar)
    
    Gam[J,] <- pmax((pbivnorm(x=(a[pairs[,1] + 1] - muJ[1])/sigma.1, y=(b[pairs[,2] + 1] - muJ[2])/sigma.2, rho=rho) - 
                       pbivnorm(x=(a[pairs[,1]] - muJ[1])/sigma.1, y=(b[pairs[,2] + 1] - muJ[2])/sigma.2, rho=rho) -
                       pbivnorm(x=(a[pairs[,1] + 1] - muJ[1])/sigma.1, y=(b[pairs[,2]] - muJ[2])/sigma.2, rho=rho) + 
                       pbivnorm(x=(a[pairs[,1]] - muJ[1])/sigma.1, y=(b[pairs[,2]] - muJ[2])/sigma.2, rho=rho)),0)[pairinv]
  }
  
  # forward algorithm
  ph <- delta
  lscale <- 0
  for (t in 1:TT) {
    v <- ph %*% Gam %*% P.M3(Y[t,], beta.1, beta.2)
    u <- sum(v)
    lscale <- lscale + log(u)
    ph <- v/u
  }
  
  return(asFinite(ifelse(is.nan(lscale) || is.infinite(lscale), -1e300, lscale)))
}




####### Parametric bootstrap CIs and bias estimates #######

gen_data_and_est <- function(phi1, phi2, sig1, sig2, bet1, bet2, rho) {
  
  Phi <- matrix(c(phi1, 0,
                  0, phi2), nrow=2, ncol=2, byrow=T)
  
  S <- matrix(c(sig1^2, rho*sig1*sig2,
                rho*sig1*sig2, sig2^2),
              nrow=2, ncol=2, byrow = T)
  
  L <- matrix(solve(diag(4) - Phi%x%Phi) %*% c(S), ncol=2, nrow=2)
  
  
  # sample new underlying MC from the true model
  X.r <- matrix(0L, nrow=2, ncol=TT)
  X.r[,1] <- t(rmvnorm(n=1, sigma=L))
  for (t in 2:TT) {
    X.r[,t] <- Phi %*% X.r[,t-1] + t(rmvnorm(n=1, sigma=S))
  }
  
  # sample observed data from the true model
  Y.r <- matrix(0L, nrow=TT, ncol=2)
  for (t in 1:TT) {
    Y.r[t,] <- c(rpois(n=1, lambda=w*bet1*exp(X.r[1,t])), rpois(n=1, lambda=w*bet2*exp(X.r[2,t])))
  }
  
  pars.r <- optim(par=c(2.23109, 2.17715, -2.32408, -1.85088, -1.68772, -2.83894, 8.39999),
                  fn = function(p) llk.M3(phi.1.t=p[1],
                                          phi.2.t=p[2],
                                          lsigma.1=p[3],
                                          lsigma.2=p[4],
                                          lbeta.1=p[5],
                                          lbeta.2=p[6],
                                          rho.t=p[7],
                                          Y=Y.r),
                  method="L-BFGS-B",
                  control=list(fnscale=-1, trace=5),
                  hessian=FALSE
  )
  
  return(pars.r)
}

B <- 100 # number of bootstrap resamples

cl <- makePSOCKcluster(50, outfile="parlog.txt") # set the number of processor cores
registerDoParallel(cl)
registerDoRNG(1729, once=FALSE)

tic()
BS <- foreach(boot = 1:B, .errorhandling = "pass", .packages=c("mvtnorm", "pbivnorm"), .options.RNG = 1729) %dopar% {
  
  sink("log.txt", append = TRUE)
  cat(paste("Starting bootstrap sample", boot, "of", B, "\n"))
  sink()
  
  pars.hat <- c(0.97694799, 0.97455621, 0.09867089, 0.15681079, 0.18535641, 0.05872134, 0.99999990)
  phi.1.hat <- pars.hat[1]
  phi.2.hat <- pars.hat[2]
  sigma.1.hat <- pars.hat[3]
  sigma.2.hat <- pars.hat[4]
  beta.1.hat <- pars.hat[5]
  beta.2.hat <- pars.hat[6]
  rho.hat <- pars.hat[7]
  
  w <- 50
  TT <- 2027
  m1 <- 40
  m2 <- 40
  m <- m1*m2
  a <- seq(from=-1.25, to=2.65, length.out=(m1+1))
  b <- seq(from=-1.75, to=3.6, length.out=(m2+1))
  a.mid <- a[-length(a)] + diff(a)/2
  b.mid <- b[-length(b)] + diff(b)/2
  
  MLE.b <- gen_data_and_est(phi.1.hat, phi.2.hat, sigma.1.hat, sigma.2.hat, beta.1.hat, beta.2.hat, rho.hat)

  sink("log.txt", append = TRUE)
  cat(paste("Finished bootstrap sample", boot, "of", B, "\n"))
  sink()
  
  MLE.b
  
}
toc()

stopImplicitCluster()
registerDoSEQ()
toc()

gc()

save.image("Model3_2_bootstrapCI_85.Rdata")



bias.BS <- colMeans(t(sapply(BS, function(x) partrans.M3(x$par)))) - pars.hat # bootstrap estimate of bias
sd.BS <- sqrt(diag(cov(t(sapply(BS, function(x) partrans.M3(x$par)))))) # bootstrap estimate of SEs

pars.hat.cor <- pars.hat - bias.BS # bias-corrected MLE
CI.L.cor <- pars.hat.cor - 1.96*sd.BS
CI.U.cor <- pars.hat.cor + 1.96*sd.BS

tab_nums.M3 <- data.frame(Estimate = pars.hat.cor,
                          Std_Err = sd.BS,
                          CI.L = CI.L.cor,
                          CI.U = CI.U.cor)

colnames(tab_nums.M3) <- c("Estimate",
                           "Standard Error",
                           "CI (Lower)",
                           "CI (Upper)")

rownames(tab_nums.M3) <- c("$\\phi_{1}$",
                           "$\\phi_{2}$",
                           "$\\sigma_{1}$",
                           "$\\sigma_{2}$",
                           "$\\beta_{1}$",
                           "$\\beta_{2}$",
                           "$\\rho$")

print(xtable(tab_nums.M3, type = "latex", digits = 6, display = rep("f", times = 5)), booktabs = TRUE, sanitize.text.function = function(x) x, math.style.exponents = F)
