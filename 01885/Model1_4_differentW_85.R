rm(list = ls())

library(tidyverse)
library(optimParallel)
library(FITSio)
library(tictoc)
library(foreach)
library(doParallel)
library(doRNG)
library(xtable)
source("HMM_Functions.R")

set.seed(1729)


####### Useful functions #######

# state-dependent distribution - w is absorbed into the betas here but extracted from the estimates after
P.M1 <- function(y, beta.1, beta.2) {
  p <- rep(0, times=m)
  for (i in 1:m) {
    p[i] <- dpois(y[1], beta.1*exp(c.mid[i]))*dpois(y[2], beta.2*exp(c.mid[i]))
  }
  return(diag(p))
}

# to transform unconstrained parameters to their natural parameter space
partrans.M1 <- function(parvec.t) {
  c(tanh(parvec.t[1]),
    exp(parvec.t[2]),
    exp(parvec.t[3])/w,
    exp(parvec.t[4])/w
  )
}


####### Fit the model #######

m <- 40
cmin <- -2.5
cmax <- 2.75
c <- seq(from=cmin, to=cmax, length.out=(m+1))
c.mid <- c[-length(c)] + diff(c)/2

llk.M1 <- function(phi.1.t, lsigma.1, lbeta.1, lbeta.2)
{
  phi.1 <- tanh(phi.1.t)
  sigma.1 <- exp(lsigma.1)
  beta.1 <- exp(lbeta.1)
  beta.2 <- exp(lbeta.2)
  
  delta <- rep(0, m)
  for (j in 1:m) {
    delta[j] <- pnorm(c[j+1], mean=0, sd=sqrt(sigma.1^2/(1-phi.1^2))) - pnorm(c[j], mean=0, sd=sqrt(sigma.1^2/(1-phi.1^2)))
  }

  Gam <- matrix(0L, nrow=m, ncol=m)
  for (i in 1:m) {
    for (j in 1:m) {
      Gam[i,j] <- pnorm(c[j+1], mean=phi.1*c.mid[i], sd=sigma.1) - pnorm(c[j], mean=phi.1*c.mid[i], sd=sigma.1)
    }
  }

  # forward algorithm
  ph <- delta
  lscale <- 0
  for (t in 1:TT) {
    v <- ph %*% Gam %*% P.M1(Y[t,], beta.1, beta.2)
    u <- sum(v)
    lscale <- lscale + log(u)
    ph <- v/u
  }
  
  return(asFinite(ifelse(is.nan(lscale), -1e300, lscale)))
}





####### w = 25 #######

EVLac85.src <- readFITS("EVLac_01885_src_dispersed_time_energy_wvl.fits", hdu=1)
EV.dat <- EVLac85.src$col

Y <- matrix(0L, ncol=3, nrow=length(unlist(EV.dat[1]))) # create matrix for observations
colnames(Y) <- c("TIME", "ENERGY", "WAVELENGTH")
Y[,1] <- unlist(EV.dat[1]) # time (in spacecraft clock seconds)
Y.t_min <- min(Y[,1]) # first observed time
Y[,1] <- Y[,1] - Y.t_min # start time at 0
Y[,2] <- unlist(EV.dat[2]) # energy (in keV)
Y[,3] <- unlist(EV.dat[3]) # wavelength (signed, in Angstroms)

Y <- as_tibble(Y)

w <- 25 # time window (in seconds) for each HMM observation
b <- 1.5 # energy threshold to split the observations

Y$t <- floor(Y$TIME/w) + 1 # discretize time
Y <- Y %>% group_by(t) %>% summarise(LOW = sum(ENERGY < b), HIGH = sum(ENERGY >= b))
TT <- max(Y$t) # number of observations after discretization

Y <- merge(
  data.frame(t = with(Y, seq(min(t), max(t), 1))),
  Y,
  all = TRUE
)
Y <- replace(Y, is.na(Y), 0) # times where no photons detected -> counts of zeros
Y <- as.matrix(dplyr::select(Y, c("LOW", "HIGH")))


cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t.25 <- optimParallel(par=c(2.20, -2.14, 1.68, -2.57),
                               fn = function(p) llk.M1(phi.1.t=p[1],
                                                       lsigma.1=p[2],
                                                       lbeta.1=p[3],
                                                       lbeta.2=p[4]),
                               method="L-BFGS-B",
                               control=list(fnscale=-1, trace=6, REPORT=1),
                               hessian=FALSE
)

stopCluster(cl)

# extract MLEs on correct scale
pars.hat.25 <- partrans.M1(pars.t.25$par)





####### w = 50 #######

EVLac85.src <- readFITS("EVLac_01885_src_dispersed_time_energy_wvl.fits", hdu=1)
EV.dat <- EVLac85.src$col

Y <- matrix(0L, ncol=3, nrow=length(unlist(EV.dat[1]))) # create matrix for observations
colnames(Y) <- c("TIME", "ENERGY", "WAVELENGTH")
Y[,1] <- unlist(EV.dat[1]) # time (in spacecraft clock seconds)
Y.t_min <- min(Y[,1]) # first observed time
Y[,1] <- Y[,1] - Y.t_min # start time at 0
Y[,2] <- unlist(EV.dat[2]) # energy (in keV)
Y[,3] <- unlist(EV.dat[3]) # wavelength (signed, in Angstroms)

Y <- as_tibble(Y)

w <- 50 # time window (in seconds) for each HMM observation
b <- 1.5 # energy threshold to split the observations

Y$t <- floor(Y$TIME/w) + 1 # discretize time
Y <- Y %>% group_by(t) %>% summarise(LOW = sum(ENERGY < b), HIGH = sum(ENERGY >= b))
TT <- max(Y$t) # number of observations after discretization

Y <- merge(
  data.frame(t = with(Y, seq(min(t), max(t), 1))),
  Y,
  all = TRUE
)
Y <- replace(Y, is.na(Y), 0) # times where no photons detected -> counts of zeros
Y <- as.matrix(dplyr::select(Y, c("LOW", "HIGH")))


cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t.50 <- optimParallel(par=c(2.20, -2.14, 1.68, -2.57),
                           fn = function(p) llk.M1(phi.1.t=p[1],
                                                   lsigma.1=p[2],
                                                   lbeta.1=p[3],
                                                   lbeta.2=p[4]),
                           method="L-BFGS-B",
                           control=list(fnscale=-1, trace=6, REPORT=1),
                           hessian=FALSE
)

stopCluster(cl)

# extract MLEs on correct scale
pars.hat.50 <- partrans.M1(pars.t.50$par)





####### w = 75 #######

EVLac85.src <- readFITS("EVLac_01885_src_dispersed_time_energy_wvl.fits", hdu=1)
EV.dat <- EVLac85.src$col

Y <- matrix(0L, ncol=3, nrow=length(unlist(EV.dat[1]))) # create matrix for observations
colnames(Y) <- c("TIME", "ENERGY", "WAVELENGTH")
Y[,1] <- unlist(EV.dat[1]) # time (in spacecraft clock seconds)
Y.t_min <- min(Y[,1]) # first observed time
Y[,1] <- Y[,1] - Y.t_min # start time at 0
Y[,2] <- unlist(EV.dat[2]) # energy (in keV)
Y[,3] <- unlist(EV.dat[3]) # wavelength (signed, in Angstroms)

Y <- as_tibble(Y)

w <- 75 # time window (in seconds) for each HMM observation
b <- 1.5 # energy threshold to split the observations

Y$t <- floor(Y$TIME/w) + 1 # discretize time
Y <- Y %>% group_by(t) %>% summarise(LOW = sum(ENERGY < b), HIGH = sum(ENERGY >= b))
TT <- max(Y$t) # number of observations after discretization

Y <- merge(
  data.frame(t = with(Y, seq(min(t), max(t), 1))),
  Y,
  all = TRUE
)
Y <- replace(Y, is.na(Y), 0) # times where no photons detected -> counts of zeros
Y <- as.matrix(dplyr::select(Y, c("LOW", "HIGH")))


cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t.75 <- optimParallel(par=c(2.20, -2.14, 1.68, -2.57),
                           fn = function(p) llk.M1(phi.1.t=p[1],
                                                   lsigma.1=p[2],
                                                   lbeta.1=p[3],
                                                   lbeta.2=p[4]),
                           method="L-BFGS-B",
                           control=list(fnscale=-1, trace=6, REPORT=1),
                           hessian=FALSE
)

stopCluster(cl)

# extract MLEs on correct scale
pars.hat.75 <- partrans.M1(pars.t.75$par)





####### w = 100 #######

EVLac85.src <- readFITS("EVLac_01885_src_dispersed_time_energy_wvl.fits", hdu=1)
EV.dat <- EVLac85.src$col

Y <- matrix(0L, ncol=3, nrow=length(unlist(EV.dat[1]))) # create matrix for observations
colnames(Y) <- c("TIME", "ENERGY", "WAVELENGTH")
Y[,1] <- unlist(EV.dat[1]) # time (in spacecraft clock seconds)
Y.t_min <- min(Y[,1]) # first observed time
Y[,1] <- Y[,1] - Y.t_min # start time at 0
Y[,2] <- unlist(EV.dat[2]) # energy (in keV)
Y[,3] <- unlist(EV.dat[3]) # wavelength (signed, in Angstroms)

Y <- as_tibble(Y)

w <- 100 # time window (in seconds) for each HMM observation
b <- 1.5 # energy threshold to split the observations

Y$t <- floor(Y$TIME/w) + 1 # discretize time
Y <- Y %>% group_by(t) %>% summarise(LOW = sum(ENERGY < b), HIGH = sum(ENERGY >= b))
TT <- max(Y$t) # number of observations after discretization

Y <- merge(
  data.frame(t = with(Y, seq(min(t), max(t), 1))),
  Y,
  all = TRUE
)
Y <- replace(Y, is.na(Y), 0) # times where no photons detected -> counts of zeros
Y <- as.matrix(dplyr::select(Y, c("LOW", "HIGH")))


cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t.100 <- optimParallel(par=c(2.20, -2.14, 1.68, -2.57),
                           fn = function(p) llk.M1(phi.1.t=p[1],
                                                   lsigma.1=p[2],
                                                   lbeta.1=p[3],
                                                   lbeta.2=p[4]),
                           method="L-BFGS-B",
                           control=list(fnscale=-1, trace=6, REPORT=1),
                           hessian=FALSE
)

stopCluster(cl)

# extract MLEs on correct scale
pars.hat.100 <- partrans.M1(pars.t.100$par)




######## Make a nice table

tab_nums.M1 <- data.frame(w25 = pars.hat.25, w50 = pars.hat.50, w75 = pars.hat.75, w100 = pars.hat.100)

colnames(tab_nums.M1) <- c("$w = 25$",
                           "$w = 50$",
                           "$w = 75$",
                           "$w = 100$")

rownames(tab_nums.M1) <- c("$\\phi_{1}$",
                           "$\\sigma_{1}$",
                           "$\\beta_{1}$",
                           "$\\beta_{2}$")

print(xtable(tab_nums.M1, type = "latex", digits = 6, display = rep("f", times = 5)), booktabs = TRUE, sanitize.text.function = function(x) x, math.style.exponents = F)
