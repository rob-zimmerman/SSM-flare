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

# state-dependent distribution
P.M3 <- function(y, beta.1, beta.2) {
  p <- dpois(y[1], beta.1*exp(0.5*(a[pairs[,2]] + a[pairs[,2] + 1]))) * dpois(y[2], beta.2*exp(0.5*(b[pairs[,1]] + b[pairs[,1] + 1])))
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
    exp(CI.t[5])/w,
    exp(CI.t[6])/w,
    tanh(CI.t[7])
  )
}


####### Fit the model #######

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


####### Fit the model #######

llk.M3 <- function(phi.1.t, phi.2.t, lsigma.1, lsigma.2, lbeta.1, lbeta.2, rho.t)
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

m1 <- 40
m2 <- 40
m <- m1*m2
a <- seq(from=-1.25, to=2.65, length.out=(m1+1))
b <- seq(from=-1.75, to=3.6, length.out=(m2+1))
a.mid <- a[-length(a)] + diff(a)/2
b.mid <- b[-length(b)] + diff(b)/2

cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t.25 <- optimParallel(par=c(2.23109, 2.17715, -2.32408, -1.85088, -1.68772, -2.83894, 8.39999),
                           fn = function(p) llk.M3(phi.1.t=p[1],
                                                   phi.2.t=p[2],
                                                   lsigma.1=p[3],
                                                   lsigma.2=p[4],
                                                   lbeta.1=p[5],
                                                   lbeta.2=p[6],
                                                   rho.t=p[7]),
                           method="L-BFGS-B",
                           control=list(fnscale=-1, trace=6),
                           hessian=FALSE
)

stopCluster(cl)

# extract MLEs on correct scale
pars.hat.25 <- partrans.M3(pars.t.25$par)
save.image("Model3_4_differentW_85.Rdata")





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

m1 <- 40
m2 <- 40
m <- m1*m2
a <- seq(from=-1.25, to=2.65, length.out=(m1+1))
b <- seq(from=-1.75, to=3.6, length.out=(m2+1))
a.mid <- a[-length(a)] + diff(a)/2
b.mid <- b[-length(b)] + diff(b)/2

cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t.50 <- optimParallel(par=c(2.23109, 2.17715, -2.32408, -1.85088, -1.68772, -2.83894, 8.39999),
                           fn = function(p) llk.M3(phi.1.t=p[1],
                                                   phi.2.t=p[2],
                                                   lsigma.1=p[3],
                                                   lsigma.2=p[4],
                                                   lbeta.1=p[5],
                                                   lbeta.2=p[6],
                                                   rho.t=p[7]),
                           method="L-BFGS-B",
                           control=list(fnscale=-1, trace=6),
                           hessian=FALSE
)

stopCluster(cl)

# extract MLEs on correct scale
pars.hat.50 <- partrans.M3(pars.t.50$par)
save.image("Model3_4_differentW_85.Rdata")




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

m1 <- 40
m2 <- 40
m <- m1*m2
a <- seq(from=-1.25, to=2.65, length.out=(m1+1))
b <- seq(from=-1.75, to=3.6, length.out=(m2+1))
a.mid <- a[-length(a)] + diff(a)/2
b.mid <- b[-length(b)] + diff(b)/2

cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t.75 <- optimParallel(par=c(2.23109, 2.17715, -2.32408, -1.85088, 2.23, 1.08, 8.39999),
                           fn = function(p) llk.M3(phi.1.t=p[1],
                                                   phi.2.t=p[2],
                                                   lsigma.1=p[3],
                                                   lsigma.2=p[4],
                                                   lbeta.1=p[5],
                                                   lbeta.2=p[6],
                                                   rho.t=p[7]),
                           method="L-BFGS-B",
                           control=list(fnscale=-1, trace=6),
                           hessian=FALSE
)

stopCluster(cl)

# extract MLEs on correct scale
pars.hat.75 <- partrans.M3(pars.t.75$par)
save.image("Model3_4_differentW_85.Rdata")




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

m1 <- 40
m2 <- 40
m <- m1*m2
a <- seq(from=-1.25, to=2.65, length.out=(m1+1))
b <- seq(from=-1.75, to=3.6, length.out=(m2+1))
a.mid <- a[-length(a)] + diff(a)/2
b.mid <- b[-length(b)] + diff(b)/2

cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t.100 <- optimParallel(par=c(2.23109, 2.17715, -2.32408, -1.85088, 2.23, 1.08, 8.39999),
                            fn = function(p) llk.M3(phi.1.t=p[1],
                                                    phi.2.t=p[2],
                                                    lsigma.1=p[3],
                                                    lsigma.2=p[4],
                                                    lbeta.1=p[5],
                                                    lbeta.2=p[6],
                                                    rho.t=p[7]),
                            method="L-BFGS-B",
                            control=list(fnscale=-1, trace=6),
                            hessian=FALSE
)

stopCluster(cl)

# extract MLEs on correct scale
pars.hat.100 <- partrans.M3(pars.t.100$par)
save.image("Model3_4_differentW_85.Rdata")



######## Make a nice table

tab_nums.M1 <- data.frame(w25 = pars.hat.25, w50 = pars.hat.50, w75 = pars.hat.75, w100 = pars.hat.100)

colnames(tab_nums.M1) <- c("$w = 25$",
                           "$w = 50$",
                           "$w = 75$",
                           "$w = 100$")

rownames(tab_nums.M1) <- c("$\\phi_{1}$",
                           "$\\phi_{2}$",
                           "$\\sigma_{1}$",
                           "$\\sigma_{2}$",
                           "$\\beta_{1}$",
                           "$\\beta_{2}$",
                           "$\\rho$")

print(xtable(tab_nums.M1, type = "latex", digits = 4, display = rep("f", times = 5)), booktabs = TRUE, sanitize.text.function = function(x) x, math.style.exponents = F)
