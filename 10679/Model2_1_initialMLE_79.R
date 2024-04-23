rm(list = ls())

library(tidyverse)
library(optimParallel)
library(FITSio)
library(tictoc)
library(foreach)
library(doParallel)
library(doRNG)
source("HMM_Functions.R")

set.seed(1729)

####### Extract data #######

EVLac79.src <- readFITS("EVLac_10679_src_dispersed_time_energy_wvl.fits", hdu=1)
EV.dat <- EVLac79.src$col

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

# to transform unconstrained parameters to their natural parameter space
partrans.M2 <- function(parvec.t) {
  c(tanh(parvec.t[1]),
    exp(parvec.t[2]),
    exp(parvec.t[3]),
    exp(parvec.t[4]),
    exp(parvec.t[5])
  )
}


####### Fit the model #######

m <- 40
cmin <- -1.5
cmax <- 2
c <- seq(from=cmin, to=cmax, length.out=(m+1))
c.mid <- c[-length(c)] + diff(c)/2

llk.M2 <- function(phi.1.t, lsigma.1, lsigma.2, lbeta.1, lbeta.2)
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
    v <- ph %*% Gam %*% P.M2(Y[t,], sigma.1, sigma.2, beta.1, beta.2)
    u <- sum(v)
    lscale <- lscale + log(u)
    ph <- v/u
  }
  
  return(asFinite(lscale))
}

cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t <- optimParallel(par=c(2.5, -2.4, -2.0, -1.9, -2.7),
                               fn = function(p) llk.M2(phi.1.t=p[1],
                                                       lsigma.1=p[2],
                                                       lsigma.2=p[3],
                                                       lbeta.1=p[4],
                                                       lbeta.2=p[5]),
                               method="L-BFGS-B",
                               control=list(fnscale=-1, trace=5),
                               hessian=FALSE
)

stopCluster(cl)

LL <- pars.t$value # log-likelihood
AIC <- 2*length(pars.t$par) - 2*LL  # AIC

# extract MLEs on correct scale
pars.t.hat <- pars.t$par
pars.hat <- partrans.M2(pars.t.hat)
phi.1.hat <- pars.hat[1]
sigma.1.hat <- pars.hat[2]
sigma.2.hat <- pars.hat[3]
beta.1.hat <- pars.hat[4]
beta.2.hat <- pars.hat[5]


######## Make state predictions to verify correctness of essential domain ######

delta.hat <- rep(0, m)
for (j in 1:m) {
  delta.hat[j] <- pnorm(c[j+1], mean=0, sd=sqrt(sigma.1.hat^2/(1-phi.1.hat^2))) - pnorm(c[j], mean=0, sd=sqrt(sigma.1.hat^2/(1-phi.1.hat^2)))
}
delta.hat <- delta.hat/sum(delta.hat)

Gam.hat <- matrix(0L, nrow=m, ncol=m)
for (i in 1:m) {
  for (j in 1:m) {
    Gam.hat[i,j] <- pnorm(c[j+1], mean=phi.1.hat*c.mid[i], sd=sigma.1.hat) - pnorm(c[j], mean=phi.1.hat*c.mid[i], sd=sigma.1.hat)
  }
}
Gam.hat <- Gam.hat/rowSums(Gam.hat)

P.M2.hat <- function(y) {asFinite(diag(P.M2(y, sigma.1.hat, sigma.2.hat, beta.1.hat, beta.2.hat)))}

X1.hat <- c.mid[local_dec(init=delta.hat, tpm=Gam.hat, data=t(Y), Fx=P.M2.hat)$preds] # predicted X1t's

ggplot(data.frame(x=X1.hat), aes(x=x)) + 
  geom_histogram(bins=30, color="black", fill="lightblue") +
  geom_vline(aes(xintercept=cmin), color="red", linetype="dashed") +
  geom_vline(aes(xintercept=cmax), color="red", linetype="dashed") +
  ylab("Count") + 
  xlab(bquote(hat(X)["t,1"])) + 
  ggtitle("EV Lac Model 2 - Histogram of Predicted States (Soft Band)")
