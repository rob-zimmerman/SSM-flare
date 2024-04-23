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



####### Useful functions #######

# state-dependent distribution
P.M1 <- function(y, beta.1, beta.2) {
  p <- rep(0, times=m)
  for (i in 1:m) {
    p[i] <- dpois(y[1], w*beta.1*exp(c.mid[i]))*dpois(y[2], w*beta.2*exp(c.mid[i]))
  }
  return(diag(p))
}

# to transform unconstrained parameters to their natural parameter space
partrans.M1 <- function(parvec.t) {
  c(tanh(parvec.t[1]),
    exp(parvec.t[2]),
    exp(parvec.t[3]),
    exp(parvec.t[4])
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

cl <- makeCluster(11, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t <- optimParallel(par=c(2.20, -2.14, 1.68, -2.57),
                               fn = function(p) llk.M1(phi.1.t=p[1],
                                                       lsigma.1=p[2],
                                                       lbeta.1=p[3],
                                                       lbeta.2=p[4]),
                               method="L-BFGS-B",
                               control=list(fnscale=-1, trace=6, REPORT=1),
                               hessian=FALSE
)

stopCluster(cl)

LL <- pars.t$value # log-likelihood
AIC <- 2*length(pars.t$par) - 2*LL  # AIC

# extract MLEs on correct scale
pars.t.hat <- pars.t$par
pars.hat <- partrans.M1(pars.t.hat)
phi.1.hat <- pars.hat[1]
sigma.1.hat <- pars.hat[2]
beta.1.hat <- pars.hat[3]
beta.2.hat <- pars.hat[4]


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

P.M1.hat <- function(y) {asFinite(diag(P.M1(y, beta.1.hat, beta.2.hat)))}

X.hat <- c.mid[local_dec(init=delta.hat, tpm=Gam.hat, data=t(Y), Fx=P.M1.hat)$preds] # predicted X1t's

ggplot(data.frame(x=X.hat), aes(x=x)) + 
  geom_histogram(bins=30, color="black", fill="lightblue") +
  geom_vline(aes(xintercept=cmin), color="red", linetype="dashed") +
  geom_vline(aes(xintercept=cmax), color="red", linetype="dashed") +
  ylab("Count") + 
  xlab(bquote(hat(X)["t"])) + 
  ggtitle("EV Lac Model 2 - Histogram of Predicted States (Soft Band)")
