rm(list = ls())

library(tidyverse)
library(optimParallel)
library(FITSio)
library(mvtnorm)
library(mnormt)
library(pbivnorm)
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

tic()
cl <- makeCluster(100, type="FORK") # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pars.t <- optimParallel(par=c(2.23109, 2.17715, -2.32408, -1.85088, -1.68772, -2.83894, 8.39999),
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
toc()

LL <- pars.t$value # log-likelihood
AIC <- 2*length(pars.t$par) - 2*LL  # AIC

# extract MLEs on correct scale
pars.t.hat <- pars.t$par
pars.hat <- partrans.M3(pars.t.hat)
phi.1.hat <- pars.hat[1]
phi.2.hat <- pars.hat[2]
sigma.1.hat <- pars.hat[3]
sigma.2.hat <- pars.hat[4]
beta.1.hat <- pars.hat[5]
beta.2.hat <- pars.hat[6]
rho.hat <- pars.hat[7]

save.image("Model3_1_initialMLE_85_fast.Rdata")


######## Make state predictions to verify correctness of essential domain ######

Phi.hat <- matrix(c(phi.1.hat, 0,
                0, phi.2.hat), nrow=2, ncol=2, byrow=T)

S.hat <- matrix(c(sigma.1.hat^2, rho.hat*sigma.1.hat*sigma.2.hat,
              rho.hat*sigma.1.hat*sigma.2.hat, sigma.2.hat^2),
            nrow=2, ncol=2, byrow = T)

L.hat <- matrix(solve(diag(4) - Phi.hat%x%Phi.hat) %*% c(S.hat), ncol=2, nrow=2)

delta.hat <- rep(0, times=m)
for (j1 in 1:m1) {
  for (j2 in 1:m2) {
    J <- pairtok(j1, j2)
    delta.hat[J] <- mvtnorm::pmvnorm(lower=c(a[j1], b[j2]), upper=c(a[j1+1], b[j2+1]), sigma=L.hat)
  }
}
delta.hat <- pmax(delta.hat, 0)
delta.hat <- delta.hat/sum(delta.hat)

Gam.hat <- matrix(0L, ncol=m1*m2, nrow=m1*m2)
for (j in 1:(m1*m2)) {
  J <- pairtok(pairs[j,1], pairs[j,2])
  
  bstar <- 0.5*c(a[pairs[j,1]] + a[pairs[j,1]+1], b[pairs[j,2]] + b[pairs[j,2]+1]) # midpoint of box
  muJ <- as.vector(Phi.hat%*%bstar)
  
  Gam.hat[J,] <- pmax((pbivnorm(x=(a[pairs[,1] + 1] - muJ[1])/sigma.1.hat, y=(b[pairs[,2] + 1] - muJ[2])/sigma.2.hat, rho=rho.hat) - 
                pbivnorm(x=(a[pairs[,1]] - muJ[1])/sigma.1.hat, y=(b[pairs[,2] + 1] - muJ[2])/sigma.2.hat, rho=rho.hat) -
                pbivnorm(x=(a[pairs[,1] + 1] - muJ[1])/sigma.1.hat, y=(b[pairs[,2]] - muJ[2])/sigma.2.hat, rho=rho.hat) + 
                pbivnorm(x=(a[pairs[,1]] - muJ[1])/sigma.1.hat, y=(b[pairs[,2]] - muJ[2])/sigma.2.hat, rho=rho.hat)),0)[pairinv]
}
Gam.hat <- Gam.hat/rowSums(Gam.hat)

P.M3.hat <- function(y) {dpois(y[1], w*beta.1.hat*exp(0.5*(a[pairs[,2]] + a[pairs[,2] + 1]))) * dpois(y[2], w*beta.2.hat*exp(0.5*(b[pairs[,1]] + b[pairs[,1] + 1])))
}

tic()
Xt.hat <- local_dec(init=delta.hat, tpm=Gam.hat, data=t(Y), Fx=P.M3.hat)$preds
toc()
Xt.hat.pair <- ktopair(Xt.hat)

X1.hat <- a.mid[Xt.hat.pair[1,]]
X2.hat <- b.mid[Xt.hat.pair[2,]]

save.image("Model3_1_initialMLE_85.Rdata")


ggplot(data.frame(x=X1.hat), aes(x=x)) +
  geom_histogram(bins=30, color="black", fill="lightblue") +
  geom_vline(aes(xintercept=min(a)), color="red", linetype="dashed") +
  geom_vline(aes(xintercept=max(a)), color="red", linetype="dashed") +
  ylab("Count") +
  xlab(bquote(hat(X)["t,1"])) +
  ggtitle("EV Lac Model 2 - Histogram of Predicted States (Soft Band)")

ggplot(data.frame(x=X2.hat), aes(x=x)) +
  geom_histogram(bins=30, color="black", fill="lightblue") +
  geom_vline(aes(xintercept=min(b)), color="red", linetype="dashed") +
  geom_vline(aes(xintercept=max(b)), color="red", linetype="dashed") +
  ylab("Count") +
  xlab(bquote(hat(X)["t,2"])) +
  ggtitle("EV Lac Model 2 - Histogram of Predicted States (Soft Band)")

preds.bihist <- ggplot(data.frame(x1 = X1.hat, x2=X2.hat), aes(x=x1, y=x2)) +
  geom_bin2d() + 
  scale_fill_gradient(low="lightblue", high="darkblue") + 
  xlab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  ylab(bquote(hat(italic(X))[italic("t,")*plain("2")])) +
  labs(fill = "ObsID 01885") +
  geom_segment(aes(x=min(a), xend=max(a), y=min(b), yend=min(b)), color="red", linetype="dashed") +
  geom_segment(aes(x=min(a), xend=max(a), y=max(b), yend=max(b)), color="red", linetype="dashed") +
  geom_segment(aes(x=min(a), xend=min(a), y=min(b), yend=max(b)), color="red", linetype="dashed") +
  geom_segment(aes(x=max(a), xend=max(a), y=min(b), yend=max(b)), color="red", linetype="dashed") +
  theme(legend.position = c(0.86, 0.4), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))

preds.bihist
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_M3_predsbihist.png"), preds.bihist, height=3, width=7)


  
