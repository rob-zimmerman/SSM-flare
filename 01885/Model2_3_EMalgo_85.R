rm(list = ls())

library(tidyverse)
library(FITSio)
library(evmix)
library(xtable)
library(tictoc)
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



####### Make state predictions #######

pars.MLE <- c(0.97964369, 0.10071215, 0.16168891, 0.19381708, 0.06241664)
phi.1.MLE <- pars.MLE[1]
sigma.1.MLE <- pars.MLE[2]
sigma.2.MLE <- pars.MLE[3]
beta.1.MLE <- pars.MLE[4]
beta.2.MLE <- pars.MLE[5]

m <- 40
cmin <- -1.25
cmax <- 2.65
c <- seq(from=cmin, to=cmax, length.out=(m+1))
c.mid <- c[-length(c)] + diff(c)/2

P.M2 <- function(y, phi.1, sigma.1, sigma.2, beta.1, beta.2) {
  p <- rep(0, times=m)
  for (i in 1:m) {
    cstar.1 <- c.mid[i]
    cstar.2 <- (sigma.2/sigma.1)*cstar.1
    p[i] <- dpois(y[1], w*beta.1*exp(cstar.1))*dpois(y[2], w*beta.2*exp(cstar.2))
  }
  return(diag(p))
}

delta.MLE <- rep(0, m)
for (j in 1:m) {
  delta.MLE[j] <- pnorm(c[j+1], mean=0, sd=sqrt(sigma.1.MLE^2/(1-phi.1.MLE^2))) - pnorm(c[j], mean=0, sd=sqrt(sigma.1.MLE^2/(1-phi.1.MLE^2)))
}
delta.MLE <- delta.MLE/sum(delta.MLE)

Gam.MLE <- matrix(0L, nrow=m, ncol=m)
for (i in 1:m) {
  for (j in 1:m) {
    Gam.MLE[i,j] <- pnorm(c[j+1], mean=phi.1.MLE*c.mid[i], sd=sigma.1.MLE) - pnorm(c[j], mean=phi.1.MLE*c.mid[i], sd=sigma.1.MLE)
  }
}
Gam.MLE <- Gam.MLE/rowSums(Gam.MLE)

P.M2.MLE <- function(y) {asFinite(diag(P.M2(y, phi.1.MLE, sigma.1.MLE, sigma.2.MLE, beta.1.MLE, beta.2.MLE)))}

X1.MLE <- c.mid[local_dec(init=delta.MLE, tpm=Gam.MLE, data=t(Y), Fx=P.M2.MLE)$preds] # predicted X1t's



####### EM algorithm for binary state classifications #######

Xq <- X1.MLE[1:750]
Xr <- X1.MLE[751:TT]

X1.jit <- X1.MLE + runif(n=TT, -0.0975/2, 0.0975/2)

Xq <- X1.jit[1:750]
Xr <- X1.jit[751:TT]

# intializations
K <- 25

Pi <- rep(1/K, times=K)

bseq <- seq(qkden(p=0.5, kerncentres=Xq), cmax, length.out=K)
int.length <- bseq[K] - bseq[K-1]

f1.hat <- function(x) {dkden(x, Xq)} # KDE for first 750 states

f2.hat <- function(x, pivec) { 
  ifelse(x <= bseq[1] || x > bseq[K], 0, pivec[findInterval(x=x, vec=bseq, left.open=FALSE, rightmost.closed=FALSE)]/int.length)
}

f2.hat <- Vectorize(f2.hat, "x")

alphap <- 0.5

err <- 10e6

phi <- 750/TT # proportion that the subset constitutes

while (err > 0.0001) {
  alphap.old <- alphap
  Pi.old <- Pi
  
  a1 <- f1.hat(Xr)*alphap/ (f1.hat(Xr)*alphap + f2.hat(Xr, Pi)*(1-alphap))
  a2 <- 1 - a1
  
  alphap <- mean(a1)
  
  suma2 <- sum(a2)
  
  for (k in 1:K) {
    Pi[k] <- sum(1*between(Xr, bseq[k], bseq[k+1])*a2)/suma2
  }
  
  print(c(alphap, Pi))
  err <- sqrt(sum( (alphap - alphap.old)^2 + (Pi - Pi.old)^2  ))
}

alpha <- alphap + phi*(1-alphap)

# make posterior flare classifications
Z.EM <- 1 - c(rep(1,750), alpha*f1.hat(Xr)/(alpha*f1.hat(Xr) + (1-alpha)*f2.hat(Xr, Pi)))

Z.01 <- 1*(Z.EM > 0.5)

gg.statepreds <- ggplot(data=data.frame(LOW=X1.MLE, t=1:TT), mapping=aes(x=t)) +
  geom_point(aes(y=LOW), size=0.5)  +
  scale_color_manual(values=c("black")) +
  ylab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  xlab(expression(paste("Time bin index ", italic(t), " (", Delta, t, " = 50 s",")"))) +
  annotate("label", label="ObsID 01885", x=1850, y=1.9, size=11/.pt, family="Times New Roman") +
  theme(legend.position = c(0.9, 0.75), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))

gg.statepreds
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_M2_statepreds.png"), gg.statepreds, height=3, width=7)



gg.postprobs <- ggplot(data=data.frame(LOW=X1.MLE, t=1:TT), mapping=aes(x=t)) +
  geom_point(aes(y=LOW, color=Z.EM), size=0.5)  +
  scale_color_gradient(low="red", high="blue") +
  ylab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  xlab(expression(paste("Time bin index ", italic(t), " (", Delta, t, " = 50 s",")"))) +
  labs(color="ObsID 01885") +
  theme(legend.position = c(0.1, 0.65), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))

gg.postprobs
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_M2_postprobs.png"), gg.postprobs, height=3, width=7)



gg.postprobsY <- ggplot(data=data.frame(LOW=Y[,"LOW"], t=1:TT), mapping=aes(x=t)) +
  geom_point(aes(y=LOW, color=Z.EM), size=0.5)  +
  scale_color_gradient(low="red", high="blue") +
  ylab(bquote(italic(Y)[italic("t,")*plain("1")])) +
  xlab(expression(paste("Time bin index ", italic(t), " (", Delta, t, " = 50 s",")"))) +
  labs(color="ObsID 01885") +
  theme(legend.position = c(0.1, 0.65), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))

gg.postprobsY
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_M2_postprobsY.png"), gg.postprobsY, height=3, width=7)



gg.componentdists <- ggplot(data=data.frame(LOW=X1.MLE), mapping=aes(x=LOW)) +
  geom_histogram(aes(y=after_stat(density)), bins=29, color="darkblue", fill="lightblue") +
  stat_function(aes(color="blue"), fun= function(x) (1-alpha)*f2.hat(x, Pi), linewidth=1) + 
  stat_function(aes(color="red"), fun= function(x) alpha*f1.hat(x), linewidth=1) +
  xlab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  ylab("Fitted density") +
  labs(color="ObsID 01885") +
  theme(legend.position = c(0.76, 0.73), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14)) +
  scale_color_manual(labels = c(expression((1-alpha) %.% hat(italic(f))[2]), expression(alpha %.% hat(italic(f))[1])), values = c("blue", "red"))

gg.componentdists
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_M2_componentdists.png"), gg.componentdists, height=3, width=7)



gg.mixturedist <- ggplot(data=data.frame(LOW=X1.MLE), mapping=aes(x=LOW)) +
  geom_histogram(aes(y=after_stat(density)), bins=29, color="darkblue", fill="lightblue") +
  stat_function(aes(color="purple"), fun= function(x) alpha*f1.hat(x) + (1-alpha)*f2.hat(x,Pi), linewidth=1) +
  xlab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  ylab("Fitted density") +
  labs(color="ObsID 01885") +
  theme(legend.position = c(0.8, 0.8), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14)) +
  scale_color_manual(labels = c(expression(alpha %.% hat(italic(f))[1] + (1-alpha) %.% hat(italic(f))[2])), values = c("purple"))

gg.mixturedist
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_M2_mixturedists.png"), gg.mixturedist, height=3, width=7)



####### Create intervals of flaring activity #######

changepts.M2 <- which(abs(diff(Z.01)) > 0) + c(1,0)
changepts.M2 <- matrix(changepts.M2, ncol=2, byrow=T)
intervals.M2 <- data.frame(t_start = changepts.M2[,1], t_end=changepts.M2[,2],
                           t_space_start = changepts.M2[,1]*w + Y.t_min,
                           t_space_end = (changepts.M2[,2]+1)*w + Y.t_min)
write.csv(intervals.M2, file=paste0(getwd(),"/Intervals/EVLac85_flareintervals_M2_w50.csv"), quote=FALSE)
