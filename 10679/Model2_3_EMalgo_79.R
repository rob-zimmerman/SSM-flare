rm(list = ls())

library(tidyverse)
library(FITSio)
library(evmix)
library(xtable)
library(tictoc)
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



####### Make state predictions #######

pars.MLE <- c(0.98904829, 0.09278349, 0.13569292, 0.14521816, 0.06126551)
phi.1.MLE <- pars.MLE[1]
sigma.1.MLE <- pars.MLE[2]
sigma.2.MLE <- pars.MLE[3]
beta.1.MLE <- pars.MLE[4]
beta.2.MLE <- pars.MLE[5]

m <- 40
cmin <- -1.5
cmax <- 2
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

# intializations

K <- 3
mu <- 0*1:K
var <- 0*1:K

km <- kmeans(X1.MLE, centers=K)

p <- sapply(X=1:K, FUN = function(k) mean(km$cluster == k))
mu <- sapply(X=1:K, FUN = function(k) mean(X1.MLE[which(km$cluster == k)]))
var <- sapply(X=1:K, FUN = function(k) var(X1.MLE[which(km$cluster == k)]))

g <- matrix(0L, nrow=K, ncol=TT)

err <- 10e6

while (err > 0.0001) {
  par.old <- c(mu, var, p)
  
  for (k in 1:K) {
    g[k,] <- p[k]*dnorm(x=X1.MLE, mean=mu[k], sd=sqrt(var[k]))
  }
  
  g <- t(apply(g, MARGIN=1, FUN=function(x) x/colSums(g)))

  for (k in 1:K) {
    mu[k] <- sum(g[k,]*X1.MLE)/sum(g[k,])
    var[k] <- sum(g[k,]*(X1.MLE - mu[k])^2)/sum(g[k,])
    p[k] <- mean(g[k,])
  }

  err <- sqrt(sum(par.old - c(mu, var, p))^2)
}

k_low <- which(mu < 0)
f.hat <- lapply(1:K, function(k) {function(x) {dnorm(x, mean=mu[k], sd=sqrt(var[k]))}})
ftot <- Vectorize(function(x) {sum(sapply(1:K, function(k) p[k]*f.hat[[k]](x)))})
fq <- Vectorize(function(x) {sum(sapply(k_low, function(k) p[k]*f.hat[[k]](x)))})

# make posterior flare classifications
Z.EM <- fq(X1.MLE)/ftot(X1.MLE)
Z.EM <- 1 - Z.EM
Z.01 <- 1*(Z.EM <= 0.5)


gg.statepreds <- ggplot(data=data.frame(LOW=X1.MLE, t=1:TT), mapping=aes(x=t)) +
  geom_point(aes(y=LOW), size=0.5)  +
  scale_color_manual(values=c("black")) +
  ylab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  xlab(bquote(italic(t))) +
  annotate("label", label="ObsID 10679", x=1800, y=1.5, size=11/.pt, family="Times New Roman") +
  theme(legend.position = c(0.9, 0.75), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))
gg.statepreds
ggsave(file=paste0(getwd(),"/Figures/EVLac_79_M2_statepreds.png"), gg.statepreds, height=3, width=7)



gg.postprobs <- ggplot(data=data.frame(LOW=X1.MLE, t=1:TT), mapping=aes(x=t)) +
  geom_point(aes(y=LOW, color=Z.EM), size=0.5)  +
  scale_color_gradient(low="red", high="blue") +
  ylab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  xlab(bquote(italic(t))) +
  labs(color="ObsID 10679") +
  theme(legend.position = c(0.9, 0.7), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))

gg.postprobs
ggsave(file=paste0(getwd(),"/Figures/EVLac_79_M2_postprobs.png"), gg.postprobs, height=3, width=7)



gg.postprobsY <- ggplot(data=data.frame(LOW=Y[,"LOW"], t=1:TT), mapping=aes(x=t)) +
  geom_point(aes(y=LOW, color=Z.EM), size=0.5)  +
  scale_color_gradient(low="red", high="blue") +
  ylab(bquote(italic(Y)[italic("t,")*plain("1")])) +
  xlab(bquote(italic(t))) +
  labs(color="ObsID 10679") +
  theme(legend.position = c(0.9, 0.7), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))

gg.postprobsY
ggsave(file=paste0(getwd(),"/Figures/EVLac_79_M2_postprobsY.png"), gg.postprobsY, height=3, width=7)



gg.componentdists <- ggplot(data=data.frame(LOW=X1.MLE), mapping=aes(x=LOW)) +
  geom_histogram(aes(y=after_stat(density)), bins=30, color="darkblue", fill="lightblue") +
  stat_function(aes(color="red"), fun= function(x) p[1]*f.hat[[1]](x), linewidth=1) +
  stat_function(aes(color="blue"), fun= function(x) p[2]*f.hat[[2]](x), linewidth=1) + 
  stat_function(aes(color="lightblue"), fun= function(x) p[3]*f.hat[[3]](x), linewidth=1) + 
  xlab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  ylab("Fitted density") +
  labs(color="ObsID 10679") +
  theme(legend.position = c(0.72, 0.7), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14)) +
  scale_color_manual(labels = c(expression(hat(alpha)[3] %.% hat(italic(f))[3]), expression(hat(alpha)[1] %.% hat(italic(f))[1]), expression(hat(alpha)[2] %.% hat(italic(f))[2])), values = c("blue", "red", "navyblue"))

gg.componentdists
ggsave(file=paste0(getwd(),"/Figures/EVLac_79_M2_componentdists.png"), gg.componentdists, height=3, width=7)



gg.mixturedist <- ggplot(data=data.frame(LOW=X1.MLE), mapping=aes(x=LOW)) +
  geom_histogram(aes(y=after_stat(density)), bins=30, color="darkblue", fill="lightblue") +
  stat_function(aes(color="purple"), fun= ftot, linewidth=1) +
  xlab(bquote(hat(italic(X))[italic("t,")*plain("1")])) +
  ylab("Fitted density") +
  labs(color="ObsID 10679") +
  theme(legend.position = c(0.8, 0.8), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14)) +
  scale_color_manual(labels = c(expression(hat(alpha)[1] %.% hat(italic(f))[1] + hat(alpha)[2] %.% hat(italic(f))[2] + hat(alpha)[3] %.% hat(italic(f))[3])), values = c("purple"))

gg.mixturedist
ggsave(file=paste0(getwd(),"/Figures/EVLac_79_M2_mixturedists.png"), gg.mixturedist, height=3, width=7)



####### Make a table with the estimated parameters #######

EM.ests <- data.frame(Component = factor(1:K), alpha = p, mu = mu, tau2 = var)
colnames(EM.ests) <- c("Component $k$", "$\\hat{\\alpha}_k$", "$\\hat{\\mu}_k$", "$\\hat{\\tau}_k^2$")
print(xtable(EM.ests, type = "latex", digits = 4, display = rep("f", times = 5)), booktabs = TRUE, sanitize.text.function = function(x) x, math.style.exponents = F, include.rownames=FALSE)


####### Create intervals of flaring activity #######

changepts.M2 <- c(0, which(abs(diff(Z.01)) > 0)) + c(1,0)
changepts.M2 <- matrix(changepts.M2, ncol=2, byrow=T)
intervals.M2 <- data.frame(t_start = changepts.M2[,1], t_end=changepts.M2[,2],
                           t_space_start = changepts.M2[,1]*w + Y.t_min,
                           t_space_end = (changepts.M2[,2]+1)*w + Y.t_min)
write.csv(intervals.M2, file=paste0(getwd(),"/Intervals/EVLac79_flareintervals_M2_w50.csv"), quote=FALSE)
