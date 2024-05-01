rm(list = ls())

library(tidyverse)
library(FITSio)
library(depmixS4)

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


####### Fit basic 2- and 3-state HMMs to the observed data ######

mod.2state <- depmix(list(LOW ~ 1, HIGH ~ 1), data=as.data.frame(Y), nstates=2, family = list(poisson(), poisson()))
fit.2state <- fit(mod.2state)
Z.2state <- factor(fit.2state@posterior$state - 1, labels=c("Flaring", "Quiescent"))

gg.2state <- ggplot(data=data.frame(LOW=Y[,"LOW"], t=1:TT), mapping=aes(x=t)) +
  geom_point(aes(y=LOW, color=Z.2state), size=0.5)  +
  scale_color_manual(values=c("blue", "red")) +
  ylab(bquote(italic(Y)[italic("t,")*plain("1")])) +
  xlab(expression(paste("Time bin index ", italic(t), " (", Delta, t, " = 50 s",")"))) +
  labs(color="ObsID 01885") +
  theme(legend.position = c(0.9, 0.8), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))

gg.2state
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_2state.png"), gg.2state, height=3, width=7)



mod.3state <- depmix(list(LOW ~ 1, HIGH ~ 1), data=as.data.frame(Y), nstates=3, family = list(poisson(), poisson()))
fit.3state <- fit(mod.3state)
Z.3state <- factor(c(3,1,2)[fit.3state@posterior$state], labels=c("Flaring 1", "Quiescent", "Flaring 2"))

gg.3state <- ggplot(data=data.frame(LOW=Y[,"LOW"], t=1:TT), mapping=aes(x=t)) +
  geom_point(aes(y=LOW, color=Z.3state), size=0.5)  +
  scale_color_manual(values=c("darkgreen", "red", "blue")) +
  ylab(bquote(italic(Y)[italic("t,")*plain("1")])) +
  xlab(expression(paste("Time bin index ", italic(t), " (", Delta, t, " = 50 s",")"))) +
  labs(color="ObsID 01885") +
  theme(legend.position = c(0.9, 0.75), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))
gg.3state
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_3state.png"), gg.3state, height=3, width=7)
