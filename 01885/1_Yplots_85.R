rm(list = ls())

library(tidyverse)
library(FITSio)

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

gg.Yboth <- ggplot(data=data.frame(LOW=Y[,"LOW"], HIGH=Y[,"HIGH"], t=1:TT), mapping=aes(x=t)) +
  geom_line(aes(y=HIGH, color="Hard band"), linewidth=0.5)  +
  geom_line(aes(y=LOW, color="Soft band"), linewidth=0.5)  +
  ylab(bquote( bold(Y[italic(t)]))) +
  xlab(bquote( italic(t))) +
  scale_color_manual(values=c(rgb(63, 141, 174, maxColorValue = 255), rgb(198, 101, 38, maxColorValue = 255))) +
  labs(color="ObsID 01885") +
  theme(legend.position = c(0.9, 0.8), legend.key = element_blank(), legend.text = element_text(size=12),
        legend.box.background = element_rect(colour = "black"), 
        axis.title = element_text(size=14.5, family="Times New Roman"),
        text=element_text(family="Times New Roman"),
        axis.text = element_text(size=14))

gg.Yboth
ggsave(file=paste0(getwd(),"/Figures/EVLac_85_Yboth.png"), gg.Yboth, height=3, width=7)
