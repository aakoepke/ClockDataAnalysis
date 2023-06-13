####################### Plots of Simulation Results ##########################
library(magrittr)

setwd("/home/cmb15/ClockDataAnalysis/Code/Paper1/")
source("../SA_ImportantFunctions.R")
################# White Noise No Gaps ###################
tmat <- readRDS(file = "Results/tmat050523_W12_K5_N2048_300sims_WhiteNoiseNoGaps.Rds")
bmat <- readRDS("Results/bmat050523_W12_K4_N2048_300sims_WhiteNoiseNoGaps.Rds")
oamat <- readRDS(file = "Results/oamat052623_N2048_300sims_WhiteNoiseNoGaps.Rds")

### Plots ###
N <- 1024
taus <- c(2^(0:9),floor(N/3))
numberOfSimulations = 300

##tidy the data
tmat %<>% t()
bmat %<>% t()
dim(bmat)
three.mat <- rbind(oamat, tmat)
method_labels <- rep(c("avar", "tr"), each = numberOfSimulations)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

#sumDat=dat %>% group_by(method,tau) %>%
# summarise(median = median(measurement),
#          lower25=quantile(measurement,prob=.25),
#            upper75=quantile(measurement,prob=.75),
#           min=min(measurement),
#          max=max(measurement))

# breaks <- 10^(-10:10)
# minor_breaks <- rep(1:9,21 )*rep(10^(-10:10), each = 9)
# 
# linedat_w=data.frame(tau=1:1000)
# linedat_w$truth=1/linedat_w$tau#+linedat$tau
# linedat_w$method=NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat_w,aes(tau,truth), linewidth = 1.2)+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("Allan","Spectral"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here
  ### basically specifies the x and y position in terms of the plot size in unit scale.
  theme_bw(base_size = 20)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  annotation_logticks()+

  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  ggtitle("White Noise Comparison, No Gaps")


################# White Noise with Gaps #################



################ Flicker Noise No Gaps ##################
tmat_flk <- readRDS(file = "Results/tmat052623_W3_K3_N2048_300sims_FlickerNoiseNoGaps.Rds")
bmat_flk <- readRDS(file = "Results/bmat052623_W3_K3_N2048_300sims_FlickerNoiseNoGaps.Rds")
oamat_flk <- readRDS(file = "Results/oamat053123_N2048_300sims_FlickerNoiseNoGaps.Rds")


### Plots ###
N <- 1024
taus <- c(2^(0:8),floor(N/3))
numberOfSimulations = 300
tmat_flk <- tmat
bmat_flk <- bmat
oamat_flk <- oamat
##tidy the data
tmat_flk %<>% t()
bmat_flk %<>% t()
dim(bmat_flk)
three.mat <- rbind(oamat_flk, bmat_flk)
method_labels <- rep(c("avar", "tr"), each = numberOfSimulations)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9,21 )*rep(10^(-10:10), each = 9)

# linedat = data.frame(tau=2:floor(N/3))
# linedat$truth = tavar_ARFIMA(floor(N/3), d = 0.25, sig.2.a = 1)
# linedat$method = NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat,aes(tau,truth), linewidth = 1.2)+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("Allan","Spectral"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here
  ### basically specifies the x and y position in terms of the plot size in unit scale.
  theme_bw(base_size = 20)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  annotation_logticks()+
  
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  ggtitle("Flicker Noise Comparison, No Gaps")


















##################################### look @ later ######################################
##tidy the data
tmat_flk %<>% t()
bmat_flk %<>% t()
dim(bmat_flk)
three.mat <- rbind(oamat_flk, tmat_flk,bmat_flk)
method_labels <- rep(c("avar", "tr", "bp"), each = numberOfSimulations)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

sumDat=dat %>% group_by(method,tau) %>%
  summarise(median = median(log10(measurement)),
            lower25=quantile(log10(measurement),prob=.25),
            upper75=quantile(log10(measurement),prob=.75),
            min=min(log10(measurement)),
            max=max(log10(measurement)))

ggplot(data = sumDat,aes(x=log10(tau),y=median,col=method,fill=method))+
  ### width controls how wide these are, can change all below if they are too wide to see tau linear relationship well with real data
  geom_errorbar(aes(ymin=min,ymax=max),width=.1,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),width=.1,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),color="black",width=.1,position = "dodge")+
  ####
  ylab(expression(log[10](sigma(tau))))+
  xlab(expression(log[10](tau)))+
  ggtitle("Flicker Noise, No Gaps")+
  ### add true line below
  geom_abline(slope = 0,intercept = 3.1,size=1)




################ Flicker Noise with Gaps ################