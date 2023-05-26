####################### Plots of Simulation Results ##########################

################# White Noise No Gaps ###################


#for plotting: bmat, tmat, amat/oamat


### Plots ###
N <- 2048
taus <- c(2^(0:9),floor(N/3))

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

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9,21 )*rep(10^(-10:10), each = 9)

linedat_w=data.frame(tau=1:1000)
linedat_w$truth=1/linedat_w$tau#+linedat$tau
linedat_w$method=NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat_w,aes(tau,truth), size = 1.2)+
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




################ Flicker Noise with Gaps ################