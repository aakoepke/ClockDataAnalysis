

library(ggplot2)
library(dplyr)

load("Data/TidyDataforBoxplots.RData")

# a=rnorm(30,10,1)
# b=rnorm(30,10,1)
# c=rnorm(30,10,1)
# dat=data.frame(method=rep(c("a","b","c"),each=30),
#                measurement=c(a,b,c),
#                tau=rep(1:5))
# 
# 
sumDat=dat %>% group_by(method,tau) %>%
  summarise(median = median(measurement),
            lower25=quantile(measurement,prob=.25),
            upper75=quantile(measurement,prob=.75),
            min=min(measurement),
            max=max(measurement))


##### don't use this anymore, found a way to use regular boxplots below
# histwidth=.1
# ggplot(data = sumDat,aes(x=tau,y=median,col=method,fill=method))+
#   ### width controls how wide these are, can change all below if they are too wide to see tau linear relationship well with real data
#   geom_errorbar(aes(ymin=min,ymax=max),width=histwidth,position = "dodge")+
#   geom_crossbar(aes(ymin=lower25,ymax=upper75),width=histwidth,position = "dodge")+
#   geom_crossbar(aes(ymin=lower25,ymax=upper75),color="black",width=histwidth,position = "dodge")+
#   ####
#   scale_y_log10()+
#   scale_x_log10()+
#   ylab("adev")+
#   xlab(expression(tau))+
#   ### add true line below
#   geom_abline(slope = .01,intercept = 1,size=1)+
#   guides(fill=guide_legend(title="title"),color=guide_legend(title="title"))
#   # scale_color_manual(values = c("Method A name","Method B name","Method C name"))#,col="legend title")
# 


linedat=data.frame(tau=1:1000)
linedat$truth=1/linedat$tau+linedat$tau
linedat$method=NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot()+
  scale_y_log10()+
  scale_x_log10()+
  ylab("adev")+
  xlab(expression(tau))+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat,aes(tau,truth))+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("Allan","Bandpass","Transfer Function"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here 
  ### basically specifies the x and y position in terms of the plot size in unit scale. 
  theme(legend.position = c(.15, .2))


# ################## boxplots with transformation first
# dat$log10tau=(log10(dat$tau))
# dat$log10meas=log10(dat$measurement)
# 
# ggplot(dat,aes(log10tau,log10meas,col=method,group=interaction(log10tau,method)))+
#   geom_boxplot()+
#   scale_color_discrete(labels= c("Allan","Bandpass","Transfer Function"),name="Method")
# ####################################

  

