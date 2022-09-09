

library(ggplot2)
library(dplyr)

a=rnorm(30,10,1)
b=rnorm(30,10,1)
c=rnorm(30,10,1)
dat=data.frame(method=rep(c("a","b","c"),each=30),
               measurement=c(a,b,c),
               tau=rep(1:5))


sumDat=dat %>% group_by(method,tau) %>%
  summarise(median = median(measurement),
            lower25=quantile(measurement,prob=.25),
            upper75=quantile(measurement,prob=.75),
            min=min(measurement),
            max=max(measurement))

           
ggplot(data = sumDat,aes(x=tau,y=median,col=method,fill=method))+
  ### width controls how wide these are, can change all below if they are too wide to see tau linear relationship well with real data
  geom_errorbar(aes(ymin=min,ymax=max),width=.05,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),width=.05,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),color="black",width=.05,position = "dodge")+
  ####
  scale_y_log10()+
  scale_x_log10()+
  ylab(expression(log[10]("some nonsense")))+
  xlab(expression(log[10](tau)))+
  ### add true line below
  geom_abline(slope = .01,intercept = 1,size=1)
