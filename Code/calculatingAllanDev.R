rm(list=ls())
library(readr)
library(ggplot2)
library(dplyr)
# primerData1 <- read_csv("Data/primerData1.txt", col_names = FALSE, skip = 10)
# colnames(primerData1)="X1tilde"
# primerData1$t = 1:dim(primerData1)[1]
# 
# y=primerData1$X1tilde
# tau=10

# y=rnorm(4000,0,1) ### increase this, what does it converge to, for which tau

avar_fn=function(y,tau){
  n=length(y)
  
  div=seq(1,n,by = tau)

  M=length(div)-1 #number of groups
  
  groupmeans = numeric(M)
  for(i in 1:M){
    groupmeans[i]=mean(y[div[i]:(div[i+1]-1)])
  }

  1/(2*(M-1)) * sum(diff(groupmeans)^2)
}


getAvars=function(N){
  y=rnorm(N,0,1) 
  taus = c(1,2,3,5,10,20)#floor(N/8),floor(N/4))#200,1000)
  taus = c(taus,seq(40,N,by=100))
  # print(taus)
  avars=numeric(length(taus))
  
  for (i in 1:length(taus)){
    avars[i]=avar_fn(y,taus[i])
  }
  res=data.frame(taus=taus,avars=avars)  
  fit=lm(log(avars)~log(taus),data = res)
  slope=as.numeric(fit$coefficients[2])
  int=as.numeric(fit$coefficients[1])
  
  res=data.frame(taus=taus,avars=avars,N=N,slope=slope,int=int)  

  ##########################################
  # get SE
  ##########################################
  keep=data.frame()
  
  new=data.frame(N=N,out=exp(int+slope*log(N)), type="AD")
  new2=data.frame(N=N,out=sd(y)/N, type="SE")
  new3=data.frame(N=N,out=1/N, type="true")
  keep=bind_rows(keep,new,new2,new3)
  
  return(list(res=res,keep=keep))
}


Ns=c(100,1000,4000,10000,20000)

allres=data.frame()
allkeep=data.frame()
for(j in 1:10){
  for(N in Ns){
    tmp=getAvars(N)
    
    allres=bind_rows(allres,tmp$res)
    allkeep=bind_rows(allkeep,tmp$keep)
  }
  
}

ggplot(allres,aes(taus,avars))+
  geom_point()+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  stat_smooth(method = "lm",se = F)+
  facet_wrap(~N,scales = "free")

ggplot(allkeep,aes(N,out,color=type))+
  geom_point()+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  geom_jitter()

# allkeep$N = as.factor(allkeep$N)
ggplot(allkeep,aes(as.factor(N),out,color=type))+
  # scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  geom_boxplot()

# 
#   # geom_abline(slope = res$slope[1],intercept = res$int[1])+
#   # coord_cartesian(xlim = c(1,1.5),ylim = c(1,1.4))
#   
#   
#   # print(exp(predict(fit))-res$avars)
#   
#   # print(exp(int+slope*log(1)))
#   
#   #extrapolate to size of data set
#   new=data.frame(N=N,out=exp(int+slope*log(N)), type="AD")
#   new2=data.frame(N=N,out=sd(y)/N, type="SE")
#   new3=data.frame(N=N,out=1/N, type="true")
#   keep=bind_rows(keep,new,new2,new3)
#   
#   #SE of the mean
#   # sd(y)/N
#   
# 
# keep
# 
# ggplot(keep,aes(N,out,color=type))+
#   geom_point()+
#   scale_x_continuous(trans = "log")+
#   scale_y_continuous(trans = "log")
