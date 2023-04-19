rm(list=ls())
library(readr)
library(ggplot2)
library(dplyr)


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

########overlapping

#overlapping_avar_fn=function(y,m){
y <- white.noise
m = 2
  M=length(y)
  
  numberOfGroups=M-2*m
  groupmeans = numeric(numberOfGroups)
  
  for(i in 1:(M-2*m)){
    groupmeans[i]= mean(y[i:(i+m)])  
  }
  
  out=0
  for(i in 1:(M-3*m)){
    out=out+(groupmeans[i+m]-groupmeans[i])^2
  }
  out=out/(2*(M-2*m+1))
  
  return(out)
#}


getAvars=function(N){
  y=rnorm(N,0,1) 
  # taus = c(1,2,3,5,10,20)
  # taus = c(taus,seq(40,N/2,by=100)) ## this goes too far, gives lots of NAs from avar_fn, fix later 
  
  maxn = ceiling((N-1)/2) 
  p = floor (log10(maxn)/log10(2)) #Number of clusters
  taus=2^(0:p)
  
  avars=numeric(length(taus))
  overlapping_avars=numeric(length(taus))
  
  for (i in 1:length(taus)){
    avars[i]=avar_fn(y,taus[i])
    overlapping_avars[i]=overlapping_avar_fn(y,taus[i])
  }
  
  m1=data.frame(taus=taus,avars=avars)  
  fit=lm(log(sqrt(avars))~log(taus),data = m1)
  slope=as.numeric(fit$coefficients[2])
  int=as.numeric(fit$coefficients[1])
  
  avarRes=data.frame(taus=taus,avars=avars, overavars=overlapping_avars,N=N,slope=slope,int=int)  

  ##########################################
  # get SE
  # ##########################################

  m2=data.frame(taus=taus,oavars=overlapping_avars)  
  fit2=lm(log(sqrt(oavars))~log(taus),data = m2)
  slope2=as.numeric(fit2$coefficients[2])
  int2=as.numeric(fit2$coefficients[1])
  
  SEests=data.frame()
  
  onew=data.frame(N=N,out=exp(int2+slope2*log(N)), type="OAD")
  new=data.frame(N=N,out=exp(int+slope*log(N)), type="AD")
  new2=data.frame(N=N,out=sd(y)/sqrt(N), type="SE")
  new3=data.frame(N=N,out=1/sqrt(N), type="true")
  SEests=bind_rows(new,new2,new3,onew)

  return(list(avarRes=avarRes,SEests=SEests))
}



### increase N, what does adev converge to, for which tau
Ns=c(100,1000,4000,10000,20000)

allavarRes=data.frame()
allSEests=data.frame()
for(j in 1:10){
  for(N in Ns){
    tmp=getAvars(N)
    
    allavarRes=bind_rows(allavarRes,tmp$avarRes)
    allSEests=bind_rows(allSEests,tmp$SEests)
  }
  
}

ggplot(allavarRes,aes(taus,sqrt(avars)))+
  geom_point()+
  geom_point(aes(taus,sqrt(overavars)),col="blue")+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  # stat_smooth(method = "lm",se = F)+
  facet_wrap(~N)#,scales = "free")

ggplot(allSEests,aes(N,out,color=type))+
  geom_point()+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")
  # geom_jitter()

ggplot(allSEests,aes(as.factor(N),out,color=type))+
  scale_y_continuous(trans = "log")+
  geom_boxplot()

### ADEV results look biased, not converging to what I expect. 
### Maybe overlapping adev helps this? 
### Or Modified Allan Deviation?


ggplot(filter(allavarRes,N==10000),aes(taus,sqrt(overavars),col=as.factor(slope)))+
  geom_point()+
  # geom_point(aes(taus,sqrt(overavars)),col="blue")+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  stat_smooth(method = "lm",se = F)+
  facet_wrap(~N)#,scales = "free")


