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

getAvars=function(N){
  y=rnorm(N,0,1) 
  taus = c(1,2,3,5,10,20)
  taus = c(taus,seq(40,N,by=100)) ## this goes too far, gives lots of NAs from avar_fn, fix later 
  
  avars=numeric(length(taus))
  
  for (i in 1:length(taus)){
    avars[i]=avar_fn(y,taus[i])
  }
  m1=data.frame(taus=taus,avars=avars)  
  fit=lm(log(sqrt(avars))~log(taus),data = m1)
  slope=as.numeric(fit$coefficients[2])
  int=as.numeric(fit$coefficients[1])
  
  avarRes=data.frame(taus=taus,avars=avars,N=N,slope=slope,int=int)  

  ##########################################
  # get SE
  ##########################################
  SEests=data.frame()
  
  new=data.frame(N=N,out=exp(int+slope*log(N)), type="AD")
  new2=data.frame(N=N,out=sd(y)/sqrt(N), type="SE")
  new3=data.frame(N=N,out=1/sqrt(N), type="true")
  SEests=bind_rows(new,new2,new3)
  
  return(list(avarRes=avarRes,SEests=SEests))
}

### increase N, what does adev converge to, for which tau
Ns=c(100,1000,4000,10000,20000)

allavarRes=data.frame()
allSEests=data.frame()
for(j in 1:1){
  for(N in Ns){
    tmp=getAvars(N)
    
    allavarRes=bind_rows(allavarRes,tmp$avarRes)
    allSEests=bind_rows(allSEests,tmp$SEests)
  }
  
}

ggplot(allavarRes,aes(taus,sqrt(avars)))+
  geom_point()+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  stat_smooth(method = "lm",se = F)+
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

########overlapping

y=rnorm(10)
m=5
overlapping_avar_fn=function(y,m){
  # tau=m?
  M=length(y)
  
  numberOfGroups=M-m+1
  groupmeans = numeric(numberOfGroups)
  
  for(i in 1:(M+1-m)){
    groupmeans[i]=mean(y[i:(i+m-1)])  
  }


  for(j in 1:(M-2*m+1)){
    for(i in j:(j+m-1))
  }
  
  1/(2*m^2*(M-2*m+1))
  

    
    
    
    
        
  
  
  div=seq(1,n,by = tau)
  
  M=length(div)-1 #number of groups
  
  groupmeans = numeric(M)
  for(i in 1:M){
    groupmeans[i]=mean(y[div[i]:(div[i+1]-1)])
  }
  
  1/(2*(M-1)) * sum(diff(groupmeans)^2)
}
