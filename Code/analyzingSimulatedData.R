# source("/home/aak3/NIST/ClockDataAnalysis/Code/analyzingSimulatedData.R")
######################################################################
### relevant libraries and cod
######################################################################
rm(list=ls())

library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library("openxlsx")

##### read in data, don't remember where it came from
library("R.matlab")

test=readMat("Data/FlickerNoise/pink_noise_2000000.mat")
plot(test$r[1,seq(1,2000000,by=1000)])


#####################make a function
source("/home/aak3/NIST/ClockDataAnalysis/Code/SA_ImportantFunctions.R")

calcAvar=function(thetau,freq,spec.hat,delta.f){
  # G_tau vector length number of frequencies 
  G_tau <- transfer.func(freq,tau = thetau) #change the tau value to get different vectors
  G_tau[1] <- 0 # this was 1 in the old code, but should be 0
  
  avar=G_tau%*%spec.hat*delta.f
  
  return(avar)
}

calculateAvars=function(x.t,t.vec,taus,N.fourier=100,numTapers=3,calcCov,myW){
  N <- length(t.vec)
  # N.fourier <- 1000#floor(N/2) + 1
  freq <- seq(0,0.5, length.out = N.fourier)
  
  delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
  
  ##calculate tapers for this data spacing
  V.mat <- get_tapers(t.vec, W = myW, K = numTapers) #W was 4/N
  
  MTSE_full <- MT_spectralEstimate_freqs(x.t, freq, V.mat$tapers) 
  
  ### calculate the covariance matrix 
  Cov.mat_chave <- matrix(NA, nrow = N.fourier, ncol = N.fourier)
  
  if(calcCov==T){
    for(i in 1:N.fourier){
      if(i %% 100 == 0){print(paste(i," of ",N.fourier))}
      j = 1
      while(j <= i){
        Cov.mat_chave[i,j] <- norm(Conj(t(V.mat$tapers*exp(-im*2*pi*freq[i]*t.vec)*(1/sqrt(numTapers))))%*%(V.mat$tapers*exp(-im*2*pi*freq[j]*t.vec)*(1/sqrt(numTapers))), type = "2") 
        j = j+1
      }
    }
    
    Cov.mat_chave[upper.tri(Cov.mat_chave)] <- t(Cov.mat_chave)[upper.tri(Cov.mat_chave)]
  }
  
  
  ### calc avar and variance of avar
  
  cov.mat=avar=numeric(length(taus))
  
  for(i in 1:length(taus)){
    #taken from above, need to calculate the covariance
    G_tau <- transfer.func(freq,tau = taus[i]) 
    G_tau[1] <- 0 
    
    avar[i]=calcAvar(taus[i],MTSE_full$freqs,MTSE_full$spectrum,delta.f)
    
    #calculate variance for the AVAR estimate at the given tau
    cov.mat[i] <- t(G_tau)%*%(Cov.mat_chave)%*%G_tau*(delta.f)^2#/(sqrt(numTapers)) #can't tell if this term should be here, details below
  }
  
  avarOut=data.frame(tau=taus,avar=avar,var=cov.mat)
  
  ############### analyze old way
  avar.calc <- getAvars(N,x.t, taus = taus)
  oldRes=data.frame(tau=avar.calc$avarRes$taus,avar=avar.calc$avarRes$avars,calculation="old")
  overRes=data.frame(tau=avar.calc$avarRes$taus,avar=avar.calc$avarRes$overavars,calculation="overlapping")
  spectralRes=data.frame(tau=avar.calc$avarRes$taus,avar=avarOut$avar,calculation="spectral")
  allRes=bind_rows(oldRes,overRes,spectralRes)
  
  return(list(x.t=x.t,t.vec=t.vec,V.mat=V.mat,MTSE_full=MTSE_full,
              freq=freq,Cov.mat_chave=Cov.mat_chave,avarOut=avarOut,allAvarRes=allRes,W=myW))
}


######################################################################
######################################################################
## 
######################################################################
######################################################################
x.t=test$r[1,seq(1,2000000,by=1000)]
t.vec <-1:length(x.t)

plot(t.vec,x.t)
acf(x.t)

###below stays the same
N=length(x.t) 
taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]

N.fourier=floor(N/2) + 1
today=format(Sys.Date(),format="%b%d")
resName=paste("pinkNoiseML",N.fourier,today,sep="_")

res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=10,calcCov = F,myW=4/N*2)
res$V.mat$e.values
# saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))


resDF=data.frame(freq=res$MTSE_full$freqs,
                 spectrum=res$MTSE_full$spectrum)
ggplot(resDF,aes(freq,spectrum))+
  geom_point()+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()



ggplot(res$allAvarRes,aes(tau,avar,color=calculation))+
  geom_point()+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()


######################another way
# 
# for(i in 1:300){
#   print(i)
#   set.seed(i)
#   #generate X.t
#   X.t <- TK95(N = 5224, alpha = 1)
#   X.t[c(100:500,1300:1400, 1700:1874, 3000:3450)] <- NA
#   X.t_sims_flk_gps[i,] <- X.t
#   
#   #calculate S.hat
#   MTSE_full <- multitaper_est(X.t, W = 0.009, K = 5)
#   
#   MTSE_saved[[i]] <- MTSE_full
# }

??TK95
library(RobPer)

# x.t=TK95(N = 5224, alpha = 1)
# t.vec <-1:length(x.t)
# 
# plot(t.vec,x.t)
# acf(x.t)
# 
# ###below stays the same
# N=length(x.t) 
# taus <- 2^(0:9)
# taus <- taus[taus<floor(N/3)]
# 
# N.fourier=floor(N/2) + 1
# today=format(Sys.Date(),format="%b%d")
# resName=paste("pinkNoiseML",N.fourier,today,sep="_")
# 
# res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=10,calcCov = F,myW=4/N*2)
# res$V.mat$e.values
# # saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# 
# resDF=data.frame(freq=res$MTSE_full$freqs,
#                  spectrum=res$MTSE_full$spectrum)
# ggplot(resDF,aes(freq,spectrum))+
#   geom_point()+
#   scale_y_log10()+
#   scale_x_log10()+
#   annotation_logticks()
# 
# ggplot(res$allAvarRes,aes(tau,avar,color=calculation))+
#   geom_point()+
#   scale_y_log10()+
#   scale_x_log10()+
#   annotation_logticks()
# 
# 
# keep=res$allAvarRes
# keep$rep=1
# 
# keep_all=data.frame()
# keep_all=bind_rows(keep_all,keep)

keep=keep_all=data.frame()

for(i in 1:5){
  print(i)
  set.seed(i)

  x.t=TK95(N = 5224, alpha = 1)
  t.vec <-1:length(x.t)
  
  ###below stays the same
  N=length(x.t) 
  taus <- 2^(0:12)
  taus <- taus[taus<floor(N/3)]
  
  N.fourier=floor(N/2) + 1
  today=format(Sys.Date(),format="%b%d")

  res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=10,calcCov = F,myW=4/N*2)
  print(res$V.mat$e.values)

  keep=res$allAvarRes
  keep$rep=i
  
  keep_all=bind_rows(keep_all,keep)
  
}


# filter(keep_all,tau!=512)
ggplot(keep_all,aes(tau,avar,color=calculation,group=interaction(tau,calculation)))+
  geom_boxplot()+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()



resDF=data.frame(freq=res$MTSE_full$freqs,
                 spectrum=res$MTSE_full$spectrum,
                 type="MT")
ggplot(resDF,aes(freq,spectrum))+
  geom_point()+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()

library("TSA")
perRes=periodogram(x.t,log="no",plot=TRUE,ylab="Periodogram", xlab="Frequency",lwd=2)

perDF=data.frame(freq=perRes$freq,
                 spectrum=perRes$spec,
                 type="periodogram")

resDF=bind_rows(resDF,perDF)

ggplot(resDF,aes(freq,spectrum,color=type))+
  geom_point()+
  geom_smooth()+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()

# look at old stats papers, see if i can understand why it looks off at low frequencies
# should it be a straight line?
# seems like the spectral estimate is really off for lower freqs, which i think messes up the allan var esti
# would a different spectral estimate work better?
  

# for(i in 1:300){
#   print(i)
#   set.seed(i)
#   #generate X.t
#   X.t <- TK95(N = 5224, alpha = 1)
#   X.t[c(100:500,1300:1400, 1700:1874, 3000:3450)] <- NA
#   X.t_sims_flk_gps[i,] <- X.t
#   
#   #calculate S.hat
#   MTSE_full <- multitaper_est(X.t, W = 0.009, K = 5)
#   
#   MTSE_saved[[i]] <- MTSE_full
# }
