# source("/home/aak3/NIST/ClockDataAnalysis/Code/analyzingRealData.R")
######################################################################
### Read in data (from readTaraData.R)
######################################################################
rm(list=ls())

library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library("openxlsx")

###############################################################################################
### 04/03/2018 analyzed on 07/10/18 frac.diff. = (measured-ideal)/ideal * 1e15
### Columns: MJD  Al vs Yb frac.diff. from(2.1628871275166585754)  Sr vs Yb frac.diff. from(1.2075070393433381221)  Al vs Sr frac.diff. from(2.6117014317814574243)
###############################################################################################
dat180403 <- read_delim("/home/aak3/NIST/ClockDataAnalysis/Data/180403 optical analysis TiS_fixed.dat",
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE, skip = 2)
colnames(dat180403)=c("MJD","AlYb","SrYb","AlSr")
dat180403long=melt(dat180403,id.vars = "MJD",variable.name = "Ratio",value.name = "FracDiff")
dat180403long=dat180403long%>%mutate(missing=is.na(FracDiff),
                                     date=convertToDateTime(MJD, origin = "1858-11-17",tz="MDT"))
dat180403long$seconds=1:dim(dat180403long)[1]
######################################################################

#####################make a function
source("/home/aak3/NIST/ClockDataAnalysis/Code/SA_ImportantFunctions.R")

### calc avar and unc
transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}

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

startTime=Sys.time()
# ######################################################################
# ## AlYb
# ######################################################################
# AlYbdat040318=filter(dat180403long,Ratio=="AlYb" & missing==F)
# t.vec <- AlYbdat040318$seconds
# x.t=AlYbdat040318$FracDiff
# resName="AlYb040318_500"
# 
# ###below stays the same
# N=length(x.t)
# taus <- 2^(0:9)
# taus <- taus[taus<floor(N/3)]
# # floor(N/2) + 1
# res=calculateAvars(x.t,t.vec,taus,N.fourier=500,numTapers=3,calcCov = T)
# saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# ######################################################################
# ## AlSr
# ######################################################################
# AlSrdat040318=filter(dat180403long,Ratio=="AlSr" & missing==F)
# t.vec <- AlSrdat040318$seconds
# x.t=AlSrdat040318$FracDiff
# resName="AlSr040318"
# 
# ###below stays the same
# N=length(x.t)
# taus <- 2^(0:9)
# taus <- taus[taus<floor(N/3)]
# 
# res=calculateAvars(x.t,t.vec,taus,N.fourier=floor(N/2) + 1,numTapers=3,calcCov = T)
# saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# ######################################################################
# ## SrYb
# ######################################################################
# SrYbdat040318=filter(dat180403long,Ratio=="SrYb" & missing==F)
# t.vec <- SrYbdat040318$seconds
# x.t=SrYbdat040318$FracDiff
# resName="SrYb040318"
# 
# ###below stays the same
# N=length(x.t)
# taus <- 2^(0:9)
# taus <- taus[taus<floor(N/3)]
# 
# res=calculateAvars(x.t,t.vec,taus,N.fourier=floor(N/2) + 1,numTapers=3,calcCov = T)
# saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# Sys.time()-startTime
# 
# 
# 
# # #####################look at res
# # resName="AlSr040318"
resName="SrYb040318"
resName="AlYb040318res"

res=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
res$V.mat$e.values

plot(log10(res$MTSE_full$freqs),log10(res$MTSE_full$spectrum),type="l")

### plot avar
ggplot(res$avarOut,aes(tau,avar,ymin=avar-var,ymax=avar+var))+
  geom_point()+
  geom_errorbar()+
  ### add true straight line below
  geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))

ggplot(res$allAvarRes,aes(tau,avar,col=calculation))+
  geom_point()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))
# 
# 
# 
# 
# 

######################################################################
## simdat
######################################################################
# AlSrdat040318=filter(dat180403long,Ratio=="AlSr" & missing==F)
# t.vec <- AlSrdat040318$seconds
# x.t=rnorm(length(t.vec))
# resName="simDat"
# 
# ###below stays the same
# N=length(x.t)
# taus <- 2^(0:9)
# taus <- taus[taus<floor(N/3)]
# # floor(N/2) + 1
# today=format(Sys.Date(),format="%b%d")
# 
# N.fourier=200
# resName=paste("simDat",N.fourier,today,sep="_")
# freq <- seq(0,0.5, length.out = N.fourier)
# res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=3,calcCov = T,myW = freq[2]/2)
# saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# N.fourier=500
# resName=paste("simDat",N.fourier,today,sep="_")
# freq <- seq(0,0.5, length.out = N.fourier)
# res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=3,calcCov = T,myW = freq[2]/2)
# saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# N.fourier=1000
# resName=paste("simDat",N.fourier,today,sep="_")
# freq <- seq(0,0.5, length.out = N.fourier)
# res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=3,calcCov = T,myW = freq[2])
# saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# N.fourier=2000
# resName=paste("simDat",N.fourier,today,sep="_")
# freq <- seq(0,0.5, length.out = N.fourier)
# res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=3,calcCov = T,myW = freq[2]*2)
# saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# 
# # change width of analysis band
# # look at tapers and eig vals, how change
# 


#########################simulate real data

######Noise Generation####
library(RobPer)

# #White Noise
# y <- TK95(N = 2000, alpha = 0)
# t <- seq(along = y)
# plot(t,y,type = "l", main = "white noise", ylab = "y", xlab = "t")
# 
# #Flicker (Pink) Noise
# y <- TK95(N = 2000, alpha = 1)
# t <- seq(along = y)
# plot(t,y,type = "l", main = "flicker noise", ylab = "y", xlab = "t")
# 
# #Random Walk
# y<- TK95(N = 2000, alpha = 2)
# t <- seq(along = y)
# plot(t,y,type = "l", main = "random walk noise", ylab = "y", xlab = "t")



AlSrdat040318=filter(dat180403long,Ratio=="AlSr" & seconds <210000 & seconds>200001)
t.vec <- filter(AlSrdat040318,missing==F)$seconds
plot(AlSrdat040318$seconds,AlSrdat040318$FracDiff)

allN=length(t.vec[1]:t.vec[length(t.vec)])

########white noise
x.t=TK95(N = allN,alpha=0)
x.t=x.t[AlSrdat040318$missing==F]
simType="whiteNoiseSimTK95"
plot(t.vec,x.t)

###below stays the same
N=length(x.t)
taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]
# floor(N/2) + 1
today=format(Sys.Date(),format="%b%d")

N.fourier=100
resName=paste(simType,N.fourier,today,sep="_")
freq <- seq(0,0.5, length.out = N.fourier)
res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=3,calcCov = T,myW = freq[2])
res$V.mat$e.values
saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))

plot(res$freq,res$MTSE_full$spectrum)

#########Flicker (Pink) Noise
x.t=TK95(N = allN,alpha=1)
x.t=x.t[AlSrdat040318$missing==F]
simType="flickerNoiseSimTK95"
plot(t.vec,x.t)

###below stays the same
N=length(x.t)
taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]
# floor(N/2) + 1
today=format(Sys.Date(),format="%b%d")

N.fourier=100
resName=paste(simType,N.fourier,today,sep="_")
freq <- seq(0,0.5, length.out = N.fourier)
res=calculateAvars(x.t,t.vec,taus,N.fourier=N.fourier,numTapers=3,calcCov = T,myW = freq[2])
res$V.mat$e.values
saveRDS(res,file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))

plot(log10(res$freq),log10(res$MTSE_full$spectrum))

