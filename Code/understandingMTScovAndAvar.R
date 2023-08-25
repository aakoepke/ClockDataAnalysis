rm(list=ls())
source("/home/aak3/NIST/ClockDataAnalysis/Code/SA_ImportantFunctions.R")

numberOfSimulations=100 #used when comparing calculated covariance to observed variance in simulated data

######################################################################
### simulate the data
######################################################################
sizeOfData=1000
t.vec <- 1:sizeOfData #time vector
omitted<-c(20:35,50:63, 100:130)
t.vec[omitted] <- NA #take out values
t.vec <- na.omit(t.vec) #vector of times with data
dist.mat <- rdist(t.vec) #distance matrix (delta_nm)

N <- length(t.vec)
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)

delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
numTapers=4

##calculate tapers for this data spacing
V.mat <- get_tapers(t.vec, W = 4/N, K = numTapers)

######################################################################
### using simulated data to check estimate of the spectrum covariance
######################################################################

simOut=data.frame()

for(i in 1:numberOfSimulations){
  ##simulate new dataset with same spacing
  x.t <- rnorm(sizeOfData) #data
  x.t[omitted] <- NA #take out values
  
  ##calculate S.hat
  # test_chave <- multitaper_est(X.t = x.t, W = 4/N, K = 3) #old method, calculates tapers and spectrum
  MTSE_full <- MT_spectralEstimate(x.t, V.mat$tapers) #new function, calculates just spectrum
  ##comparing these 2 estimates
  # plot(MTSE_full$freqs,MTSE_full$spectrum)
  # points(freq,test_chave$spectrum,col="red")
  ## ??? unclear to me why these are different, should not be, maybe small differences in the code?
  # oneSimOut=data.frame(freq=freq,spec.hat=test_chave$spectrum)
  oneSimOut=data.frame(freq=freq,spec.hat=MTSE_full$spectrum)
  
  simOut=bind_rows(simOut,oneSimOut)
}

### calculate the covariance matrix 
Cov.mat_chave <- matrix(NA, nrow = N.fourier, ncol = N.fourier)

for(i in 1:N.fourier){
  j = 1
  while(j <= i){
    Cov.mat_chave[i,j] <- norm(Conj(t(V.mat$tapers*exp(-im*2*pi*freq[i]*t.vec)*(1/sqrt(numTapers))))%*%(V.mat$tapers*exp(-im*2*pi*freq[j]*t.vec)*(1/sqrt(numTapers))), type = "2") 
    ####I added the 1/sqrt(numTapers) factor, seems necessary to get the simulations and the formula to match and makes sense given the Bronez normalization 
    j = j+1
  }
}

Cov.mat_chave[upper.tri(Cov.mat_chave)] <- t(Cov.mat_chave)[upper.tri(Cov.mat_chave)]
# specDF_chave=data.frame(freq=freq,spec.hat=test_chave$spectrum,variance=diag(Cov.mat_chave),type="Chave")
specDF_chave=data.frame(freq=freq,spec.hat=MTSE_full$spectrum,variance=diag(Cov.mat_chave),type="Chave")

simSum=simOut %>% group_by(freq) %>%
  summarise(lowerCI=quantile(spec.hat,.025),upperCI=quantile(spec.hat,.975),mean=mean(spec.hat))

ggplot()+
  geom_point(data = specDF_chave,mapping = aes(freq,spec.hat,col=type))+
  geom_errorbar(data = specDF_chave,aes(freq,spec.hat,ymin=spec.hat-2*sqrt(variance),ymax=spec.hat+2*sqrt(variance),col=type))+
  geom_hline(yintercept = 1)+
  # geom_errorbar(data = simSum, mapping = aes(x = freq,ymin=lowerCI,ymax=upperCI))+
  geom_point(data = simSum, mapping = aes(x = freq,y=mean))
  # geom_point(data = simOut,mapping = aes(freq,spec.hat))

## ??? all the spec.hat i generate are > 0, should my intervals not contain 0? Maybe 2sigma intervals are not appropriate here. 

#############################################################################################
#compare var of simulated points to calc var
#############################################################################################

simSum=simOut %>% group_by(freq) %>% summarise(var=var(spec.hat))
mean(simSum$var)
diag(Cov.mat_chave)

plot(freq,simSum$var)
points(freq,diag(Cov.mat_chave))


#############################################################################################
## do simulation test with avar covariance calculation
#############################################################################################

transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}

## testing the transfer function, looks like the plot in Reappraisal paper
plot(freq,G_tau <- transfer.func(freq,tau = 1),type="l")
points(freq,G_tau <- transfer.func(freq,tau = 4),type="l")
points(freq,G_tau <- transfer.func(freq,tau = 16),type="l")
abline(v=c(1/16,1/8,.125))


calcAvar=function(thetau,freq,spec.hat){
  # G_tau vector length number of frequencies 
  G_tau <- transfer.func(freq,tau = thetau) #change the tau value to get different vectors
  G_tau[1] <- 0 # this was 1 in the old code, but should be 0
  
  avar=G_tau%*%spec.hat*delta.f
  
  return(avar)
}



taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]
# taus<-c(taus,N)

avarOut=data.frame()
for(i in 1:numberOfSimulations){
  x.t <- rnorm(sizeOfData) #data
  x.t[omitted] <- NA #take out values
  
  spec.hat <- MT_spectralEstimate(x.t, V.mat$tapers)
  
  avar=numeric(length(taus))
  for(j in 1:length(taus)){
    avar[j]=calcAvar(taus[j],freq,spec.hat$spectrum)
  }

  avarOut=bind_rows(avarOut,data.frame(tau=taus,avar=avar))
}

avarSimSum=avarOut %>% group_by(tau) %>% summarise(var=var(avar),mean=mean(avar),adev=mean(sqrt(avar)))

##do just one to calculate the covmat
cov.mat=avar=numeric(length(taus))

for(i in 1:length(taus)){
  G_tau <- transfer.func(freq,tau = taus[i]) #change the tau value to get different vectors
  G_tau[1] <- 0 # this was 1 in the old code, but should be 0
  # G_tau vector length number of frequencies 
  
  avar[i]=G_tau%*%MTSE_full$spectrum*delta.f
  
  #calculate variance for the AVAR estimate at the given tau
  # cov.mat.calc[i] <- t(G_tau)%*%(Cov.mat_chave)%*%G_tau*(0.5/256)^2
  cov.mat[i] <- t(G_tau)%*%(Cov.mat_chave)%*%G_tau*(delta.f)^2
  # sample.var <- var(tmat[i,])
  # print(abs(cov.mat.calc-sample.var)/sample.var)
}
sqrt(avar[8])
sqrt(1/N) #doesn't match anymore
round(cov.mat,4)
round(avarSimSum$var,4)



# plot code below taken from SimulationPlots.R
ggplot(avarOut,aes(tau,avar,group=(tau)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  # scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  # scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
    ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau)) 


### calculate avars the old way

oldavars <- data.frame()
for(i in 1:numberOfSimulations){
  set.seed(i)
  print(i)
  #generate X.t
  X.t <- rnorm(N,mean = 0, sd = 1)
  
  # avar.calc <- getAvars(N,X.t, taus = taus[-length(taus)])
  avar.calc <- getAvars(N,X.t, taus = taus)
  temp=data.frame(tau=avar.calc$avarRes$taus,avar=avar.calc$avarRes$avars)
  
  oldavars=bind_rows(oldavars,temp)
}

ggplot(oldavars,aes(tau,avar,group=(tau)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  # scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  # scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau)) 


## compare both

oldavars$method="old"
avarOut$method="spectrum"

allRes=bind_rows(oldavars,avarOut)
ggplot(allRes,aes(tau,avar,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  # scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  # scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau)) 


oldMethodSum=oldavars %>% group_by(tau) %>% summarise(mean=mean(avar),var=var(avar))
simRes=avarOut %>% group_by(tau) %>% summarise(mean=mean(avar),var=var(avar))

round(cov.mat,4)

cov.mat-simRes$var
simRes$var
oldMethodSum$var-simRes$var


allvar=data.frame(tau=rep(taus,3),
           method=rep(c("calculated","SpectralSim","oldSim"),each=length(taus)),
           var=c(cov.mat,simRes$var,oldMethodSum$var))
ggplot(allvar,aes(tau,var,col=method))+
  geom_point()+
  scale_y_log10()+
  scale_x_log10()
  
