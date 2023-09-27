rm(list=ls())
source("/home/aak3/NIST/ClockDataAnalysis/Code/SA_ImportantFunctions.R")

numberOfSimulations=1000 #used when comparing calculated covariance to observed variance in simulated data

######################################################################
### simulate the data
######################################################################
sizeOfData=100
t.vec <- 1:sizeOfData #time vector
omitted<-c(20:35,50:63, 100:130)
t.vec[omitted] <- NA #take out values
t.vec <- na.omit(t.vec) #vector of times with data

N <- length(t.vec)
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)

delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
numTapers=3

##calculate tapers for this data spacing
V.mat <- get_tapers(t.vec, W = 4/N, K = numTapers)

### calculate the covariance matrix; this will be the same for all data with this spacing
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

#############################################################################################
## simulation test with avar covariance calculation
#############################################################################################

transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}

calcAvarAndVar=function(thetau,freq,spec.hat,delta.f,Cov.mat){
  # G_tau vector length number of frequencies 
  G_tau <- transfer.func(freq,tau = thetau) #change the tau value to get different vectors
  G_tau[1] <- 0 # this was 1 in the old code, but should be 0
  
  avar=G_tau%*%spec.hat*delta.f
  
  avar.var <- t(G_tau)%*%(Cov.mat)%*%G_tau*(delta.f)^2#/(sqrt(numTapers)) #can't tell if this term should be here, details below
  
  return(data.frame(tau=thetau,avar=avar,avar.var=avar.var))
}


taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]

avarOut=data.frame()
for(i in 1:numberOfSimulations){
  x.t <- rnorm(sizeOfData) #data
  x.t[omitted] <- NA #take out values
  
  spec.hat <- MT_spectralEstimate(x.t, V.mat$tapers)
  
  avarTemp=avars=data.frame()
  
  for(j in 1:length(taus)){
    avarTemp=calcAvarAndVar(taus[j],freq,spec.hat$spectrum,delta.f,Cov.mat_chave)
    avars=bind_rows(avars,avarTemp)
  }

  avarOut=bind_rows(avarOut,avars)
}

avarSimSum=avarOut %>% group_by(tau) %>% 
  summarise(mean=mean(avar),adev=mean(sqrt(avar)),simVar=var(avar),var=mean(avar.var)) %>% 
  mutate(varRatio=var/simVar)

avarSimSum

sqrt(numTapers)

### Our calculated variance for avar seems to be off by a factor of sqrt(numTapers), ass seen in avarSimSum$varRatio, 
### but it gets closer as tau increases. I can't explain this. 
### If I add a 1/sqrt(numTapers) to the cov.mat for the avar, the variance ratio is much closer to 1, except for higher tau
### Is there intuition for why the spectral avar would be more variable for higher tau than expected? Like with the old method (effectively less data)

