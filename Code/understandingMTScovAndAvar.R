rm(list=ls())
source("/home/aak3/NIST/ClockDataAnalysis/Code/SA_ImportantFunctions.R")

numberOfSimulations=100 #used when comparing calculated covariance to observed variance in simulated data

######################################################################
### simulate the data
######################################################################
sizeOfData=100
t.vec <- 1:sizeOfData #time vector
omitted<-c(20:35,50:63, 100:130)
t.vec[omitted] <- NA #take out values
t.vec <- na.omit(t.vec) #vector of times with data
dist.mat <- rdist(t.vec) #distance matrix (delta_nm)

N <- length(t.vec)
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)

delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
numTapers=3

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

### data frame with last spectral estimate and its covariance
specDF_chave=data.frame(freq=freq,spec.hat=MTSE_full$spectrum,variance=diag(Cov.mat_chave),type="Chave")

# simSum=simOut %>% group_by(freq) %>%
#   summarise(mean=mean(spec.hat),variance=var(spec.hat))

# ggplot()+
#   # geom_point(data = specDF_chave,mapping = aes(freq,spec.hat,col=type))+
#   # geom_errorbar(data = specDF_chave,aes(freq,spec.hat,ymin=spec.hat-sqrt(variance),ymax=spec.hat+sqrt(variance),col=type))+
#   geom_hline(yintercept = 1)+
#   geom_errorbar(data = simSum, mapping = aes(x = freq,ymin=mean-sqrt(var.Chave),ymax=mean+sqrt(var.Chave),col="Chave"))+
#   geom_errorbar(data = simSum, mapping = aes(x = freq,ymin=mean-sqrt(variance),ymax=mean+sqrt(variance)))+
#   geom_point(data = simSum, mapping = aes(x = freq,y=mean))
#   # geom_point(data = simOut,mapping = aes(freq,spec.hat))

#############################################################################################
#compare var of simulated points to calc var
#############################################################################################

simSum=simOut %>% group_by(freq) %>% summarise(var=var(spec.hat))
simSum$var.Chave=diag(Cov.mat_chave)

mean(simSum$var)
diag(Cov.mat_chave)

plot(freq,simSum$var)
points(freq,diag(Cov.mat_chave),col="red")

simSum=simSum %>% 
mutate(varianceDifference=var.Chave-var)

ggplot(simSum, aes(freq,varianceDifference))+
geom_point()

#### These tests indicated to me that there is a (1/sqrt(numTapers)) term missing from our covariance calculation notation

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
points(freq,G_tau <- transfer.func(freq,tau = N),type="l")
abline(v=c(1/16,1/8,.125))


calcAvar=function(thetau,freq,spec.hat,delta.f){
# G_tau vector length number of frequencies 
G_tau <- transfer.func(freq,tau = thetau) #change the tau value to get different vectors
G_tau[1] <- 0 # this was 1 in the old code, but should be 0

avar=G_tau%*%spec.hat*delta.f

return(avar)
}


taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]


avarOut=data.frame()
for(i in 1:numberOfSimulations){
x.t <- rnorm(sizeOfData) #data
x.t[omitted] <- NA #take out values

spec.hat <- MT_spectralEstimate(x.t, V.mat$tapers)

avar=numeric(length(taus))
for(j in 1:length(taus)){
  avar[j]=calcAvar(taus[j],freq,spec.hat$spectrum,delta.f)
}

avarOut=bind_rows(avarOut,data.frame(tau=taus,avar=avar))
}

avarSimSum=avarOut %>% group_by(tau) %>% summarise(var=var(avar),mean=mean(avar),adev=mean(sqrt(avar)))

##do just one to calculate the covmat
cov.mat=avar=numeric(length(taus))

for(i in 1:length(taus)){
#taken from above, need to calculate the covariance
G_tau <- transfer.func(freq,tau = taus[i]) 
G_tau[1] <- 0 

avar[i]=calcAvar(taus[i],freq,spec.hat$spectrum,delta.f)

#calculate variance for the AVAR estimate at the given tau
# cov.mat.calc[i] <- t(G_tau)%*%(Cov.mat_chave)%*%G_tau*(0.5/256)^2
cov.mat[i] <- t(G_tau)%*%(Cov.mat_chave)%*%G_tau*(delta.f)^2#/(sqrt(numTapers)) #can't tell if this term should be here, details below
# sample.var <- var(tmat[i,])
# print(abs(cov.mat.calc-sample.var)/sample.var)
}


###################################################################################
### Does the calculated variance match the simulated variance?
###################################################################################

round(cov.mat,4)
round(avarSimSum$var,4)
cov.mat/avarSimSum$var
sqrt(numTapers)
## these numbers are close, not quite the same, close enough? is there a (1/sqrt(numTapers)) missing from the calculation?

### plot code below taken from SimulationPlots.R
### This plot looks at the distribution of the simulated data spectral avar estimates
ggplot(avarOut,aes(tau,avar,group=(tau)))+
geom_boxplot(lwd = 1.2)+
### add true straight line below
geom_abline(slope = -1,intercept = 0,size=1)+
theme(legend.position = c(.15, .2))+
scale_y_log10()+
scale_x_log10()+
annotation_logticks()+
ylab(expression(sigma^2*(tau)))+
xlab(expression(tau)) 
### looks like they follow expected the line well

################################
### calculate avars the old way
################################

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

### This plot looks at the distribution of the simulated data avar estimates via old method
ggplot(oldavars,aes(tau,avar,group=(tau)))+
geom_boxplot(lwd = 1.2)+
### add true straight line below
geom_abline(slope = -1,intercept = 0,size=1)+
theme(legend.position = c(.15, .2))+
scale_y_log10()+
scale_x_log10()+
annotation_logticks()+
ylab(expression(sigma^2*(tau)))+
xlab(expression(tau)) 


### plot them both to see how they compare

oldavars$method="old"
avarOut$method="spectrum"

allRes=bind_rows(oldavars,avarOut)
ggplot(allRes,aes(tau,avar,col=method,group=interaction(tau,method)))+
geom_boxplot(lwd = 1.2)+
### add true straight line below
geom_abline(slope = -1,intercept = 0,size=1)+
scale_y_log10()+
scale_x_log10()+
annotation_logticks()+
ylab(expression(sigma^2*(tau)))+
xlab(expression(tau)) 

### results look similar, spectral looks tighter especially for higher tau

oldMethodSum=oldavars %>% group_by(tau) %>% summarise(mean=mean(avar),var=var(avar))
simRes=avarOut %>% group_by(tau) %>% summarise(mean=mean(avar),var=var(avar))

allvar=data.frame(tau=rep(taus,3),
     method=rep(c("calculated","SpectralSim","oldSim"),each=length(taus)),
     var=c(cov.mat,simRes$var,oldMethodSum$var))
ggplot(allvar,aes(tau,var,col=method))+
geom_point()+
scale_y_log10()+
scale_x_log10()

##############################################################################
### outstanding questions
##############################################################################

##############################################################################
### is the avar covariance right?

round(cov.mat,4)

cov.mat/simRes$var
sqrt(numTapers)
### seems to be off by a factor of sqrt(numTapers) (?), but gets closer as tau increases? maybe? 
### I can't explain this or tell if this is real. After looking at it with higher numberOfSimulations
### and sizeOfData, I'm not sure it's related to K. 
### If I add a 1/sqrt(numTapers) to the cov.mat for the avar, the variance ratio is much closer to 1, except for higher tau
### Is there intuition for why the spectral avar would be more variable for higher tau than expected? Like with the old method (effectively less data)
oldMethodSum$var-simRes$var

###############################################################################
# CAIT ignore below here, just playing around
###############################################################################
### can I calculate avar for tau =N?
## testing the transfer function, looks like the plot in Reappraisal paper
plot(freq,G_tau <- transfer.func(freq,tau = 1),type="l")
points(freq,G_tau <- transfer.func(freq,tau = N),type="l")
points(freq,G_tau <- transfer.func(freq,tau = N/1.5),type="l")
points(freq,G_tau <- transfer.func(freq,tau = N/2),type="l")

## this is my problem, G_tau is essentially 0 here

G_tau <- transfer.func(freq,tau = N) #change the tau value to get different vectors
G_tau[1] <- 0

avar_N=G_tau%*%spec.hat$spectrum

sqrt(avar_N)
sqrt(1/N) 

### this seems to work, removing the delta.f. I have no clue why. Tried for N=100 and 1000, both match

calcAvar(N,freq = freq,spec.hat = spec.hat$spectrum,delta.f = delta.f)
calcAvar(N-1,freq = freq,spec.hat = spec.hat$spectrum,delta.f = delta.f)
calcAvar(N/1.5,freq = freq,spec.hat = spec.hat$spectrum,delta.f = delta.f)
1/N

### Is there a minimum tau that this seems to work with?


##############################################################################
##############################################################################
### bandpass variance
##############################################################################
##############################################################################

#calculate bandpass variance
tau=4
temp_bp <- integrate(approxfun(freq, MTSE_full$spectrum), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 10000)
bpvar.vec <- 4*temp_bp$value



getBPvar=function(freq,spectrum, tau){
  temp_bp <- integrate(approxfun(freq, spectrum), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 1000)
  bpvar.vec <- 4*temp_bp$value
}
  
bpvars <- data.frame()
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
