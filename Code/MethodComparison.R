########## Method Comparison ###############
## This script compares the traditional   ##
## AVAR methods of analysis with spectral ##
## approach and includes some EDA tools   ##
############################################
library(purrr) #discrete uniform dist.
library(tidyverse)
library(RSpectra)
#create data
N <-  2048
X.t <- X.t_missing <-  rnorm(N, mean = 0, sd = 3)
#create data with gaps

###random gaps#####
###set.seed(20)
###endpoints <- sort(rdunif(4,2048,1))

###placed gaps####
endpoints <- c(450,700,1000,1200)


X.t_missing[c(endpoints[1]:endpoints[2], endpoints[3]:endpoints[4])] <- NA


#look at the data
par(mfrow = c(2,1))
plot(X.t, type = "l")
plot(X.t_missing, type = "l")


#### EDA on X.t and X.t_missing ####

#histogram
hist(X.t)
hist(X.t_missing)

#still looks normally distributed
acf(X.t)
acf(na.omit(X.t_missing))



############# Two Sample Variance Methods  ####################
### **see steps in google drive file and add here when ready** ####

#calculate AVAR of X.t without any missing data
AVAR.WN <- getAvars(N,X.t)

plot(log10(AVAR.WN$avarRes$taus), log10(AVAR.WN$avarRes$avars), ylab = "log10(AVAR)", xlab = "log10(tau)",pch = 19 )

#calculate AVAR of X.t with missing gaps, pushed together
AVAR.WN_missing <- getAvars(length(na.omit(X.t_missing)),na.omit(X.t_missing)) #1435 data points

points(log10(AVAR.WN_missing$avarRes$taus), log10(AVAR.WN_missing$avarRes$avars), col = "red", xlab = "log10(AVAR)", ylab = "log10(tau)",pch = 19  )
abline(a = 1, b = -1, col = "blue")


############# Multitaper Spectral Estimate (MTSE)  ###############
### **see steps in google drive file and add here when ready** ####

###########################################################
##### calculate MTSE for X.t without any missing data #####
###########################################################

f.nyquist <- 1/2
N.prime <- N
f.j.prime <- seq(0,f.nyquist, length.out = N)
X.tilde.prime <- X.t


#function to create taper h'_k,N,t
h.k.prime <- function(k,N,N.prime){
  t <- 0:(N-1)
  h.k.N <- (2/(N+1))^(1/2)*sin((k + 1)*pi*(t + 1)/(N + 1)) #value of h' taper for t = 0, .., N-1
  if(N.prime > N){
    return(c(h.k.N, rep(0, times = N.prime-N)))#add in the 0's to pad the end to get up to 4096
  }
  else{
    return(h.k.N)
  }
}

plot(h.k.prime(k = 0, N = 4000, N.prime = 4096), type = "l") #look at taper for k = 0 with zero padding
plot(h.k.prime(k = 1, N = 4000, N.prime = 4096), type = "l") #look at taper for k = 1 with zero padding


#building multitaper S_x(f) estimator
t <- 0:(N.prime-1)
S.x.hat <- rep(NA, times = N.prime) #where S_x.hat(f'_j) will be saved

for(j in 0:(N.prime-1)){
  k.vec <- rep(NA,times = 6)
  for(k in 0:5){
    W.t <- h.k.prime(k = k, N=N, N.prime = N.prime)*X.tilde.prime
    inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*t*j/N.prime))
    k.vec[k + 1] <- abs(inner.sum)^2
  }
  S.x.hat[j+1] <- (1/length(k.vec))*sum(k.vec)
}

plot(log10(f.j.prime), log10(S.x.hat), type = "l") #plot of multitaper spectral estimator
abline(h = 0)



##Calculate \int_{-f_N}^{f_N} S.x.hat(f)df to get a variance estimate
delta.f <- f.j.prime[2] - f.j.prime[1]

sqrt(2*delta.f*sum(na.omit(S.x.hat)))


###########################################################
#####     calculate MTSE for X.t with missing data    #####
###########################################################

#create missing data tapers

t.n <- 1:N
t.n[c(endpoints[1]:endpoints[2], endpoints[3]:endpoints[4])] <- NA
NW <- 7
W <- NW/length(t.n)
dist.mat <- rdist(t.n)
K = 6 #number of sequences we want

#create the A' matrix (Chave 2019 equation (22))
A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
A.prime[row(A.prime) == col(A.prime)] <- W*2
eigdec <- eigs_sym(A.prime, k = K, which = "LM")


eig_vecs <- eigdec$vectors #get only the vectors

#some sign maintenance
    for(i in seq(1,K,by = 2)){
      if (mean(Re(eig_vecs[,i]))<0){
        eig_vecs[,i] <- -eig_vecs[,i]
      }
    }
    
    for(i in seq(2,K-1,by = 2)){
      if (Re(eig_vecs[2,i] - eig_vecs[1,i])<0){
        eig_vecs[,i] <- -eig_vecs[,i]
      }
    }

##plot the tapers
colors = c("blue", "red", "green", "magenta", "cyan")
k = 0
s = 1
mdss <- eig_vecs[,s]
mdss_long <- X.t_missing #this contains NAs in the correct spots
mdss_long[!is.na(mdss_long)] <- mdss

plot(mdss_long,type = "l",col = colors[1]) #sign(sum(mdss1))*

for(i in 2:6){
  mdss <- eig_vecs[,i] #uk.mat$u[,i]#
  mdss_long <- X.t_missing #this contains NAs in the correct spots
  mdss_long[!is.na(mdss_long)] <- mdss
  lines(mdss_long,type = "l",col = colors[i-5*k])
}


##use tapers to generate spectral estimate
N.short <- floor(length(t.n)/2)
S.x.hat_MD <- rep(NA, times = N)

for(j in 0:(N.prime-1)){
  k.vec <- rep(NA,times = 6)
  for(k in 0:5){
    v_k <- insert(eig_vecs[,k+1], ats=c(rep(endpoints[1], times = endpoints[2] - endpoints[1] + 1),rep(endpoints[3], times = endpoints[4] - endpoints[3] + 1)))
    W.t <- v_k*X.t_missing
    inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*t.n*j/N.prime), na.rm = TRUE)
    k.vec[k + 1] <- abs(inner.sum)^2
  }
  S.x.hat_MD[j+1] <- mean(k.vec)
}

freqs <- seq(0,f.nyquist, length.out = N)
lines(log10(freqs), log10(S.x.hat_MD), type = "l", col = "red") #plot of multitaper spectral estimator
plot(freqs, S.x.hat_MD, type = "l") #plot of multitaper spectral estimator

S.x.hat_MD <- rep(NA, times = N)

#for(j in 0:(N-1)){
j = 0
  k.vec <- rep(NA,times = 6)
  #for(k in 0:5){
    k = 0
    v_k <- insert(eig_vecs[,k+1], ats=c(rep(endpoints[1], times = endpoints[2] - endpoints[1] + 1),rep(endpoints[3], times = endpoints[4] - endpoints[3] + 1)))
    W.t <- v_k*X.t_missing
    inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*t.n*j/N), na.rm = TRUE)
    k.vec[k + 1] <- abs(inner.sum)^2
  #}
  S.x.hat_MD[j+1] <- mean(k.vec)
#}




############# MTSE for different W and K ##############





  
  
############# Clock Data Analysis  ##################
load("Data/ClockData.RData")
plot(clock_df$time,clock_df$AlYb)

hist(na.omit(clock_df$AlYb))
acf(na.omit(clock_df$AlYb)) #definitely not white noise







