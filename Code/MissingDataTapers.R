######## Missing Data Tapers #######
## In this script we generalize  ###
## the Missing Data Tapers from ####
## Chave 2019 and Bronez 1985    ###
####################################
library(fields)
library(tidyverse)
library(RSpectra)
library(R.matlab)

#The integrated Spectrum estimate is given by: 
# S_hat(f) = (1/K)sum_{k = 0}^{K-1}|sum_{i = 0}^{N-1}v_i^kx_i exp(-i2pift_i)|^2 (eqn 25 of Chave 2019)

#Steps to getting the multitaper estimate:
### 1. choose W
### 2. 2NW tapers will be close to 1, so this selects K
### 3. create estimate, compare to truth if you have it

#############################################
######## MTSE for gappy data function #######
#############################################

multitaper_est <- function(X.t, W, K){
  X.t <- X.t - mean(X.t, na.rm = TRUE) #demean
  N.long <- length(X.t)
  t.n <- 1:N.long
  missing.indices <- which(is.na(X.t))
  t.n[which(is.na(X.t))] <- NA
  
  
  dist.mat <- rdist(na.omit(t.n))
  
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  print("A matrix computed")
  
  eigdec <- eigs(A.prime, k = K, which = "LM")
  eig_vecs <- eigdec$vectors #get only the vectors
  print("tapers computed")
  
  if(K ==1){
    if (mean(Re(eig_vecs))<0){
      eig_vecs <- -eig_vecs
    }
  }
  
  if(K == 2  || K == 3){
    
      if (mean(Re(eig_vecs[,1]))<0){
        eig_vecs[,1] <- -eig_vecs[,1]
      }
      if (Re(eig_vecs[2,2] - eig_vecs[1,2])<0){
        eig_vecs[,2] <- -eig_vecs[,2]
      }
  
    if(K == 3){
      if (mean(Re(eig_vecs[,3]))<0){
        eig_vecs[,3] <- -eig_vecs[,3]
      }
    }
  }
  if(K >=4){
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
  }
  
  print("sign maintenance done")
  
  ##use tapers to generate spectral estimate
  N <- length(na.omit(t.n))
  S.x.hat_MD <- rep(NA, times = floor(N/2) + 1)
  
  for(j in 0:floor(N/2)){
    k.vec <- rep(NA,times = K)
    for(k in 0:(K-1)){
      W.t <- eig_vecs[,k+1]*na.exclude(X.t)
      inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*na.omit(t.n)*j/N), na.rm = TRUE)
      k.vec[k + 1] <- abs(inner.sum)^2
    }
    S.x.hat_MD[j+1] <- mean(k.vec)
  }
  
  
  return(list("tapers" = eig_vecs, "e.values" = eigdec$values, "spectrum" = S.x.hat_MD))
}


#############################################################
######## Check that this equals the Chave, 2019 tapers ######
#############################################################
#generate a time series with the same gaps as the Chave 2019 data for Fig 1 (doesn't have to be the same data)
X.t <- rnorm(14500)
X.t[c(4745:5447, 8378:9545,12823:13051)] <- NA

##set same W as ChaveReceation.R file and ChaveTapers.m matlab file
W <- 0.0007

#Run this
chave.check <- multitaper_est(X.t, W, K = 5)

#these are the tapers/eigenvalues/A matrix from the Matlab code he has
uvec.mat <- readMat("C:/Users/cmb15/OneDrive - UCB-O365/NIST/ClockDataAnalysis/Data/uvec_mat_0.0007.mat")

A.prime.matlab <- uvec.mat$a
lambda.matlab <- uvec.mat$lam
u.matlab <- uvec.mat$uvec

#look at difference between eigenvalues
lambda.matlab
chave.check$e.values

#look at maximum of difference between eigenvectors
max(abs(u.matlab - chave.check$tapers))

plot(log10(seq(0,0.5, length.out = length(chave.check$spectrum))), log10(chave.check$spectrum))




#####check against true spectrum #####

X_process <- rnorm(5000)


plot(spectral.estimate$spectrum)
spec.pgram(X_process)
abline(h = 1, col = "blue")

?arima.sim

X_process <- arima.sim(n = 1000, list(ar = 0.5, sd = 1))


ar.spectrum <- function(omega, phi, sig.e){
  (sig.e^2/(2))*(1/(1 + phi^2 - 2*phi*cos(omega)))
}

spec.pgram(X_process)
lines(seq(0,0.5, length.out = 501), ar.spectrum(seq(0,pi, length.out = 501), phi = 0.5, sig.e = 1), col = "blue")
spectral.estimate <- multitaper_est(X.t = X_process, W = 0.01, K = 5)
plot(seq(0,pi, length.out = 501)/(2*pi), spectral.estimate$spectrum)
