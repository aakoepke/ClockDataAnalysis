#######################################################################################
##### Calculate sigma.hat_tau (the Allan Variance) using the Spectral Estimate ########
#######################################################################################

###################
#### libraries ####
###################

library(tidyverse)
library(RSpectra)
library(fields)

###########################
#### things we'll need ####
###########################

## i
im <- complex(real = 0, imaginary = 1)

## MTSE for gappy data function 
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


######################################################################################################
###### Calculate sigma.hat_tau using transfer function method (see "reappraisal...", page 70) ########
######################################################################################################

#G_tau(f) function
transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}

###########   Allan Variance Calculation    #############
######## (eq'n at bottom of p. 70 of reappraisal) #######
#input: spectral estimate (as a vector), taus (as a vector) where you want the AVAR calculated
#output: a vector of the AVAR estimates

AVAR_trfunc <- function(spectral_est, taus){
out <- rep(NA, times = length(taus))
f <- seq(0,0.5,length.out = length(spectral_est))

for(i in 1:length(taus)){
G.vec <- transfer.func(f, taus[i])
G.vec[1] <- 1
out[i] <- f[2]*sum(G.vec*spectral_est)
}

return(out)
}



######## Calculating the variance of this estimate of the AVAR ########


### to test out function below, run the following lines:

#make a white noise process and add some NA values
X.t_missing <- rnorm(1000, mean = 0, sd = 1)
X.t_missing[c(8:20,51:70, 200:220)] <- NA

#get spectral estimate
spec.est <- multitaper_est(X.t_missing,W=0.0125, K = 5) #0.09
tapers <- spec.est$tapers
t.vec <- 1:300
t.n <- t.vec[which(!is.na(X.t_missing))]
X.t <-  na.omit(X.t_missing)
taus <-  2^(0:6)
plot(tapers[,1])

#this is where the function that calculates the variance of the
#allan variance estimate, which is in progress:

#var_AVAR_trfunc <- function(tapers, t.n, X.t, taus){

saved <- acf(X.t_missing, na.action = na.pass, lag.max = 100)
#distance matrix of time points
dist.mat <- rdist(t.n)
#R.x matrix
X.t <- X.t - mean(X.t)
L <- length(X.t)
R.x <- matrix(saved$acf[dist.mat+1], nrow = L, ncol = L)
#f vector
f <- seq(0,0.5, length.out = floor(length(t.n))/2 + 1)


#Calculate Cov matrix C
N <- length(f)
Cov.mat <- matrix(NA, nrow = N, ncol = N)
W_matlist <- list()

for(i in 1:N){
  print(i)
  W_matlist[[i]] <- t(tapers)%*%exp(-im*2*pi*f[i]*dist.mat)
}

for(i in 1:N){
  print(i)
  j = i
  while(j>=i & j <=N){
    W.star <- t(Conj(W_matlist[[i]]))
    W <- W_matlist[[j]]
    #Cov.mat[i,j] <- sum(abs(W.star%*%R.x%*%W)^2) #frobenius norm 
    Cov.mat[i,j] <- sqrt(sum(abs(W.star%*%diag(x = 1, nrow = length(X.t), ncol = length(X.t))%*%W)^2)) #2-norm
    j = j + 1
  }
}

Cov.mat[lower.tri(Cov.mat)] <- t(Cov.mat)[lower.tri(Cov.mat)]

out <- rep(NA, times = length(taus))
#for(k in 1:length(taus)){
k = 1
#Calculate transfer function vector G_tau
G_tau <- transfer.func(f,taus[k])
G_tau[1] <- 1
#G_tau^T * C * G_tau
out[k] <- t(G_tau)%*%Cov.mat%*%G_tau 
#}
return(out)
#}





#Calculate sigma.hat_tau using bandpass variance

#input: spectral estimate, taus where you want the AVAR calculated
#output: a vector of the AVAR estimates
AVAR_bpvar <- function(spectral_est, taus){
  out <- rep(NA, times = length(taus))
  f <- seq(0,0.5,length.out = length(spectral_est))
  delta.f <- f[2]
  
  for(i in 1:length(taus)){
  tau <- taus[i]
  if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
    f.min.index <- which(f == 1/(4*tau))
    f.max.index <- which(f == 1/(2*tau))
    out[i] <- 4*delta.f*(sum(spectral_est[f.min.index:(f.max.index - 1)]))
  }
  else{
    f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
    f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
    if(f.min.index == f.max.index){
      out[i] <- 4*(spectral_est[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + spectral_est[f.min.index]*(1/(2*tau) - f[f.min.index]))
    }
    else{
      out[i] <- 4*delta.f*(sum(spectral_est[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
    }
  }
}

  return(out)
}


#Variance of sigma.hat_tau using transfer function

var_AVAR_trfunc <- function(tapers, t.n, X.t, taus){
  #distance matrix of time points
  dist.mat <- rdist(t.n)
  #R.x matrix
  X.t <- X.t - mean(X.t)
  L <- length(X.t)
  R.x <- matrix(NA, nrow = L, ncol = L)
 
  
  #f vector
  f <- seq(0,0.5, length.out = floor(length(t.n))/2 + 1)
 
  
  #Calculate Cov matrix C
  N <- dim(tapers)[1]
  
  for(i in 1:N){
    for(j in 1:N){
      W.star <- Conj(t(tapers)%*%exp(-im*2*pi*f[i]*dist.mat))
      W <- t(tapers)%*%exp(-im*2*pi*f[j]*dist.mat)
      C[i,j] <- norm(W.star%*%R.x[i,j]%*%W, type = "F")
    }
  }
    
for(k in 1:length(taus)){
#Calculate transfer function vector G_tau
  G_tau <- transfer.func(f,tau)

#G_tau^T * C * G_tau
  out[k] <- t(G_tau)%*%C%*%G_tau 
}
  return(out)
}

#Variances of sigma.hat_tau using bandpass variance

