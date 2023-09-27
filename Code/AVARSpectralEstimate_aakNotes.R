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
      inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*na.omit(t.n)*j/N), na.rm = TRUE)# question: where does the 1/N come from?
      k.vec[k + 1] <- abs(inner.sum)^2
    }
    S.x.hat_MD[j+1] <- mean(k.vec)
  }
  
  
  return(list("tapers" = eig_vecs, "e.values" = eigdec$values, "spectrum" = S.x.hat_MD))
}

myN=1000

#make a white noise process and add some NA values
myWNvar=1
X.t_missing <- rnorm(myN,mean = 0,sd = sqrt(myWNvar))
# X.t_missing[c(8:20,51:70)] <- NA

#get spectral estimate
spec.est <- multitaper_est(X.t_missing,W=0.01, K = 5)
spec.est$e.values
plot(spec.est$spectrum,type="l") # 

# ############# real dat
# 
# ggplot(alldat,aes(MJD,AlYb))+
#   geom_point()
# 
# spec.est <- multitaper_est(alldat$AlYb[43460:(43460+1000)],W=0.09, K = 5)
# # X.t=filter(alldat,MJD<58190)$AlYb#[43460:(43460+1000)]
# X.t=alldat$AlYb[43460:(43460+1000)]
#############
tapers <- spec.est$tapers

plot(tapers[,1],type = "l")
t(tapers[,2])%*%tapers[,2]
IdenMat=diag(1,myN)
t(tapers[,2])%*% IdenMat %*%tapers[,2]

t.vec <- 1:myN
t.n <- t.vec[which(!is.na(X.t_missing))]
X.t <-  na.omit(X.t_missing)
taus <-  2^(0:6)


#this is where the function that calculates the variance of the
#allan variance estimate, which is in progress:

#var_AVAR_trfunc <- function(tapers, t.n, X.t, taus){

saved <- acf(X.t_missing, na.action = na.pass, lag.max = 1000)
#distance matrix of time points
dist.mat <- rdist(t.n)
#R.x matrix
X.t <- X.t - mean(X.t)
L <- length(X.t)
R.x <- matrix(saved$acf[dist.mat+1], nrow = L, ncol = L)
#f vector
f <- seq(0,0.5, length.out = floor(length(t.n))/2 + 1)

plot(f,spec.est$spectrum)
abline(h=myWNvar)
#for white noise, spectrum is roughly a straight line at the variance of my gaussian data 

#Calculate Cov matrix C
N <- length(f)
Cov.mat <- matrix(NA, nrow = N, ncol = N)
W_matlist <- list()

for(i in 1:N){
  print(i)
W_matlist[[i]] <- t(tapers)%*%exp(-im*2*pi*f[i]*dist.mat)
}
length(W_matlist)

for(i in 1:N){ #going over frequencies
  print(i)
  j = i
  while(j>=i & j <=N){ #going over frequencies
W.star <- Conj(t(W_matlist[[i]]))
    W <- W_matlist[[j]]
    #Cov.mat[i,j] <- sum(abs(W.star%*%R.x%*%W)^2) #frobenius norm 
Cov.mat[i,j] <- sqrt(sum(abs(W%*%diag(x = 1, nrow = myN, ncol = myN)%*%W.star)^2)) #2-norm
    j = j + 1
}
}
Cov.mat[lower.tri(Cov.mat)] <- t(Cov.mat)[lower.tri(Cov.mat)]
which.max(Cov.mat)
norm(W%*%diag(x = 1, nrow = myN, ncol = myN)%*%W.star, type = "F")


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

