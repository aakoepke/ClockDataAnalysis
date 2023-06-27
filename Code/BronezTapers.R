library(RSpectra) #eigensolving
library(geigen)
library(tidyverse)
library(fields)

################ Bronez tapers ################

get_k_eigs=function(Ra,thek,f.w){
  evs=eigs(Ra,k = thek)
  
  return(list("weights" = evs$vectors, "eigenvalues" = evs$values))
}

get_k_geigen=function(Ra,Rb,thek,f.w){
  evs=geigen(Ra,Rb, symmetric = TRUE)
  
  N=dim(evs$vectors)[1]
  return(list("weights" = evs$vectors[,(N-thek + 1):N], "eigenvalues" = sort(evs$values, decreasing = TRUE)[1:thek]))
}


get.weights_bronez <- function(input.list){
  t.n=input.list$t.n 
  K=input.list$K
  f.c=input.list$f.c
  f.w=input.list$f.w 
  dist.mat = input.list$tn_m
  index = input.list$index
  print(index)
  #dist.mat <- rdist(t.n)
    #outer(t.n,t.n,"-") #this option doesn't create a symmetric matrix
  
  # B = [-pi,pi] for omega or [-1/2,1/2] for f
  # R.b <- 1/(pi*(dist.mat))*(sin(dist.mat*pi))
  # R.b[row(R.b) == col(R.b)] <- 1
  # B = identity
  R.b <- diag(1,nrow = length(t.n),ncol = length(t.n))
  
  j = complex(real = 0, imaginary = 1)
  R.a <- exp(j*2*pi*f.c*dist.mat)*(sin(2*pi*f.w*(dist.mat))/(pi*dist.mat))
  R.a[row(R.a) == col(R.a)] <- f.w*2
  
  #Solve the eigenvalue problem
  # out <- tryCatch(get_k_eigs(R.a,K,f.w),error=function(err){get_k_geigen(R.a,R.b,K,f.w)})
  out <- get_k_geigen(R.a,R.b,K,f.w)
  
  return(out)
}

# dim(evs$vectors)
# length(evs$values)
N <- 512
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)

tapers_matlist <- list()
e.vals_list <- list()

# for(i in 1:length(freq)){
#   print(i)
#   temp <- get.weights_bronez(t.n = 1:256, K = 3, f.c = freq[i], f.w = 4/256)
#   #tapers_matlist[[i]] <- temp$weights
#   #e.vals_list[[i]] <- temp$eigenvalues
# }



library(parallel)

input.list=list()
for(i in 1:length(freq)){
  input.list[[i]]=list("t.n"=1:N,
                       "K"=3,
                       "f.c"=freq[i],
                       "f.w"=4/N,
                       "tn_m"=rdist(1:N),
                       "index" = i)
}
## function.to.use takes a one-row dataframe as input,
## and returns a dataframe

startTime = Sys.time()
parResult = mclapply(input.list, get.weights_bronez, mc.cores = 3)
print(Sys.time()-startTime)

#This creates a list where each element corresponds to a frequency

## MTSE for gappy data function 
multitaper_est_bronez <- function(X.t, eig_vecs, K){
  X.t <- X.t - mean(X.t, na.rm = TRUE) #demean
  N.long <- length(X.t)
  t.n <- 1:N.long
  missing.indices <- which(is.na(X.t))
  t.n[which(is.na(X.t))] <- NA

  ##use tapers to generate spectral estimate
  N <- length(na.omit(t.n))
  S.x.hat_MD <- rep(NA, times = floor(N/2) + 1)
  #freqs <- seq(0,0.5, length.out = floor(N/2) + 1)
  
  for(j in 0:floor(N/2)){
    k.vec <- rep(NA,times = K)
    for(k in 0:(K-1)){
      W.t <- eig_vecs[,k+1]*na.exclude(X.t)
      inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*na.omit(t.n)*j/N), na.rm = TRUE)
      k.vec[k + 1] <- abs(inner.sum)^2
    }
    S.x.hat_MD[j+1] <- mean(k.vec)
  }
  
  
  return(list("spectrum" = S.x.hat_MD))
}


#dat=rnorm(N)
#multitaper_est_bronez(dat,parResult[[1]]$weights,3)