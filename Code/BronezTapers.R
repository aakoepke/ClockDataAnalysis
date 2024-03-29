library(RSpectra) #eigensolving
library(geigen)
library(tidyverse)
library(fields)
library(parallel)

################ Bronez tapers ################

# get_k_eigs=function(Ra,thek,f.w){
#   evs=eigs(Ra,k = thek)
#   
#   return(list("weights" = evs$vectors, "eigenvalues" = evs$values))
# }

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

N <- 256
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)


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

# saveRDS(parResult,"/home/aak3/NIST/ClockDataAnalysis/Code/BronezResults256_070423.Rds")