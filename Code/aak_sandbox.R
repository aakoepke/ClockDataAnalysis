library(tidyverse)
library(RSpectra)
library(fields)

setwd("/home/aak3/NIST/ClockDataAnalysis/")  #save this to your file path, or get rid of it if you don't need it

source("Code/readTaraData.R")

clock_df <- dat180403$SrYb
clock_df <- clock_df[min(which(!is.na(clock_df))):length(clock_df)]

###split function to get runs in data between missing values
t.n <- 1:length(clock_df)
t.n <- t.n[-c(which(is.na(clock_df)))]

splits <- split(t.n, cumsum(c(1,diff(t.n) != 1)))
length(splits)
splits[1:10]


count = 0
for(i in 1:length(splits)){
  if(length(splits[[i]])< 10){
    count = count + 1
  }
}
count

#get indices of elements that are part of runs < 10 in length
ind_to_rmv <- unlist(splits[which(lengths(splits)<10)], use.names = FALSE)

clock_df_omitted <- clock_df[-ind_to_rmv]
length(clock_df) - length(clock_df_omitted)
t.n <- 1:length(clock_df_omitted)
t.n <- t.n[-c(which(is.na(clock_df_omitted)))]


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





clock_MTSE_omitted <- multitaper_est(clock_df_omitted, W = 2, K = 3)


# saveRDS(clock_MTSE_omitted,file = "clock_MTSE_omitted.Rds") #had to do this so i could run it on tethys and then play with it locally
# clock_MTSE_omitted=readRDS("clock_MTSE_omitted.Rds")

plot(clock_MTSE_omitted$tapers[,1])
clock_MTSE_omitted$e.values

f <- seq(0, 0.5, length.out = length(clock_MTSE_omitted$spectrum))
# f.all <- seq(0, 0.5, length.out = length(clock_MTSE$spectrum))



plot(log10(f),log10(clock_MTSE_omitted$spectrum))
# plot(log10(f.all),log10(clock_MTSE$spectrum), type = "l", ylab = "log10(S(f))", xlab = "log10(f)")
abline(v = -1.4, lwd = 2, col = "red")
line1.f <- seq(-1.4, log10(0.5), length.out = 100)
a1 = -1.59
b1 = -0.811
lines(line1.f, a1 + b1*line1.f, col = "blue", lwd = 2)
a2 = -0.04
b2 = 0.3
line2.f <- seq(-4.527321, -1.4, length.out = 100)
lines(line2.f, a2 + b2*line2.f, col = "blue", lwd = 2)

# freqs <- log10(f.all[which(log10(f.all)>-1.4)])
# spec.hat <- log10(clock_MTSE$spectrum[which(log10(f.all)> -1.4)])
# lin.fit <- lm(spec.hat ~ freqs)
# 
# freqs1 <- log10(f.all[which(log10(f.all)< -1.4)])[-1]
# spec.hat1 <- log10(clock_MTSE$spectrum[which(log10(f.all)< -1.4)])[-1]
# lin.fit1 <- lm(spec.hat1 ~ freqs1)
# 
# summary(lin.fit1)
# abline(a = -0.04, b = 0.3, col = "blue")
# plot(lin.fit)
# 
