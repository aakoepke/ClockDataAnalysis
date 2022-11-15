##### Chave 2019 ####
library(multitaper)
library(R.matlab)
library(tidyverse)
library(RSpectra)
library(fields)
#time sequence with gaps
t.n <- 1:14500
t.n <- t.n[-c(4745:5447, 8378:9545,12823:13051)]
W <- 0.0007 #0.00097
dist.mat <- rdist(t.n)

A.prime <- sin(2*pi*W*dist.mat)/(pi*dist.mat)
A.prime[row(A.prime) == col(A.prime)] <- W*2
eigdec <- RSpectra::eigs(A.prime, k = 5, which = "LM")

#eigdec <- eigen(A.prime, symmetric = TRUE)

##changing signs of the vectors
K = 5 #number of sequences we want
eig_vecs <- eigdec$vectors[,1:K] #get only those vectors
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




##from matlab:vk.mat, uk.mat
#uk.mat <- readMat("C:/Users/cmb15/OneDrive - UCB-O365/NIST/ClockDataAnalysis/Data/uk.mat")
#dim(uk.mat$u)

uvec.mat <- readMat("C:/Users/cmb15/OneDrive - UCB-O365/NIST/ClockDataAnalysis/Data/uvec_mat_0.0007.mat")

A.prime.matlab <- uvec.mat$a
lambda.matlab <- uvec.mat$lam
u.matlab <- uvec.mat$uvec

### Fig 1
par(mfrow = c(2,1))
colors = c("blue", "red", "green", "magenta", "cyan")
k = 0
s = 1
mdss1 <-eig_vecs[,s]#uk.mat$u[,s]#u.matlab[,s]  #
mdss1_long <- c(mdss1[1:4744],rep(0,times = 5447-4745 + 1), mdss1[4745:(4745+8378-5447-1)],rep(0,times = 9545-8378 + 1), mdss1[7676:(7676 + 12823-9545-1)],rep(0,times = 13051-12823+1),mdss1[10954:length(mdss1)])
#p <- sqrt(sum(mdss1))
plot(mdss1_long,type = "l",ylim = c(-0.03,0.03),col = colors[1]) #sign(sum(mdss1))*

for(i in 2:5){
  mdss1 <- eig_vecs[,i] #uk.mat$u[,i]#u.matlab[,i] 
  mdss1_long <- c(mdss1[1:4744],rep(0,times = 5447-4745 + 1), mdss1[4745:(4745+8378-5447-1)],rep(0,times = 9545-8378 + 1), mdss1[7676:(7676 + 12823-9545-1)],rep(0,times = 13051-12823+1),mdss1[10954:length(mdss1)])
  lines(mdss1_long,type = "l",col = colors[i-5*k])
  
}

max(abs(eig_vecs - u.matlab))



### Fig 3

u.mat <- t(uk.mat$u)
dim(u.mat)

energy.conc <- apply(u.mat^2, FUN = sum, MARGIN = 2)
length(energy.conc)
energy.conc_long <- c(energy.conc[1:4744],rep(0,times = 5447-4745 + 1), energy.conc[4745:(4745+8378-5447-1)],rep(0,times = 9545-8378 + 1), energy.conc[7676:(7676 + 12823-9545-1)],rep(0,times = 13051-12823+1),energy.conc[10954:length(energy.conc)])
plot(energy.conc_long,type = "l",col = colors[1])


### OSS
OSS <- dpss(14500,15,nw=12)
oss.mat <- t(OSS$v)
matplot(OSS$v)
energy.conc.OSS <- apply(oss.mat^2, FUN = sum, MARGIN = 2)
length(energy.conc.OSS)
lines(energy.conc.OSS,type = "l",col = "red")
