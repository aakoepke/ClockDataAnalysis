########### What's the deal? ###########
## C_ij = ||W*(f_j)R_xW(f_i)||^2
library(fields)
library(RSpectra)

### Look at Bronez method

t.vec <- 1:60 #vector of time points
t.vec[c(3:8,14, 18, 25, 33:37)] <- NA
t.vec <- na.omit(t.vec)

dist.mat <- rdist(t.vec)
W <- 0.19
A <- sin(2*pi*W*dist.mat)/(pi*dist.mat)
A[col(A) == row(A)] <- 2*W  #A matrix

eig.vecs <- eigs(A, k = 6, which = "LM")
eig.vecs <- eigen(A, symmetric = TRUE)

eig.vecs$values

vectors <- eig.vecs$vectors

dim(vectors)

dim(vectors%*%Conj(t(vectors)))
dim(Conj(t(vectors))%*%vectors)

Conj(t(vectors))%*%vectors

vectors%*%Conj(t(vectors))

eigen(matrix(c(1,1,2,1), ncol = 2))

eigen(matrix(2*c(1,1,2,1), ncol = 2))


#### Covariance calculation #####

t.vec <- 1:300
t.vec[c(20:35,50:63, 100:130)] <- NA
t.vec <- na.omit(t.vec)
dist.mat <- rdist(t.vec)
W <- 0.04

##Generate A
A <- sin(2*pi*W*dist.mat)/(pi*dist.mat)
A[col(A) == row(A)] <- 2*W  #A matrix
  
##Solve eigenvalue problem for A
eig.vecs <- eigs(A, k = 6, which = "LM")

##Get V, V* matrix
V <- eig.vecs$vectors

V_star <- Conj(t(V))

##R_x matrix (identity for white noise)


##Covariance
N <- length(na.omit(t.vec))
fourier.length <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N)
Cov.mat <- matrix(NA, ncol = N, nrow = N)

for(i in 1:length(freq)){
  print(i)
  j = 1
  while(j <= i){
    #Cov.mat[i,j] <- exp(freq[i] + freq[j])*norm(V_star%*%V, type = "2") - exp(2*freq[i] + 2*freq[j])*sum(diag(V_star%*%V))
    Cov.mat[i,j] <- norm(V_star%*%V, type = "2") - sum(diag(V_star%*%V))
    j = j+1
  }
}

Cov.mat[upper.tri(Cov.mat)] <- t(Cov.mat)[upper.tri(Cov.mat)]
norm(Cov.mat, type= "2")

G_tau <- transfer.func(freq,8)
G_tau[1] <- 1

t(G_tau)%*%Cov.mat%*%G_tau*(freq[2])^2



X_t <- rnorm(1000)

test <- multitaper_est(X_t, W = 0.09, K = 5)
plot(test$tapers[,1])
plot(test$spectrum)

test2 <- multitaper_est(X_t, W = 0.01, K = 5)
plot(test2$tapers[,1])
plot(test2$spectrum, type = "l")

test$e.values
test2$e.values

###### eigen values of 4 #####
myWNvar=1
X.t_missing <- rnorm(myN,mean = 0,sd = sqrt(myWNvar))
# X.t_missing[c(8:20,51:70)] <- NA

#get spectral estimate
spec.est <- multitaper_est(X.t_missing,W=0.09, K = 5)

