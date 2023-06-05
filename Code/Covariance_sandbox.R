########### What's the deal? ###########
## C_ij = ||W*(f_j)R_xW(f_i)||^2

source("Code/SA_ImportantFunctions.R")

im <- complex(real = 0, imaginary = 1)

#### Covariance calculation #####

t.vec <- 1:300 #time vector
X.t <- rnorm(300) #data
t.vec[c(20:35,50:63, 100:130)] <- NA #take out values
X.t[c(20:35,50:63, 100:130)] <- NA #take out values
t.vec <- na.omit(t.vec) #vector of times with data
dist.mat <- rdist(t.vec) #distance matrix (delta_nm)

##get multitaper spectral estimate:
MTSE_test <- multitaper_est(X.t, W = 4/239, K = 5) #calculate tapers


##Get V, V* matrix
V <- MTSE_test$tapers
dim(V) #columns of V are the eigenvectors, dimension 239 X 5
V_star <- Conj(t(V)) #dim: 5 x 239, should be real valued so the conjugate isn't really necessary here
dim(V_star)

##R_x matrix (identity for white noise)


##Covariance
N <- length(na.omit(t.vec))
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)
Cov.mat <- matrix(NA, ncol = N.fourier, nrow = N.fourier)

for(i in 1:N.fourier){
  print(i)
  j = 1
  while(j <= i){
    Cov.mat[i,j] <- norm(V_star%*%exp(im*2*pi*freq[i]*dist.mat)%*%exp(-im*2*pi*freq[j]*dist.mat)%*%V, type = "2")
    j = j+1
  }
}

Cov.mat[upper.tri(Cov.mat)] <- t(Cov.mat)[upper.tri(Cov.mat)]

#Calculate G_tau vector
G_tau <- transfer.func(freq,tau = 50) #change the tau value to get different vectors
G_tau[1] <- 1
#calculate variance for the AVAR estimate at the given tau
t(G_tau)%*%Cov.mat%*%G_tau*(freq[2])^2





