########### What's the deal? ###########
## C_ij = ||W*(f_j)R_xW(f_i)||^2

source("Code/SA_ImportantFunctions.R")

im <- complex(real = 0, imaginary = 1)

#### Covariance calculation #####
###gaps####
t.vec <- 1:300 #time vector
X.t <- rnorm(300) #data
t.vec[c(20:35,50:63, 100:130)] <- NA #take out values
X.t[c(20:35,50:63, 100:130)] <- NA #take out values
t.vec <- na.omit(t.vec) #vector of times with data
dist.mat <- rdist(t.vec) #distance matrix (delta_nm)

#### no gaps #######
t.vec <- 1:300 #time vector
X.t <- rnorm(300) #data
dist.mat <- rdist(t.vec) #distance matrix (delta_nm)

##get multitaper spectral estimate:
MTSE_test <- multitaper_est(X.t, W = 4/300, K = 5) #calculate tapers


##Get V, V* matrix
V <- MTSE_test$tapers
dim(V) #columns of V are the eigenvectors, dimension 239 X 5
V_star <- Conj(t(V)) #dim: 5 x 239, should be real valued so the conjugate isn't really necessary here
dim(V_star)



A.size <- 4/300
numTapers <- 5


##Covariance
N <- length(na.omit(t.vec))
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)
Cov.mat <- matrix(NA, ncol = N.fourier, nrow = N.fourier)

for(i in 1:N.fourier){
  print(i)
  j = 1
  while(j <= i){
    Cov.mat[i,j] <- norm(V_star%*%exp(im*2*pi*freq[i]*dist.mat)%*%exp(-im*2*pi*freq[j]*dist.mat)%*%V, type = "2")*(A.size/numTapers)^2
    j = j+1
  }
}





############# treating exp like a vector #############
N = 2048
t.vec <- 1:N #time vector
X.t <- rnorm(N) #data
dist.mat <- rdist(t.vec) #distance matrix (delta_nm)

##get multitaper spectral estimate:
W = 4/N
K = 5
MTSE_test_tapers <- get_tapers(X.t, W = W, K = K) #calculate tapers
MTSE_test_spec <- MT_spectralEstimate(X.t = X.t, V.mat = MTSE_test_tapers)

t(V[,1])%*%V[,1]
##Get V, V* matrix
V <- MTSE_test_tapers*(sqrt(2*W/K))

sum(abs(t(V[,1]*exp(im*2*pi*freq[10]*1:N))%*%(V[,1]*exp(-im*2*pi*freq[129]*1:N)))^2)

N <- dim(V)[1]
time.vec <- dist.mat[1,]
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)
Cov.mat <- matrix(NA, ncol = N.fourier, nrow = N.fourier)

for(i in 1:N.fourier){
V.exp.i <- V*exp(-im*2*pi*freq[i]*1:N)
  print(i)
  j = 1
  while(j <= i){
    V.exp.j <- V*exp(-im*2*pi*freq[j]*1:N)
    Cov.mat[i,j] <- norm(Conj(t(V.exp.i))%*%V.exp.j, type = "2") #norm(V_star%*%exp(im*2*pi*freq[i]*dist.mat)%*%exp(-im*2*pi*freq[j]*dist.mat)%*%V, type = "2")*(A.size/numTapers)^2
    j = j+1
  }
}


Cov.mat[upper.tri(Cov.mat)] <- t(Cov.mat)[upper.tri(Cov.mat)]
size.cov <- norm(Cov.mat)
Cov.mat[1:10,1:10]
#Calculate G_tau vector

G_tau <- transfer.func(freq,tau = 1) #change the tau value to get different vectors
G_tau[1] <- 1
#calculate variance for the AVAR estimate at the given tau
t(G_tau)%*%(Cov.mat)%*%G_tau*(freq[2])^2

s##to delete
tmat <- readRDS(file = "Code/Paper1/Results/tmat050523_W4_K6_N2048_300sims_WhiteNoiseNoGaps.Rds") #readRDS(file = "Code/Paper1/Results/bmat053123_W5_K9_N2048_300sims_FlickerNoiseNoGaps.Rds") #
tmat %<>% t()
dim(tmat) #each column is 300 tau estimates
taus <- c(2^(0:9), floor(2048/3))

f <- seq(0,0.5, length.out = 2048/2 + 1)
#Calculate G_tau vector
i  = 2
G_tau <- transfer.func(f,tau = taus[i]) #change the tau value to get different vectors
G_tau[1] <- 1

t(G_tau)%*%G_tau*(f[2])^2
var(tmat[,1])

i = j = 1
exp(im*2*pi*freq[i]*dist.mat)%*%exp(-im*2*pi*freq[j]*dist.mat)
dist.mat
freq[1]
exp(im*2*pi*freq[i]*dist.mat)
exp(-im*2*pi*freq[j]*dist.mat)
