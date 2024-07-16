########## Covariance Matrix test ########

#### calculating C_ij = ||V_i*RV_j||^2 with parallelization #####

#### libraries #####
library(future) #parallel
library(future.apply) #parallel
plan(multicore) # can set the number of cores with "workers = X" argument

### needed functions ####
get_tapers <- function(t.n, W, K){
  dist.mat <- fields::rdist(stats::na.omit(t.n))
  
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  print("A matrix computed")
  
  eigdec <- RSpectra::eigs(A.prime, k = K, which = "LM")
  eig_vecs <- eigdec$vectors #get only the vectors
  print("tapers computed")
  
  return(list("tapers" = eig_vecs, "e.values" = eigdec$values))
} 


## time vector
t.vec <- 1:1000
N <- length(t.vec)

## tapers
W = 5/N
K = 5
taperObject <- get_tapers(t.vec, W = W, K = K)
taperMat <- taperObject$tapers

## R_mat

R_mat <- diag(1, nrow = N) #to start

##### uncomment if data is not white noise ######
#acf.lag=4
#s_acf <- stats::acf(X.t, plot=FALSE, lag.max=acf.lag,na.action = stats::na.exclude)$acf

# Create a Toeplitz matrix from the autocorrelation values
#R_mat <- matrix(0, nrow = N, ncol = N)
#R_mat <- stats::toeplitz(c(s_acf, rep(0, N - acf.lag - 1)))

## frequency vector

freq <- seq(0,0.5, length.out = N)

compute_entry <- function(i,j, N = N, K = K, taperMat = taperMat, R_mat = R_mat, freq = freq, t.vec = t.vec){
  
  V_star_mat <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  
  V_mat <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
  
  return((1/K)*norm(V_star_mat%*%R_mat%*%V_mat, type = "2"))
  
}



my_list = future.apply::future_mapply(compute_entry, 1:N, rep(1:N, N))


