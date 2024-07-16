########## Covariance Matrix test ########

#### calculating C_ij = ||V_i*RV_j||^2 with parallelization #####

#### libraries #####
library(future) #parallel
library(future.apply) #parallel
plan(multicore, workers = 12) # can set the number of cores with "workers = X" argument
#plan(multisession, workers = 2)

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
t.vector <- 1:200
dat_length <- length(t.vector)

## tapers
setW = 5/N
setK = 5
taperObject <- get_tapers(t.vector, W = setW, K = setK)
taperMatrix <- taperObject$tapers

## R_mat

R_matrix <- diag(1, nrow = dat_length) #to start

##### uncomment if data is not white noise ######
#acf.lag=4
#s_acf <- stats::acf(X.t, plot=FALSE, lag.max=acf.lag,na.action = stats::na.exclude)$acf

# Create a Toeplitz matrix from the autocorrelation values
#R_mat <- matrix(0, nrow = N, ncol = N)
#R_mat <- stats::toeplitz(c(s_acf, rep(0, N - acf.lag - 1)))

## frequency vector

freq_vec <- seq(0,0.5, length.out = N)

compute_entry <- function(i,j, N = dat_length, K = setK, taperMat = taperMatrix, R_mat = R_matrix, freq = freq_vec, t.vec = t.vector){
  V_star_mat <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  
  V_mat <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
  
  return((1/K)*norm(V_star_mat%*%R_mat%*%V_mat, type = "2"))
  
}


startTime = Sys.time()
print(startTime)
my_list = future.apply::future_mapply(compute_entry, 1:dat_length, rep(1:dat_length, dat_length))
print(Sys.time() - startTime)

