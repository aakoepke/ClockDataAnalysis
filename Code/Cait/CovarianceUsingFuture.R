########## Covariance Matrix test ########

#### calculating C_ij = ||V_i*RV_j||^2 with parallelization #####

#### libraries #####
library(future)
library(future.apply)
plan(multicore) # can set the number of cores with "workers = X" argument


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


