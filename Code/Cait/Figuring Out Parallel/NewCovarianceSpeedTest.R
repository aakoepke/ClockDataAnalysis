########## New Covariance in C++ vs. Old R code speed comparison #########
setwd("C:/Users/cmb15/OneDrive - UCB-O365/NIST/ClockDataAnalysis/Code/Cait")
library(future) #parallel
library(future.apply) #parallel
#plan(multicore, workers = 12) # can set the number of cores with "workers = X" argument
plan(multisession, workers = 8)

mtse <- modules::use("Functions.R")

compute_entry_slow <- function(i,j, taperMat = taperMatrix, setK = K, R_mat = R_matrix){
  #message("Worker PID: ", Sys.getpid())
  V_star_mat <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  V_mat <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
  return((1/setK)*norm(V_star_mat%*%R_mat%*%V_mat, type = "2"))
}

compute_entry <- function(i,j, taperMat = taperMatrix, setK = K, setN = N, c = c_vec){
  V_star <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  V <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
  return(est_entry_FFT(V_star, V, c, setK, setN))
}

library(Rcpp)
library(fields)
library(RSpectra)

sourceCpp('est_entry_FFT.cpp')

# test compute entry cpp

N=400
t.vec <- 1:N  # Time vector

## tapers
setW = 8/N
K = 5
taperObject <- mtse$get_tapers(t.vec, W = setW, K = setK)
taperMatrix <- taperObject$tapers

###c vector
max.lag.acf = 10
sample.acf <- seq(1,0.1, length.out = max.lag.acf) #toy example
c_vec <- c(sample.acf,rep(0, times = N-max.lag.acf), 0, rep(0, times = N-max.lag.acf), rev(sample.acf[-1]))
length(c_vec)


## frequency vector
freq <- seq(0, 0.5, length.out = N)  # Frequency vector

tN = 10
R_matrix <- toeplitz(c(seq(1,0.1, length.out = tN), rep(0, times = N-tN))) #to start




# Test a single computation


start_time_fast = Sys.time()

compute_entry(1,10)

total_time_fast = Sys.time() - start_time_fast

start_time_slow = Sys.time()
compute_entry_slow(1,10)
total_time_slow = Sys.time() - start_time_slow

as.numeric(total_time_slow) / as.numeric(total_time_fast)




### future
startTime = Sys.time()
print(startTime)
my_list = future.apply::future_mapply(compute_entry, 1:N, rep(1:N, N)) #
print(Sys.time() - startTime)

# future kinda sucks bro, so we are gonna like try some foreach stuff? is that sigma or skibidi?

### no future
k = 1
startTime = Sys.time()
print(startTime)
for(i in 1:N){
  for(j in 1:N){
    my_list[[k]] = compute_entry_slow(i,j)
    k = k+1
  }
}
print(Sys.time() - startTime)







