
library(future)
plan(multisession)
mtse <- modules::use("Functions.R")

library(Rcpp)

sourceCpp('est_entry.cpp')

# test compute entry cpp

N=5000
tN = 10
R_mat <- toeplitz(c(seq(1,0.1, length.out = tN), rep(0, times = N-tN))) #to start
# Toy example parameters
N.fourier <- N  # Small number for testing

t.vec <- 1:N  # Time vector
freq <- seq(0, .5, length.out = N)  # Frequency vector

setW = 8/N
setK = 5
taperObject <- get_tapers(t.vec, W = setW, K = setK)
taperMat <- taperObject$tapers

V_star_mat <- t(taperMat*exp(1i*2*pi*freq[1]*t.vec))
V_mat <- taperMat*exp(-1i*2*pi*freq[2]*t.vec)

cpp_start_time = Sys.time()
print("Starting Cpp function")
test_entry = est_entry(V_star_mat, V_mat, R_mat, K)
print(Sys.time()-cpp_start_time)

print(test_entry)

C.mat <- matrix(NA, nrow = N.fourier, ncol = N.fourier)


compute_entry <- function(i,j, N = N, K = K, taperMat = taperMat, R_mat = R_mat, freq = freq, t.vec = t.vec){
  V_star_mat <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  V_mat <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
  return((1/K)*norm(V_star_mat%*%R_mat%*%V_mat, type = "2"))
}


if (0) {
start.time=Sys.time()

for(i in 1:N.fourier){
  if(i %% 10 == 0){print(paste(i," of ",N.fourier))}
  j = 1
  V_star_mat <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  while(j <= i){
    V_mat <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
    (1/K)*norm(V_star_mat%*%R_mat%*%V_mat, type = "2") 
    j = j+1
  }
}
C.mat[upper.tri(C.mat)] <- t(C.mat)[upper.tri(C.mat)]
print(Sys.time()-start.time)
}

#Check the timing for each of the operations in the standard loop

# 1. Compute the V_star_mat
start.time=Sys.time()
t(taperMat*exp(1i*2*pi*freq[1]*t.vec))
print(Sys.time()-start.time)

# 2. Compute the V_mat
start.time=Sys.time()
taperMat*exp(-1i*2*pi*freq[2000]*t.vec)
print(Sys.time()-start.time)

# 3. Compute the norm
start.time=Sys.time()
norm(V_star_mat%*%R_mat%*%V_mat, type = "2")
print(Sys.time()-start.time)

if (0) {

my_list <- list()

k = 1

for(i in 1:N.fourier){
  if(i %% 10 == 0){print(paste(i," of ",N.fourier))}
  j = 1
  while(j <= i){
    my_list[[k]] <-  future({ compute_entry(i,j, N = N, K = K, taperMat = taperMat, R_mat = R_mat, freq = freq, t.vec = t.vec) })
    j = j+1
    k = k+1
  }
}

fast.start.time=Sys.time()
my_list = lapply(my_list, value)
print(Sys.time()-fast.start.time)

print(my_list)

}