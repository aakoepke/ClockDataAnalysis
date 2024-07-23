############################## Parallelization ##################################
##### Title: Parallel_Covariance_MacLinux.R                                 #####
##### Description: This file contains code to parallelize the calculation   #####
#####              of the covariance matrix of the spectral estimate using  #####
#####              c++ function and mcapply() for use on Mac/Linux machines #####
#################################################################################-

###################-
#### libraries ####
###################-

library(parallel)
library(Rcpp)
library(RcppArmadillo)

###################-
#### functions ####
###################-

## upper triangle indices
generate_upper_triangle_indices <- function(N) {
  # Create a sequence of row indices
  row_indices <- rep(1:N, each = N)
  # Create a sequence of column indices
  col_indices <- rep(1:N, times = N)
  # Combine row and column indices into a matrix
  all_indices <- cbind(row_indices, col_indices)
  # Filter out indices to get only the upper triangle (excluding the diagonal)
  upper_triangle_indices <- all_indices[row_indices < col_indices, ]
  return(upper_triangle_indices)
}

## compute the entry C.mat_ij
compute_entry_parallel <- function(ij, taperMat = taperMatrix, setK = K, setN = N, c = c_vec){
  #i,j, taperMat = taperMatrix, setK = K, setN = N, c = c_vec
  i = ij[1]
  j = ij[2]
  
  V_star <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  V <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
  return(est_entry_FFT(V_star, V, c, setK, setN))
}

## fill in the upper triangle of a matrix from a vector
fill_upper_triangle <- function(vec, N, incl_diag = TRUE) {
  # Initialize an NxN matrix with NA
  mat <- matrix(NA, nrow = N, ncol = N)
  
  # Get the indices of the upper triangle (excluding the diagonal)
  upper_tri_indices1 <- which(upper.tri(mat, diag = incl_diag), arr.ind = TRUE)
  upper_tri_indices <- upper_tri_indices1[order(upper_tri_indices1[,1]),]
  
  # Check if the length of the vector matches the number of upper triangle elements
  if (length(vec) != nrow(upper_tri_indices)) {
    stop("Length of vector does not match the number of upper triangle elements excluding the diagonal.")
  }
  
  # Fill the matrix upper triangle with the elements of the vector
  mat[upper_tri_indices] <- vec
  
  return(mat)
}

matrix_to_list <- function(mat) {
  # Convert each row of the matrix to an element of the list
  row_list <- split(mat, row(mat))
  
  return(row_list)
}


setwd("/home/cmb15/ClockDataAnalysis/Code/Cait")
## mtse module
mtse <- modules::use("Functions.R")

###################-
#####  Method #####
###################-

sourceCpp('est_entry_FFT.cpp')


## Parallelization ####

### pick number of cores 
### can use detectCores() to see how many are available
numCores <- 16

### 1. Calculate predetermined variables ####

N=1000 # length of data
t.vec <- 1:N  # time vector

#### tapers
setW = 8/N
K = 5
taperObject <- mtse$get_tapers(t.vec, W = setW, K = K)
taperMatrix <- taperObject$tapers

#### c vector
max.lag.acf = 10
sample.acf <- seq(1,0.1, length.out = max.lag.acf) #toy example
c_vec <- c(sample.acf,rep(0, times = N-max.lag.acf), 0, rep(0, times = N-max.lag.acf), rev(sample.acf[-1]))
length(c_vec)

#### frequency vector
freq <- seq(0, 0.5, length.out = N)  # Frequency vector

#### list of indices 
upper_triangle_indices <- generate_upper_triangle_indices(N)
row_list <- matrix_to_list(as.matrix(upper_triangle_indices))

## Run code on multiple cores ####

##### (run this code chunk all at once ###
#####   for most accurate timing)      ###
  start_time_fast = Sys.time() #for timing, start time
  print("started mclapply")
  my_list <- vector(mode = "list", length = length(row_list)) #list() #a vector to hold the entries
  my_list <- mclapply(row_list,compute_entry_parallel,  mc.cores = numCores)
  
  total_time_fast = Sys.time() - start_time_fast
  total_time_fast 
  
#######--------------------------------###
#   
# ## Make the Covariance matrix ####
# ### 1. Fill in upper triangle of C.mat ####
# C.mat <- fill_upper_triangle(vec = my_list, N, incl_diag = FALSE)
#   
# ### 2. Fill in lower triangle of C.mat ####
# C.mat[lower.tri(C.mat)] <- t(C.mat)[lower.tri(C.mat)]
# 
# ### 3. Fill in Diagonal of C.mat ####
# tN = max.lag.acf
# R_mat <- toeplitz(c(seq(1,0.1, length.out = tN), rep(0, times = N-tN))) #to start
# diag(C.mat) <- norm(t(taperMatrix)%*%R_matrix%*%taperMatrix/K, type = "2")
# 
# ## look at result ##
# C.mat[1:10,1:10]


