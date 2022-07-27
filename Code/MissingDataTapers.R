######## Missing Data Tapers #######
## In this script we generalize  ###
## the Missing Data Tapers from ####
## Chave 2019 and Bronez 1985    ###
####################################


#The integrated Spectrum estimate is given by: 
# S_hat(f) = (1/K)sum_{k = 0}^{K-1}|sum_{i = 0}^{N-1}v_i^kx_i exp(-i2pift_i)|^2 (eqn 25 of Chave 2019)

#Steps to getting the multitaper estimate:
### 1. choose W
### 2. 2NW tapers will be close to 1, so this selects K
### 3. create estimate, compare to truth if you have it

########## Example 1: White Noise  #####################
N <- 2048
X.t <- rnorm(N)
