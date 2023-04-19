######### Multitaper Spectral Estimation of Power Law Processes ##############

### This script checks the variance of the multitaper spectral estimate of a power law process



######## Example #1 #############

## White noise with variance = 1

#set parameters
sig.e <- 1

#generate 1000 simulations of the process
sims.mat <- matrix(NA, nrow = 1000, ncol = 100)
for(i in 1:1000){
  set.seed(i)
  sims.mat[i,] <- rnorm(100)
}

#Calculate S.hat for each of the 1000 simulations
###periodogram

spec.periodogram <-matrix(NA, nrow = 1000, ncol = 50)
spec.multitaper <-matrix(NA, nrow = 1000, ncol = 51)
for(i in 1:1000){
  spec.estimate <- spec.pgram(sims.mat[i,], plot = FALSE)
  spec.periodogram[i,] <- spec.estimate$spec
}

###MT spectral estimator
for(i in 1:1000){
  spec.estimate <- multitaper_est(sims.mat[i,], W = 0.09, K = 5)
  spec.multitaper[i,] <- spec.estimate$spectrum
}


#look at distribution of S.hat(f)

par(mfrow = c(2,1))
hist(spec.periodogram[,1], breaks = 100, freq = FALSE)
hist(spec.multitaper[,1], breaks = 100)

hist(spec.periodogram[,10], breaks = 100, freq = FALSE)
hist(spec.multitaper[,10], breaks = 100, freq = FALSE)
spec.estimate <- spec.pgram(sims.mat[1,], plot = FALSE)
lines(seq(0,8,by = 0.1),3*dchisq(seq(0,8,by = 0.1), df = 1)/2)



######## Example #2 #############

## a fractional difference process X_t with d = 0.25 (so stationary) and sigma^2_epsilon = 1

#set parameters
d <- 0.25
sig.e <- 1

#generate 1000 simulations of the process


#Calculate S.hat for each of the 1000 simulations


#look at distribution of S.hat(f)







