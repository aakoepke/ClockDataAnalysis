### setwd
setwd("/home/cmb15/ClockDataAnalysis/Code/Paper1")

source("../SA_ImportantFunctions.R")
################## libraries ###############################

library(tidyverse)

############################################################


#source("/home/cmb15/ClockDataAnalysis/Code/Paper1/FlickerNoise_noGaps.R")

numberOfSimulations = 1000
N = 2048

## keeping track of how long this all takes
startTime=Sys.time()

## saving the date to label file outputs
runDate=format(Sys.Date(),"%m%d%y")

### add in determination of W and K for this data pattern?
###run1
# setWnum = 3
# setW = setWnum/N
# setK = 3

###run2
# setWnum = 3
# setW = setWnum/N
# setK = 4

###run3
# setWnum = 3
# setW = setWnum/N
# setK = 5

###run4
# setWnum = 4
# setW = setWnum/N
# setK = 3

###run5
# setWnum = 4
# setW = setWnum/N
# setK = 4

###run6
setWnum = 6
setW = setWnum/N
setK = 8

V.mat <- get_tapers(1:N, W = setW, K = setK)
taperMatrix <- V.mat$tapers
V.mat$e.values
plot(taperMatrix[,setK])
###run7
# setWnum = 4
# setW = setWnum/N
# setK = 7


    
print(setWnum)
print(setK)

#######################################
###### Study 2a: Flicker Noise ########
######      with no gaps     ##########
#######################################

##initialize vectors/matrices to save results into

#avar estimates using transfer function, bandpass variance, and LS periodogram
trfunc.vec <- bpvar.vec <- trfunc.vec_LS <- rep(NA, times = numberOfSimulations) #1 x numberofSimulations

#values of tau to calculate sigma^2(tau) 
taus <- c(2^(0:9), floor(N/3))

#matrices for saving transfer function, bandpass variance, and ls periodogram based estimates
tmat <- bmat <- lmat <- matrix(NA, nrow = numberOfSimulations, ncol = length(taus)) #number of Simulations x number of taus

#matrix for saving the numberofSimulations MTSEs
MTSE_mat <- matrix(NA, ncol = N/2 + 1, nrow = numberOfSimulations)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2] #delta of frequencies

##calculate tapers
V.mat <- get_tapers(1:N, W = setW, K = setK)
taperMatrix <- V.mat$tapers
V.mat$e.values
plot(taperMatrix[,6])

for(i in 1:numberOfSimulations){
  print(i)
  set.seed(i)
  
  #generate X.t
  X.t <- generate_pink_noise(length = N, fs = 1) #arfima.sim(N, model = list(dfrac = 0.4999999))
  #X.t_missing[c(100:150,1300:1400, 1700:1752)] <- NA
  
  #calculate S.hat, save to matrix row
  MTSE_mat[i,] <- MT_spectralEstimate_fft(X.t, taperMatrix)$spectrum
  #LS periodogram
  #lsperio <- lomb_scargle(X.t_missing, f[-1])
  
}

g = 0
for(k in taus){
  print(k)
  g = g + 1
  tau = k
  
  for(i in 1:numberOfSimulations){
    
    #calculate bandpass variance
    temp_bp <- integrate(approxfun(f, MTSE_mat[i,]), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 1000)
    bpvar.vec[i] <- 4*temp_bp$value
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 0
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_mat[i,])
    
  }
  
  tmat[,g] <- trfunc.vec
  bmat[,g] <- bpvar.vec
}

##AVAR calculation
### only need to run once #####
# 
amat <- oamat <- matrix(NA, nrow = numberOfSimulations, ncol = length(taus))

for(i in 1:numberOfSimulations){
  print(i)
  set.seed(i)
  X.t <- generate_pink_noise(length = N, fs = 1) #arfima.sim(N, model = list(dfrac = 0.4999999))

  avar.calc <- getAvars(N,X.t, taus = taus)
  amat[i,] <- avar.calc$avarRes$avars
  oamat[i,] <- avar.calc$avarRes$overavars
}

# likely need to save tmat and bmat to work with outside of the titans
saveRDS(tmat,paste("Code/Paper1/Results/tmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_FlickerNoiseNoGaps.Rds",sep=""))
saveRDS(bmat,paste("Code/Paper1/Results/bmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_FlickerNoiseNoGaps.Rds",sep=""))
saveRDS(amat,paste("Code/Paper1/Results/amat",runDate,"_N",N,"_",numberOfSimulations,"sims_FlickerNoiseNoGaps.Rds",sep=""))
saveRDS(oamat,paste("Code/Paper1/Results/oamat",runDate,"_N",N,"_",numberOfSimulations,"sims_FlickerNoiseNoGaps.Rds",sep=""))

# }

