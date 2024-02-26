########### AVAR estimate comparisons for first differenced data ##########



### setwd
setwd("/home/cmb15/ClockDataAnalysis/Code/Paper1")

source("Code/SA_ImportantFunctions.R")
##### functions ####

WNFN_theoreticalAVAR <- function(taus, e_white, flicker_floor){
  ### calculate theoretical allan variance of white noise + flicker process ###
  #taus = vector of taus you would like calculated
  #e_white = standard deviation of the white noise process
  #flicker_floor = the flicker floor
  
  ##white noise AVAR + flicker floor
  out <- e_white^2/taus + flicker_floor
  return(out)
}

################## libraries ###############################

library(tidyverse)

############################################################
#missings for 1000 missing: c(100:500,1300:1400, 1700:1874, 2400:2722)
#missings for 205 missing: c(100:150,1300:1400, 1700:1752)
#missings for 410 missing: c(100:201,1300:1501, 1700:1805)
#2662:  c(100:300,1300:1552, 1700:1805, 1913:1966)

#source("/home/cmb15/ClockDataAnalysis/Code/Paper1/FlickerNoise_noGaps.R")

#1024 missing: c(100:400, 750:864 ,1300:1600, 1700:1905, 2013:2113)
#614 missing: c(100:300,1300:1552, 1700:1805, 1913:1966)
#410 missing: c(100:201,1300:1501, 1700:1805)
#205 missing: c(100:150,1300:1400, 1700:1752)

numberOfSimulations = 1000
N.long = 2048 + 205
t.n_missing <- 1:N.long
missing.window <- c(100:150,1300:1400, 1700:1752)
t.n_missing[missing.window] <- NA
N <- length(na.omit(t.n_missing))
N
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
setWnum = 12
setW = setWnum/N
setK = 12

##calculate tapers
V.mat <- get_tapers(t.n_missing, W = setW, K = setK)
taperMatrix <- V.mat$tapers
V.mat$e.values

###run7
# setWnum = 4
# setW = setWnum/N
# setK = 7



print(setWnum)
print(setK)

##initialize vectors/matrices to save results into

#avar estimates using transfer function, bandpass variance, and LS periodogram
trfunc.vec <- bpvar.vec <- trfunc.vec_LS <- rep(NA, times = numberOfSimulations) #1 x numberofSimulations

#values of tau to calculate sigma^2(tau) 
taus <- c(2^(0:9), floor(N/3), 1000)

#matrices for saving transfer function, bandpass variance, and ls periodogram based estimates
tmat <- bmat <- lmat <- matrix(NA, nrow = numberOfSimulations, ncol = length(taus)) #number of Simulations x number of taus

#matrix for saving the numberofSimulations MTSEs
MTSE_mat <- matrix(NA, ncol = N/2 + 1, nrow = numberOfSimulations)
perio_mat <- matrix(NA, ncol = N/2 + 1, nrow = numberOfSimulations)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2] #delta of frequencies




for(i in 1:numberOfSimulations){
  print(i)
  set.seed(i)
  
  X.t_missing <- rnorm(N.long) + 0.001*TK95(N = N.long, alpha = 1) 
  X.t_missing[missing.window] <- NA
  
  #calculate S.hat, save to matrix row
  MTSE_mat[i,] <- MT_spectralEstimate(X.t_missing, taperMatrix)$spectrum
  #periodogram (because there are no gaps)
  #perio_mat[i,2:(N/2+1)] <- spec.pgram(combo, plot = FALSE)$spec
  
}

g = 0
for(k in taus){
  print(k)
  g = g + 1
  tau = k
  
  for(i in 1:numberOfSimulations){
    
    #calculate bandpass variance
    #temp_bp <- integrate(approxfun(f, MTSE_mat[i,]), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 1000)
    #bpvar.vec[i] <- 4*temp_bp$value
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 0
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_mat[i,])
    
  }
  
  tmat[,g] <- trfunc.vec
  #bmat[,g] <- bpvar.vec
}

##AVAR calculation
### only need to run once #####
# 
amat <- oamat <- matrix(NA, nrow = numberOfSimulations, ncol = length(taus))

for(i in 1:numberOfSimulations){
  print(i)
  set.seed(i)
  
  X.t_missing <- rnorm(N.long) + 0.001*TK95(N = N.long, alpha = 1) 
  X.t_missing[missing.window] <- NA
  
  avar.calc <- getAvars(N.long,na.omit(X.t_missing), taus = taus)
  amat[i,] <- avar.calc$avarRes$avars
  oamat[i,] <- avar.calc$avarRes$overavars
}

# likely need to save tmat and bmat to work with outside of the titans
saveRDS(tmat,paste("Code/Paper1/Results/tmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_ComboNoiseWithGaps_10per.Rds",sep=""))
saveRDS(oamat,paste("Code/Paper1/Results/oamat",runDate,"_N",N,"_",numberOfSimulations,"sims_ComboNoiseWithGaps_10per.Rds",sep=""))

# }

