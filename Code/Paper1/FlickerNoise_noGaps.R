### setwd
setwd("/home/cmb15/ClockDataAnalysis/Code/Paper1")

source("../SA_ImportantFunctions.R")
################## libraries ###############################

library(tidyverse)

############################################################


#source("/home/cmb15/ClockDataAnalysis/Code/Paper1/FlickerNoise_noGaps.R")

numberOfSimulations = 300
N <- 2048
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
# setWnum = 4
# setW = setWnum/N
# setK = 5

###run7
# setWnum = 4
# setW = setWnum/N
# setK = 7

Wnums <- 6:10
Ks <- 2*Wnums -1


for(w in 1:length(Wnums)){
    setWnum = Wnums[w]
    setW = setWnum/N
    setK = Ks[w]
    
print(setWnum)
print(setK)

#######################################
###### Study 2a: Flicker Noise ########
######      with no gaps     ##########
#######################################

trfunc.vec <- bpvar.vec <- rep(NA, times = numberOfSimulations)

tmat <- bmat <- matrix(NA, ncol = numberOfSimulations, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]


X.t_sims_flk <- readRDS("Results/FlickerSims_noGaps_N2048_300.Rds")

g = 0
for(k in c(2^(0:9), floor(N/3))){
  g = g + 1
  tau = k

  for(i in 1:numberOfSimulations){
    print(i)
    print(paste("g = ", g))
    #get ith simulation
    X.t <- X.t_sims_flk[i,]
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, W = setW, K = setK)
    
    #calculate bandpass variance
    temp_bp <- integrate(approxfun(f, MTSE_full$spectrum), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 1000)
    bpvar.vec[i] <- 4*temp_bp$value
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat[g,] <- trfunc.vec
  bmat[g,] <- bpvar.vec
}


##AVAR calculation
### only need to run once #####
# 
# amat <- oamat <- matrix(NA, nrow = numberOfSimulations, ncol = 11)
# 
# for(i in 1:numberOfSimulations){
#   #get X.t from simulation matrix
#   print(i)
#   X.t <- X.t_sims_flk[i,]
# 
#   avar.calc <- getAvars(N,X.t, taus = c(2^(0:9), floor(N/3)))
#   amat[i,] <- avar.calc$avarRes$avars
#   oamat[i,] <- avar.calc$avarRes$overavars
# }

# likely need to save tmat and bmat to work with outside of the titans
saveRDS(tmat,paste("Results/tmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_FlickerNoiseNoGaps.Rds",sep=""))
saveRDS(bmat,paste("Results/bmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_FlickerNoiseNoGaps.Rds",sep=""))
# saveRDS(amat,paste("Results/amat",runDate,"_N",N,"_",numberOfSimulations,"sims_FlickerNoiseNoGaps.Rds",sep=""))
# saveRDS(oamat,paste("Results/oamat",runDate,"_N",N,"_",numberOfSimulations,"sims_FlickerNoiseNoGaps.Rds",sep=""))

}
