# source("/home/aak3/NIST/ClockDataAnalysis/Code/Paper1/WhiteNoise_noGaps.R")
# source("/home/cmb15/ClockDataAnalysis/Code/Paper1/WhiteNoise_noGaps.R")
#source("C:/Users/cmb15/OneDrive - UCB-O365/NIST/ClockDataAnalysis/Code/SA_ImportantFunctions.R")
source("/home/cmb15/ClockDataAnalysis/Code/SA_ImportantFunctions.R")
library(lomb)
#test
##############################################
##############################################
### read in the file with functions

# setwd("/home/aak3/NIST/ClockDataAnalysis/Code/Paper1/")
# setwd("/home/cmb15/ClockDataAnalysis/Code/Paper1/")

#source("../SA_ImportantFunctions.R")
##############################################
##############################################

numberOfSimulations = 500
N.long = 2048 + 1000
t.n_missing <- 1:N.long
t.n_missing[c(100:500,1300:1400, 1700:1874, 2400:2722)] <- NA
N <- length(na.omit(t.n_missing))

## keeping track of how long this all takes
startTime=Sys.time()

## saving the date to label file outputs
runDate=format(Sys.Date(),"%m%d%y")

### add in determination of W and K for this data pattern?
###run1
# setWnum = 12
# setW = setWnum/2048
# setK = 3

###run2
# setWnum = 12
# setW = setWnum/2048
# setK = 4

###run3
# setWnum = 12
# setW = setWnum/2048
# setK = 5

###run4
# setWnum = 4
# setW = setWnum/2048
# setK = 3

###run5
# setWnum = 4
# setW = setWnum/2048
# setK = 4

###run6
setWnum = 5
setW = setWnum/N
setK = 8

print(setWnum)
print(setK)
######################################
###### Study 2: White Noise ##########
######   WN(0,1), gaps    #########
######################################

trfunc.vec <- bpvar.vec <- trfunc.vec_LS <- rep(NA, times = numberOfSimulations)
taus <- c(2^(0:9), floor(N/3))
tmat <- bmat <- lmat <- matrix(NA, ncol = numberOfSimulations, nrow = length(taus))

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]

##calculate tapers
V.mat <- get_tapers(t.n_missing, W = setW, K = setK)

r = 0
for(k in taus){
  r = r + 1
  tau = k
  print(paste("r = ", r))
  
  for(i in 1:numberOfSimulations){
    print(i)
    set.seed(i)
    #generate X.t
    X.t_missing <- rnorm(N.long,mean = 0, sd = 1)
    X.t_missing[c(100:500, 1300:1400, 1700:1874, 2400:2722)] <- NA
    #calculate S.hat
    MTSE_full <- MT_spectralEstimate(X.t_missing, V.mat$tapers)
    lsperio <- lsp(x = X.t_missing, times = t.n_missing, plot = FALSE)
    
    #calculate bandpass variance
    temp_bp <- integrate(approxfun(f, MTSE_full$spectrum), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 1000)
    bpvar.vec[i] <- 4*temp_bp$value
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 0
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    trfunc.vec_LS[i] <- f[2]*sum(G.vec[-1]*lsperio$power)
    
  }
  tmat[r,] <- trfunc.vec
  bmat[r,] <- bpvar.vec
  lmat[r,] <- trfunc.vec_LS
}


###also Calculate AVAR###

amat <- oamat <- matrix(NA, nrow = numberOfSimulations, ncol = length(taus))

for(i in 1:numberOfSimulations){
  set.seed(i)
  print(i)
  #generate X.t
  X.t_missing <- rnorm(N.long,mean = 0, sd = 1)
  X.t_missing[c(100:500, 1300:1400, 1700:1874, 2400:2722)] <- NA
  
  avar.calc <- getAvars(N,na.omit(X.t), taus = taus)
  amat[i,] <- avar.calc$avarRes$avars
  oamat[i,] <- avar.calc$avarRes$overavars
}


### print time to run this first part
print(startTime-Sys.time())


# likely need to save tmat and bmat to work with outside of the titans
saveRDS(tmat,paste("/home/cmb15/ClockDataAnalysis/Code/Paper1/Results/tmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseGaps.Rds",sep=""))
saveRDS(bmat,paste("/home/cmb15/ClockDataAnalysis/Code/Paper1/Results/bmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseGaps.Rds",sep=""))
saveRDS(lmat,paste("/home/cmb15/ClockDataAnalysis/Code/Paper1/Results/lmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseGaps.Rds",sep=""))
saveRDS(amat,paste("/home/cmb15/ClockDataAnalysis/Code/Paper1/Results/amat",runDate,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseGaps.Rds",sep=""))
saveRDS(oamat,paste("/home/cmb15/ClockDataAnalysis/Code/Paper1/Results/oamat",runDate,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseGaps.Rds",sep=""))



