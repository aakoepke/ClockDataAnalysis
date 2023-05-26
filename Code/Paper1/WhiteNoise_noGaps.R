# source("/home/aak3/NIST/ClockDataAnalysis/Code/Paper1/WhiteNoise_noGaps.R")
# source("/home/cmb15/ClockDataAnalysis/Code/Paper1/WhiteNoise_noGaps.R")
#test
##############################################
##############################################
### read in the file with functions

# setwd("/home/aak3/NIST/ClockDataAnalysis/Code/Paper1/")
# setwd("/home/cmb15/ClockDataAnalysis/Code/Paper1/")

source("../SA_ImportantFunctions.R")
##############################################
##############################################

numberOfSimulations = 300

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
setWnum = 4
setW = setWnum/2048
setK = 6

print(setWnum)
print(setK)
######################################
###### Study 1: White Noise ##########
######   WN(0,1), no gaps    #########
######################################
# N <- 7200
N = 2048 
trfunc.vec <- bpvar.vec <- rep(NA, times = numberOfSimulations)

tmat <- bmat <- matrix(NA, ncol = numberOfSimulations, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]


r = 0
for(k in c(2^(0:9), floor(N/3))){
  r = r + 1
  tau = k
  print(paste("r = ", r))
  
  for(i in 1:numberOfSimulations){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- rnorm(N,mean = 0, sd = 1)
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, W = setW, K = setK)
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau))) # some error here for some N choices, not sure what the issue is yet, but they return Inf and -Inf
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(MTSE_full$spectrum[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + MTSE_full$spectrum[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }
      else{
        bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat[r,] <- trfunc.vec
  bmat[r,] <- bpvar.vec
}


###also Calculate AVAR###

amat <- oamat <- matrix(NA, nrow = numberOfSimulations, ncol = 11)

for(i in 1:numberOfSimulations){
  set.seed(i)
  print(i)
  #generate X.t
  X.t <- rnorm(N,mean = 0, sd = 1)
  
  avar.calc <- getAvars(N,X.t, taus = c(2^(0:9), floor(N/3)))
  amat[i,] <- avar.calc$avarRes$avars
  oamat[i,] <- avar.calc$avarRes$overavars
}


### print time to run this first part
print(startTime-Sys.time())


# likely need to save tmat and bmat to work with outside of the titans
saveRDS(tmat,paste("Results/tmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseNoGaps.Rds",sep=""))
saveRDS(bmat,paste("Results/bmat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseNoGaps.Rds",sep=""))



