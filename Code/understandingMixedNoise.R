rm(list=ls())
source("/home/aak3/NIST/ClockDataAnalysis/Code/SA_ImportantFunctions.R")
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
numberOfSimulations = 10
N.long = 20000
t.n_missing <- 1:N.long
sort(sample(1:20000,size = 8))
missing.window <- c(453:680,5000:7700,9771:13147, 16539:18231)

# c(329:3611, 5812:6278, 7905:8684)
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
  # set.seed(i)
  
  X.t_missing <- rnorm(N.long) + 0.001*TK95(N = N.long, alpha = 1) 
  X.t_missing[missing.window] <- NA
  
  #calculate S.hat, save to matrix row
  MTSE_mat[i,] <- MT_spectralEstimate(X.t_missing, taperMatrix)$spectrum
  #periodogram (because there are no gaps)
  #perio_mat[i,2:(N/2+1)] <- spec.pgram(combo, plot = FALSE)$spec
  
}

tdf=data.frame()
for(k in taus){
  print(k)
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
  
  tdftemp <- data.frame(tau=tau,avar=trfunc.vec,Method="Spectral")
  tdf=bind_rows(tdf,tdftemp)
  
}

##AVAR calculation
### only need to run once #####
# 
amat <- oamat <- data.frame()#matrix(NA, nrow = numberOfSimulations, ncol = length(taus))

for(i in 1:numberOfSimulations){
  print(i)
  # set.seed(i)
  
  X.t_missing <- rnorm(N.long) + 0.001*TK95(N = N.long, alpha = 1) 
  X.t_missing[missing.window] <- NA
  
  avar.calc <- getAvars(N.long,na.omit(X.t_missing), taus = taus)
  # amat[i,] <- avar.calc$avarRes$avars
  oamattemp <- data.frame(tau=taus,avar=avar.calc$avarRes$overavars,Method="Current")
  oamat=bind_rows(oamat,oamattemp)
}


# tmat
# oldavars=data.frame(tau=taus,avar=oamat[1,],Method="Current")
# myavars=data.frame(tau=taus,avar=tmat[1,],Method="Spectral")

# allavars=bind_rows(oldavars,myavars)
allavars=bind_rows(tdf,oamat)
truth=data.frame(tau=taus,
                 avar=WNFN_theoreticalAVAR(taus = taus,
                                           e_white = 1,
                                           flicker_floor = .001),
                 Method="Truth")

ggplot(allavars,aes(tau,avar,col=Method))+
  geom_point(size=2)+
  ### add true straight line below
  geom_line(data = truth)+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))

ggplot(allavars,aes(tau,avar,col=Method,group=interaction(tau,Method)))+
  geom_boxplot()+
  ### add true straight line below
  geom_line(data = truth,mapping = aes(tau,avar,group=1))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))


Sys.time()-startTime
