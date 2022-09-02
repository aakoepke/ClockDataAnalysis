########## Method Comparison ###############
## This script compares the traditional   ##
## AVAR methods of analysis with spectral ##
## approach and includes some EDA tools   ##
############################################
library(purrr) #discrete uniform dist.
library(tidyverse)
library(RSpectra)


####       Simulated Data     #####
## A. white noise                ##
## B. flicker noise?             ##
## C. ARFIMA/combo of ARFIMA     ##
###################################

#create data
N <-  2048
X.t <- X.t_missing <-  rnorm(N, mean = 0, sd = 0.05)
#create data with gaps

###random gaps#####
###set.seed(20)
###endpoints <- sort(rdunif(4,2048,1))

###placed gaps####
endpoints <- c(450,700,1000,1200)


X.t_missing[c(endpoints[1]:endpoints[2], endpoints[3]:endpoints[4])] <- NA


#look at the data
par(mfrow = c(2,1))
plot(X.t, type = "l")
plot(X.t_missing, type = "l")


#### Quick little look at X.t and X.t_missing ####

#histogram
hist(X.t)
hist(X.t_missing)

#acfs
acf(X.t)
acf(na.omit(X.t_missing))




#############                    Two Sample Variance Method                        ####################
### 1.	Concatenate to remove gaps                                                                  ###
### 2.	Calculate Allan Deviation (overlapping)                                                     ###
### 3.	For clock ratios containing Al start cut off to use ADEV                                    ###
###      points above 10-20 seconds and for Yb or Sr maybe start with 20 seconds to be safe         ###
### 4.	Fit a line for ADEV points for tau = 10-20 seconds (depending on step 3 above)              ###
###      through tau = 1/3 the length of the dataset where you set the slope of the line            ###
###      to follow 1/sqrt(tau) and solve for the intercept                                          ###
### 5.	Extrapolate that line out to tau = length of the dataset (or until you detect a floor perhaps)#
#######################################################################################################


#### 1. Concatenate gappy data  ####
X.t_concat <- na.omit(X.t_missing)

plot(X.t_concat, type = "l")
N.short <- length(X.t_concat)


#### 2. Calculate ADEV          ####
AVAR.WN <- getAvars(N,X.t, taus = c(1:20, 2^(5:9)))

plot(log10(AVAR.WN$avarRes$taus), log10(AVAR.WN$avarRes$avars), ylab = "log10(AVAR)", xlab = "log10(tau)",pch = 19, ylim = c(-6,-1))

#calculate AVAR of X.t with missing gaps, pushed together
AVAR.WN_missing <- getAvars(N.short,X.t_concat, taus = c(1:20, 2^(5:9)))

points(log10(AVAR.WN_missing$avarRes$taus), log10(AVAR.WN_missing$avarRes$avars), col = "red", xlab = "log10(AVAR)", ylab = "log10(tau)",pch = 19  )


#### 3. + 4. Fit a line for points between tau = 10-1/3 of data set ####
floor(N/3) #682
N.short/3 #532

AVAR.WN <- getAvars(N,X.t, taus = seq(10,floor(N/3), by = 2))

log.avar.full <- log10(AVAR.WN$avarRes$avars)
log.taus.full <- log10(AVAR.WN$avarRes$taus)

lin.fit.full <- lm(log.avar.full ~ log.taus.full)
summary(lin.fit.full)


AVAR.WN.short <- getAvars(N.short,X.t_concat, taus = seq(10,floor(N.short/3), by = 2))

log.avar.short <- log10(AVAR.WN.short$avarRes$avars)
log.taus.short <- log10(AVAR.WN.short$avarRes$taus)

lin.fit.short <- lm(log.avar.short ~ log.taus.short)
summary(lin.fit.short)

plot(lin.fit.short)

### show the fits ###

## Full Data ##

plot(log.taus.full, log.avar.full, ylab = "log10(AVAR)", xlab = "log10(tau)",pch = 19, ylim = c(-6.5,-3)) #avar of full data with no gaps

full.intercept <- as.numeric(lin.fit.full$coefficients[1]) #intercept estimate from linear fit of full data with no gaps
full.slope <- as.numeric(lin.fit.full$coefficients[2]) #slope from linear fit of full data with no gaps

abline(a = full.intercept, b = full.slope, col = "blue", lwd = 2) #linear fit

## Partial Data ##

plot(log.taus.short, log.avar.short, ylab = "log10(AVAR)", xlab = "log10(tau)",pch = 19, ylim = c(-6.5,-3)) #plot of avar for concatenated gappy data

short.intercept <- as.numeric(lin.fit.short$coefficients[1]) #intercept estimate from linear fit of concatenated gappy data
short.slope <- as.numeric(lin.fit.short$coefficients[2]) #slope estimate from "" ""

abline(a = lin.fit.short$coefficients[1], b = lin.fit.short$coefficients[2], col = "blue", lwd = 2) #linear fit line


#### 5.	Extrapolate that line out to tau = length of the dataset (or until you detect a floor perhaps) ####


AVAR.hat.full <- exp(full.intercept + full.slope*log10(N)) #estimate of AVAR by extrapolating out to full data set

sqrt(AVAR.hat.full)

AVAR.hat.short <- exp(short.intercept + short.slope*log10(N.short)) #estimate of AVAR by extrapolating out to lenght of concatenated gappy data

sqrt(AVAR.hat.short)


############# Multitaper Spectral Estimate (MTSE)  ###############
### **see steps in google drive file and add here when ready** ####

###########################################################
##### calculate MTSE for X.t without any missing data #####
###########################################################


MTSE_full <- multitaper_est(X.t, NW = 2, K = 3) #compute the multitaper spectral estimate for the full data set

f <- seq(0,0.5,length.out = length(MTSE_full$spectrum)) #grid of frequencies

plot(log10(f[-1]), log10(MTSE_full$spectrum)[-1], type = "l") #plot of f vs. S.hat(f) on log-log scale

MTSE_full$e.values #looking that the eigenvalues are close to 1



## Bandpass Variance or Transfer Function comparison to AVAR
delta.f <- 0.5/(length(MTSE_full$spectrum)-1) #delta f, distance between each frequency in the frequency grid

#bandpass variance
tau <-  floor(N/3) #1/3 the length of the dataset

bp.var <- 2*(MTSE_full$spectrum[1]*(f[2] - 1/(4*tau)) + MTSE_full$spectrum[2]*(1/(2*tau) - f[2]))


#transfer function
transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}
G.vec <- transfer.func(f, tau)
G.vec[1] <- 1

f[2]*sum(G.vec*MTSE_full$spectrum)


###########################################################
#####     calculate MTSE for X.t with missing data    #####
###########################################################

MTSE_short <- multitaper_est(X.t_missing, NW = 2, K = 3)
lines(log10(seq(0,0.5,length.out = length(MTSE_short$spectrum)))[-1], log10(MTSE_short$spectrum)[-1], type = "l", col = "red")






############# MTSE for different W and K ##############





  
  
############# Clock Data Analysis  ##################
load("Data/ClockData.RData")
plot(clock_df$time,clock_df$AlYb)

hist(na.omit(clock_df$AlYb))
acf(na.omit(clock_df$AlYb)) #definitely not white noise







