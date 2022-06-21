######### Spectral Analysis for "Clock" Data ########

##### Libraries #####
library(tidyverse)


#We follow the steps for multitaper analysis outlined in "Clock Statistics Tutorial"

########### Steps for Spectral Analysis of Clock Data ###########
### 1. Estimate Spectrum                                      ###
###     (Functions can be found in multitaperClockAnalysis.R) ###
### 2. Fit a linear(?) regression to the estimated spectrum   ###
### 3. Look at the slope of the regression                    ###
#################################################################

#### Periodogram Example ####

#generate data
N = 512
x.t <- arima.sim(n = N, list(order = c(0,1,0))) #random walk
x.t <- rnorm(N) #white noise


plot(1:N, x.t, type = "l")

#periodogram spectral estimate
S.x.hat <- spec.pgram(x.t, plot = FALSE)

plot(log10(S.x.hat$freq),log10(S.x.hat$spec))

#linear regresssion

ln.reg <- lm(log10(S.x.hat$spec) ~ log10(S.x.hat$freq))
ln.reg
summary(ln.reg)

abline(a = ln.reg$coefficients[1], b = ln.reg$coefficients[2])

plot(ln.reg)



