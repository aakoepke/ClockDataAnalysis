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

plot(log10(S.x.hat$freq),log10(S.x.hat$spec), type = "l")

#linear regresssion

ln.reg <- lm(log10(S.x.hat$spec) ~ log10(S.x.hat$freq))
ln.reg
summary(ln.reg)

abline(a = ln.reg$coefficients[1], b = ln.reg$coefficients[2])

plot(ln.reg)


#Random Walk (ARFIMA(0,1,0))
f <- seq(-0.5,0.5, by = 0.01)
plot(f,abs(1-exp(-j*2*pi*f))^2, type = "l")
lines(f, abs(2*sin(pi*f))^2, col = "blue")

spec.RW <- function(f,sigma.e = 1){
  return(sigma.e/abs(2*sin(pi*f))^2)
}

lines(log10(f),log10(spec.RW(f)), col = "blue")
par(mfrow = c(2,1))


#White Noise
spec.WN <- function(sigma.e = 1,N){
  return(rep(sigma.e,times = N))
}
lines(log10(f), log10(spec.WN(1,N = length(f))), col = "blue")



