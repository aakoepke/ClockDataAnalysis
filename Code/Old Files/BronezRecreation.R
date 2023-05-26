######################### Bronez, 1985 ################################
## This script recreates the analysis described in Bronez, 1985     ###
##  for calculating weighting sequences for irregularly spaced data ###
#######################################################################

#### Libraries #####
library(tidyverse)
library(fields)
library(geigen)
library(RSpectra)
#####################     GPSS Method: Chart 1, p. 45     ###############################
## Steps:                                                                              ##
##  1. Get sample times t.n and process values x.t                                     ##
##  2. Select the analysis band A (can be refined later)                               ##
##  3. Calculate R.a, R.b                                                              ##
##  4. Solve eigenvalue problem: R.a * w.k = lambda.k * R.b * w.k                      ##
##  5. Select the K weighting sequences you will pick from 4.                          ##
##      -they should be close to unity, can control this by resizing A                 ##
##  6. Normalize w.k                                                                   ##
##  7. Calculate the integrated spectrum estimate: P.A = sum_{k = 1}^{k} |w.k * x.t|^2 ##
#########################################################################################

#### Recreating Bronez ideas ####
j <- complex(real = 0, imaginary = 1) #used as imaginary i in this dissertation


n <- t.n <- 1:50 #regularly spaced data
t.n <- 0.625 + 0.3625*n + 0.0125*n^2 #arithmetically spaced data
m <- 1:60
t.n <- 5*m[-c(1,5,17,18,23,27,32,39,53,56)]/6 #missing data


#sampling kernel: K(w) = sum_{k = 1}^{N} exp(-j*w*t.n)
L <-  5000
w <- seq(-pi,pi, length.out = L)
W <- matrix(w, nrow = length(t.n), ncol = L, byrow = TRUE)
K.w <- apply(exp(-j*t.n*W), MARGIN = 2, FUN = sum)
K.w.mag <- abs(K.w)
length(K.w)
plot(w,K.w.mag, ylim = c(-80,45), type = "l")
#lines(w,K.w.mag, ylim = c(-80,45), type = "l", col = "red")
lines(w,K.w.mag, ylim = c(-80,45), type = "l", col = "blue")


#fixed center analysis band A: |f| < f_w for 0 < f < 0.1

f.w <- seq(0,0.1*2*pi,length.out = 41)
dist.mat <- rdist(t.n)

#B = [-pi,pi]?
R.b <- 1/(pi*(dist.mat))*(sin(dist.mat*pi))
R.b[row(R.b) == col(R.b)] <- 1

top.eight <- matrix(NA, nrow = length(f.w), ncol = 8)

for(i in 1:length(f.w)){
R.a <- 1/(pi*(dist.mat))*(sin(dist.mat*f.w[i]))

R.a[row(R.a) == col(R.a)] <- f.w[i]/(pi)

#Solve the generalized eigenvalue problem
evs <- geigen(R.a,R.b, symmetric = TRUE)
lambdas <- sort(evs$values, decreasing = TRUE)

top.eight[i,] <- lambdas[1:8]

}

matplot(f.w/(2*pi),top.eight) 



#sidelobe energy

gamma.k <- 10*log10(1-top.eight)

matplot(f.c/(2*pi),gamma.k, ylim = c(-60,0)) 


##Figure 4/5/6 Bronez 1985
par(mfrow = c(2,1), mar = c(4,4,3,3))

matplot(f.w/(2*pi),top.eight,lty = 1, ylab = "Eigenvalues", xlab = "f_w", main = "Figure 4") 
matlines(f.w/(2*pi),top.eight,lty = 1) 
matplot(f.w/(2*pi),gamma.k, ylim = c(-60,0),ylab = "Sidelobe Energy (dB)", xlab = "f_w") 
matlines(f.w/(2*pi),gamma.k, ylim = c(-60,0), lty = 1) 


#fixed width analysis band A: |f-f_c| < 0.05 for -0.45 < f_c < 0.45

f.c <- seq(-0.45*2*pi,0.45*2*pi,length.out = 41)
f.w <- 0.05*2*pi
dist.mat <- rdist(t.n)

#B = [-pi,pi] or [-1/2,1/2]?
R.b <- 1/(pi*(dist.mat))*(sin(dist.mat*pi))
R.b[row(R.b) == col(R.b)] <- 1

top.eight <- matrix(NA, nrow = length(f.c), ncol = 8)

for(i in 1:length(f.c)){
  #R.a <- (1/(2*j*pi*dist.mat))*(exp(j*(f.c[i] + f.w)*dist.mat) - exp(j*(f.c[i] - f.w)*dist.mat))
  R.a <- (j/(2*pi*dist.mat))*(exp(j*(f.c[i] - f.w)*dist.mat) - exp(j*(f.c[i] + f.w)*dist.mat))
  R.a[row(R.a) == col(R.a)] <- f.w/(pi)
  R.a[upper.tri(R.a)] <- Conj(R.a[lower.tri(R.a)])
  
  #Solve the generalized eigenvalue problem
  evs <- geigen(R.a,R.b, symmetric = TRUE)
  lambdas <- sort(evs$values, decreasing = TRUE)
  
  top.eight[i,] <- lambdas[1:8]
  
}

matplot(f.c/(2*pi),top.eight ) 



#sidelobe energy

gamma.k <- 10*log10(1-top.eight)

matplot(f.c/(2*pi),gamma.k, ylim = c(-70,0)) 



##Figure 7/8/9 Bronez 1985
par(mfrow = c(2,1), mar = c(4,4,3,3))

matplot(f.c/(2*pi),top.eight,lty = 1, ylab = "Eigenvalues", xlab = "center frequency", main = "Figure 8") 
matlines(f.c/(2*pi),top.eight,lty = 1) 
matplot(f.c/(2*pi),gamma.k, ylim = c(-60,0),ylab = "Sidelobe Energy (dB)", xlab = "center frequency") 
matlines(f.c/(2*pi),gamma.k, ylim = c(-60,0), lty = 1) 


#####          Functions                     #####
### Here we turn the work above into functions ###
###  which output the best K weighting sequences #
##################################################

get.weights_bronez <- function(t.n = 1:50, K = 1, f.c = 0, f.w){
  N = length(t.n)
  dist.mat <- rdist(t.n)
  
  #B = [-pi,pi] for omega or [-1/2,1/2] for f
  R.b <- 1/(pi*(dist.mat))*(sin(dist.mat*pi))
  R.b[row(R.b) == col(R.b)] <- 1
  
  R.a <- exp(j*2*pi*f.c*dist.mat)*(sin(2*pi*f.w*(dist.mat))/(pi*dist.mat))
  R.a[row(R.a) == col(R.a)] <- f.w*2
  
  #Solve the generalized eigenvalue problem
  evs <- geigen(R.a,R.b, symmetric = TRUE)
  
  return(list("weights" = evs$vectors[,(N-K + 1):N]*sqrt(2*f.w), "eigenvalues" = sort(evs$values, decreasing = TRUE)[1:K]))
}


##Example 1: best weighting sequence for regularly spaced data
#(Figures 10 + 11)
#f_c = 0, f_w = 0.05
fig10 <- get.weights_bronez(f.w = 0.05)
plot(1:50,abs(fig10$weights)^2, ylim = c(0,0.01), type = "h")

#plot(abs(fft(fig10$weights)), type = "l")

#f_c = -0.3, f_w = 0.05
fig11 <- get.weights_bronez(f.w = 0.05, f.c = -0.3)
plot(1:50,abs(fig11$weights)^2, ylim = c(0,0.01), type = "h")

#plot(seq(0,1, length.out = length(fig11$weights)),abs(fft(fig11$weights)), type = "l")


##Example 2: best weighting sequence for arithmetically spaced data
#(Figures 12 + 13)
#f_c = 0, f_w = 0.05
n <- 1:50
t.n <- 0.625 + 0.3625*n + 0.0125*n^2 
fig12 <- get.weights_bronez(t.n = t.n, f.w = 0.05)
plot(n,abs(fig12$weights)^2,ylim = c(0,0.01),type = "h")

#plot(seq(0,1, length.out = length(fig12$weights)),abs(fft(fig12$weights)), type = "l")


#f_c = -0.3, f_w = 0.05
fig13 <- get.weights_bronez(f.w = 0.05, f.c = -0.3, t.n = t.n)
plot(1:50,0.1*abs(fig13$weights),ylim = c(0,0.05), type = "h")

plot(seq(0,1, length.out = length(fig13$weights)),abs(fft(fig13$weights)), type = "l")
abline(v = 0.7)


