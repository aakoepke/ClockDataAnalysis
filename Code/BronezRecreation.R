######################### Bronez, 1985 ################################
## This script recreates the analysis described in Bronez, 1985     ###
##  for calculating weighting sequences for irregularly spaced data ###
#######################################################################


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

#### Example 1: Regularly Spaced Data  ####

t.n <- 1:50

j <- complex(real = 0, imaginary = 1)

#sampling kernel: K(w) = sum_{k = 1}^{N} exp(-j*w*t.n)
L <-  40
w <- seq(-1,1, length.out = L)
W <- matrix(w, nrow = length(t.n), ncol = L, byrow = TRUE)
K.w <- apply(exp(-j*t.n*W), MARGIN = 2, FUN = sum)
K.w.mag <- abs(K.w)
length(K.w)
plot(w,K.w.mag)


#fixed center analysis band A: |f| < f_w for 0 < f < 0.1
library(fields)
library(geigen)

f.w <- seq(0,0.1*2*pi,length.out = 41)
dist.mat <- rdist(t.n)

#B = [-pi,pi]?
R.b <- 1/(2*pi*j*(dist.mat))*(exp(j*pi*dist.mat) - exp(-j*pi*dist.mat))
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


#what is going on with the f.w's? Why can't I get the same scale as in the paper?
#maybe it has something to do with 2*pi * f = omega


#sidelobe energy

gamma.k <- 10*log10(1-top.eight)

matplot(f.w/(2*pi),gamma.k, ylim = c(-60,0)) 


##Figure 4 Bronez 1985
par(mfrow = c(2,1), mar = c(4,4,3,3))

matplot(f.w/(2*pi),top.eight,lty = 1, ylab = "Eigenvalues", xlab = "f_w", main = "Figure 4") 
matlines(f.w/(2*pi),top.eight,lty = 1) 
matplot(f.w/(2*pi),gamma.k, ylim = c(-60,0),ylab = "Sidelobe Energy (dB)", xlab = "f_w") 
matlines(f.w/(2*pi),gamma.k, ylim = c(-60,0), lty = 1) 





