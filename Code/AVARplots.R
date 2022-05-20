##### AVAR Plots ######
library(tidyverse)
library(forecast)
library(avar)


###Gaussian White Noise####

N = 2048
set.seed(10)
y=rnorm(N)

AVAR.WN <- getAvars(N,y)

plot(log10(AVAR.WN$avarRes$taus), log10(AVAR.WN$avarRes$avars), xlim = c(0,5),ylim = c(-5,1) )

lin.fit <- lm(log10(AVAR.WN$avarRes$avars[1:1000]) ~ log10(AVAR.WN$avarRes$taus[1:1000]))
abline(a = lin.fit$coefficients[1], b = lin.fit$coefficients[2])
lin.fit
length(log10(AVAR.WN$avarRes$taus))
plot(lin.fit)

#theoretical allan variance
avar_WN <- function(N,sigma.2){
  taus <- seq(1,N/2, by = 1)
  return(sigma.2/taus)
}


##look at a bunch of simulations

avar_saved_WN <- matrix(NA, nrow = 100, ncol = 112)

for(i in 1:100){
  set.seed(i)
  y <- rnorm(N)
  avar_saved_WN[i,] <- getAvars(N,y)$avarRes$avars
}

colors = rainbow(100)
plot(log10(seq(1,N/2, by = 1)), log10(avar_WN(N,1)), type = "l", lwd = 2, ylim = c(-5,0), ylab = "log10(AVAR)", xlab = "log10(tau)", main = "Theoretical AVAR vs. 100 Simulations for WN(0,1) Process")
for(i in 1:100){
  points(log10(AVAR.WN$avarRes$taus), log10(avar_saved_WN[i,]),col = "red")
}
lines(log10(seq(1,N/2, by = 1)), log10(avar_WN(N,1)), type = "l", lwd = 2, ylab = "AVAR", xlab = "tau")

#Show averages for each tau
avar_WN_means <- apply(log10(avar_saved_WN), MARGIN = 2, FUN = "mean")

points(log10(AVAR.WN$avarRes$taus), avar_WN_means, pch = 19, col = "blue")

legend(x = 0, y = -1, legend = c("Theoretical", "Estimated", "Mean"), col = c("black", "red", "blue"), lty = c(1,NA,NA),pch = c(NA,1,19), lwd = c(2,1,1))



###AR(1) Noise####
set.seed(10)
phi = 0.5
y <- arima.sim(n = N, list(ar = c(phi)))

AVAR.AR1 <- getAvars(N,y)

plot(log10(AVAR.AR1$avarRes$taus), log10(AVAR.AR1$avarRes$avars), xlim = c(0,5),ylim = c(-5,1) )

#theoretical allan variance
avar_AR1 <- function(phi,N,sigma.e){
  taus <- seq(1,N/2, by = 1)
  return((taus - 3*phi - taus*phi^2 + 4*phi^(taus + 1) - phi^(2*taus + 1))/(taus^2*(1-phi)^2*(1-phi^2)) )
}


plot(log10(seq(1,N/2, by = 1)), log10(avar_AR1(phi, N, sigma.e = 1)),ylim = c(-6,0), xlim = c(0,4), type = "l")
points(log10(AVAR.AR1$avarRes$taus), log10(AVAR.AR1$avarRes$avars), col = "red")
points(log10(AVAR.AR1_2$avarRes$taus), log10(AVAR.AR1_2$avarRes$avars), col = "blue")
lines(log(seq(1,N/2, by = 1)), log(av_ar1(seq(1,N/2, by = 1),phi, 1)), col = "green")


## look at a bunch of simulations

avar_saved_AR1 <- matrix(NA, nrow = 100, ncol = 112)

for(i in 1:100){
  set.seed(i)
  y <- arima.sim(n = N, list(ar = c(0.5)))
  
  avar_saved_AR1[i,] <- getAvars(N,y)$avarRes$avars
}

colors = rainbow(100)
plot(log10(seq(1,N/2, by = 1)), log10(av_ar1(seq(1,N/2,by = 1),phi,1)), type = "l", lwd = 2, ylim = c(-5,0), ylab = "log10(AVAR)", xlab = "log10(tau)", main = "Theoretical AVAR vs. 100 Simulations for AR(1) Process")
for(i in 1:100){
  points(log10(AVAR.AR1$avarRes$taus), log10(avar_saved_AR1[i,]),col = "red")
}
lines(log10(seq(1,N/2, by = 1)), log10(av_ar1(seq(1,N/2,by = 1),phi,1)), type = "l", lwd = 2, ylab = "AVAR", xlab = "tau")

#Show averages for each tau
avar_AR1_means <- apply(log10(avar_saved_AR1), MARGIN = 2, FUN = "mean")

points(log10(AVAR.AR1$avarRes$taus), avar_AR1_means, pch = 19, col = "blue")

legend(x = 0, y = -1, legend = c("Theoretical", "Estimated", "Mean"), col = c("black", "red", "blue"), lty = c(1,NA,NA),pch = c(NA,1,19), lwd = c(2,1,1))



####MA(1) Noise ######
N = 2048
theta = 0.5
avar_MA1 <- function(theta,N,sigma.e){
  taus <- seq(1,N/2, by = 1)
  return((1 + theta^2)*((taus - (2*taus - 3)*(theta/(1 + theta^2)))/taus^2)*sigma.e)
}

plot(log10(seq(1,N/2, by = 1)), log10(avar_MA1(theta, N, sigma.e = 1)))

y <- arima.sim(n = N, list(ma = c(0.5)), sd = 1)

AVAR.MA1 <- getAvars(N,y)
points(log10(AVAR.MA1$avarRes$taus), log10(AVAR.MA1$avarRes$avars), col = "red")


#look at a bunch of simulations
theta = 0.5
avar_saved_MA1 <- matrix(NA, nrow = 100, ncol = 112)

for(i in 1:100){
  set.seed(i)
  y <- arima.sim(n = N, list(ma = c(theta)))
  
  avar_saved_MA1[i,] <- getAvars(N,y)$avarRes$avars
}


plot(log10(seq(1,N/2, by = 1)), log10(avar_MA1(-theta,N,1)), type = "l", lwd = 2, ylim = c(-5,1),ylab = "log10(AVAR)", xlab = "log10(tau)", main = "Theoretical AVAR vs. 100 Simulations for MA(1) Process")

for(i in 1:100){
  points(log10(AVAR.AR1$avarRes$taus), log10(avar_saved_MA1[i,]),col = "red")
}
lines(log10(seq(1,N/2, by = 1)), log10(avar_MA1(-theta,N,1)), type = "l", lwd = 2, ylab = "AVAR", xlab = "tau")

#Show averages for each tau
avar_MA1_means <- apply(log10(avar_saved_MA1), MARGIN = 2, FUN = "mean")

points(log10(AVAR.AR1$avarRes$taus), avar_MA1_means, pch = 19, col = "blue")

legend(x = 0, y = -2.5, legend = c("Theoretical", "Estimated", "Mean"), col = c("black", "red", "blue"), lty = c(1,NA,NA),pch = c(NA,1,19), lwd = c(2,1,1))



####Random Walk######

y <- arima.sim(n = N, list(order = c(0,1,0)))
AVAR.RW <- getAvars(N,y)

points(log10(AVAR.RW$avarRes$taus), log10(AVAR.RW$avarRes$avars))

plot(log10(seq(1,N/2,by = 1)), log10(av_rw(1,seq(1,N/2,by = 1))))


##look at a bunch of simulations

avar_saved_RW <- matrix(NA, nrow = 100, ncol = 112)

for(i in 1:100){
  set.seed(i)
  y <- arima.sim(n = N, list(order = c(0,1,0)))
  
  avar_saved_RW[i,] <- getAvars(N,y)$avarRes$avars
}


plot(log10(seq(1,N/2, by = 1)), log10(av_rw(1,seq(1,N/2, by = 1) )), type = "l", lwd = 2,ylab = "log10(AVAR)", xlab = "log10(tau)", main = "Theoretical AVAR vs. 100 Simulations for RW")

for(i in 1:100){
  points(log10(AVAR.AR1$avarRes$taus), log10(avar_saved_RW[i,]),col = "red")
}
lines(log10(seq(1,N/2, by = 1)), log10(av_rw(1,seq(1,N/2, by = 1))), type = "l", lwd = 2, ylab = "AVAR", xlab = "tau")

#Show averages for each tau
avar_RW_means <- apply(log10(avar_saved_RW), MARGIN = 2, FUN = "mean")

points(log10(AVAR.AR1$avarRes$taus), avar_RW_means, pch = 19, col = "blue")

legend(x = 0, y = 2, legend = c("Theoretical", "Estimated", "Mean"), col = c("black", "red", "blue"), lty = c(1,NA,NA),pch = c(NA,1,19), lwd = c(2,1,1))


###First Difference of RW####

avar_saved_RW_FD <- matrix(NA, nrow = 100, ncol = 112)

for(i in 1:100){
  set.seed(i)
  y <- arima.sim(n = N, list(order = c(0,1,0)))
  y <- diff(y)
  avar_saved_RW_FD[i,] <- getAvars(N,y)$avarRes$avars
}

plot(log10(AVAR.AR1$avarRes$taus), log10(avar_saved_RW_FD[1,]),col = "red", ylab = "log10(AVAR)", xlab = "log10(tau)", main = "First Difference of RW Experiment")

for(i in 1:100){
  points(log10(AVAR.AR1$avarRes$taus), log10(avar_saved_RW_FD[i,]),col = "red")
}
lines(log10(seq(1,N/2, by = 1)), log10(avar_WN(N,1)), type = "l", lwd = 2)

#averages
avar_RW_FD_means <- apply(log10(avar_saved_RW_FD), MARGIN = 2, FUN = "mean")

points(log10(AVAR.AR1$avarRes$taus), avar_RW_FD_means, pch = 19, col = "blue")


legend(x = 0, y = -2, legend = c("Theoretical WN", "Estimated", "Mean"), col = c("black", "red", "blue"), lty = c(1,NA,NA),pch = c(NA,1,19), lwd = c(2,1,1))


#####AVAR computation function #####

avar_fn=function(y,tau){
  n=length(y)
  
  div=seq(1,n,by = tau)
  
  M=length(div)-1 #number of groups
  
  groupmeans = numeric(M)
  for(i in 1:M){
    groupmeans[i]=mean(y[div[i]:(div[i+1]-1)])
  }
  
  1/(2*(M-1)) * sum(diff(groupmeans)^2)
}

getAvars=function(N, y){
  taus = c(1:10,15,20)
  taus = c(taus,seq(30,N/2,by=10)) 
  
  avars=numeric(length(taus))
  
  for (i in 1:length(taus)){
    avars[i]=avar_fn(y,taus[i])
  }
  
  m1=data.frame(taus=taus,avars=avars)  
  fit=lm(log(sqrt(avars))~log(taus),data = m1)
  slope=as.numeric(fit$coefficients[2])
  int=as.numeric(fit$coefficients[1])
  
  avarRes=data.frame(taus=taus,avars=avars,N=N,slope=slope,int=int)  
  
  ##########################################
  # get SE
  ##########################################
  SEests=data.frame()
  
  new=data.frame(N=N,out=exp(int+slope*log(N)), type="AD")
  new2=data.frame(N=N,out=sd(y)/sqrt(N), type="SE")
  new3=data.frame(N=N,out=1/sqrt(N), type="true")
  SEests=bind_rows(new,new2,new3)
  
  return(list(avarRes=avarRes,SEests=SEests))
}



