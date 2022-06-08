##### AVAR Plots ######
library(tidyverse)
library(forecast)
library(avar)


###Gaussian White Noise####

N = 4096
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

#theoretical confidence bounds:
avar_CI <- function(y,CI.level,noise_type = "white noise", avar_type, len){
  M <- len
  N <- length(y) + 1
  s.2 <- rep(NA, times = M)
  a <- (1-CI.level)/2
  
  if(avar_type == "overlapping"){
  for(k in 1:(M)){
    s.2[k] <- overlapping_avar_fn(y,m = k)
  }
  }
  if(avar_type == "regular"){
    for(k in 1:(M)){
    s.2[k] <- avar_fn(y,tau = k)
    }
  }
  
  if(avar_type == "theoretical W FM"){
      s.2 <- avar_WN(N,1)[1:M]
  }
  
  edf <- rep(NA, times = M)
  
  if(noise_type == "white noise"){
    for(i in 1:(M)){
      edf[i] <- ((3*(N-1)/(2*i)) - (2*(N-2)/N))*(4*i^2)/(4*i^2 + 5)
    }
  }
  
  if(avar_type == "regular"){
    CI.limits <- bind_rows("lower" = s.2 - s.2/N, "upper" = s.2 +  s.2/N)
  }
  else{
  CI.limits <- bind_rows("lower" = s.2*edf/qchisq(1-a,edf),"upper" = s.2*edf/qchisq(a, edf) )
  }
  return(CI.limits)
}

N = 4096
y <- rnorm(N)
#use the means of overlapping avar for WN to create bounds
y <- oavar_WN_means
test <- avar_CI(y,CI.level = 0.90, noise_type = "white noise", avar_type = "theoretical W FM", len = 300)
test2 <- avar_CI(y,CI.level = 0.90, noise_type = "white noise", avar_type = "overlapping", len = 300)
test3 <- avar_CI(y,CI.level = 0.90, noise_type = "white noise", avar_type = "regular", len = 300)


plot(log10(seq(1,300, by = 1)), log10(avar_WN(N,1)[1:300]), type = "l", lwd = 2, ylab = "log10(AVAR)", xlab = "log10(tau)", main = "Theoretical AVAR vs. 100 Simulations for WN(0,1) Process")
lines(log10(1:300), log10(test$lower), lty = 2)
lines(log10(1:300), log10(test$upper), lty = 2)
lines(log10(1:300), log10(test2$lower), lty = 2, col = "red") #overlapping 
lines(log10(1:300), log10(test2$upper), lty = 2, col = "red") #overlapping
#lines(log10(1:300), log10(test3$lower), lty = 2, col = "blue") #regular?
#lines(log10(1:300), log10(test3$upper), lty = 2, col = "blue") #regular?

##look at a bunch of simulations
avar_saved_WN <- matrix(NA, nrow = 100, ncol = 214)
ms <- c(1:10,15,20,seq(25,N/2, by = 50))
oavar_saved_WN <- matrix(NA, nrow = 100, ncol = length(ms))
modavar_saved_WNPM <- matrix(NA, nrow = 100, ncol = 80)
oavar.WN <- rep(NA, times = length(ms))

for(i in 1:100){
  set.seed(i)
  y <- rnorm(N)
  #x <- diffinv(y)
  
  for(k in 1:length(ms)){
    print(k)
    oavar.WN[k] <- overlapping_avar_fn(y,m = ms[k])
  }
  oavar_saved_WN[i,] <-  oavar.WN   # original avar: getAvars(N,y)$avarRes$avars
  #avar_saved_WN[i,] <- getAvars(N,y)$avarRes$avars
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

modavar_WNPM_means <- apply(log10(modavar_saved_WNPM), MARGIN = 2, FUN = "mean")
points(log10(ms),modavar_WNPM_means, pch = 19, col = "dark green")

oavar_WN_means <- apply(log10(oavar_saved_WN), MARGIN = 2, FUN = "mean")
points(log10(ms),oavar_WN_means, pch = 19, col = "dark green")

legend(x = 2.5, y =0, legend = c("Theoretical", "AVAR", "OAVAR"), col = c("black", "blue", "dark green"), lty = c(1,NA,NA),pch = c(NA,19,19), lwd = c(2,1,1))




#difference in mean to theoretical
diffs <- log10(avar_WN(N,1)[AVAR.WN$avarRes$taus])- avar_WN_means
plot(log10(AVAR.WN$avarRes$taus), diffs, ylab = "log(Theoretical) - log(Means)", xlab = "log(tau)", pch = 19)
 

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
taus = c(1:10,15,20,seq(30,N/2,by=10)) #taus from getAvars function
avar_saved_AR1 <- matrix(NA, nrow = 100, ncol = 214)
oavar_saved_AR1 <- matrix(NA, nrow = 100, ncol = 214)
oavar.AR1 <- rep(NA, times = length(taus))

for(i in 1:100){
  set.seed(i)
  print(i)
  y <- arima.sim(n = N, list(ar = c(0.5)))
  
  avar_saved_AR1[i,] <- getAvars(N,y)$avarRes$avars
  
  for(k in 1:length(taus)){
    oavar.AR1[k] <- overlapping_avar_fn(y, m = taus[k])
  }
  oavar_saved_AR1[i,] <- oavar.AR1
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




#######ARFIMA(0,d,0) d < 0.5 ###########
##Formulae taken from Zhang, 2008


#use R package function for the theoretical acvf of ARFIMA
library(arfima)

tavar_ARFIMA <- function(N,d, sig.2.a){
  rho.vec <- tacvfARFIMA(phi = 0, theta = 0, dfrac = d, maxlag = 2*N)
  corr.vec <- rho.vec/max(rho.vec) #normalize
  taus <- 2:N

  sum.vec <- rep(NA, times = N-1)
  for(k in 2:N){
    print(k)
    total <- 0
    for(i in 1:(k - 1)){
      print(i) 
      total = total + i*(2*corr.vec[k-i + 1] - corr.vec[i + 1] - corr.vec[2*k-i + 1])
      }
      sum.vec[k-1] <- total
    }
  numerator <-(taus*(rep(corr.vec[1],times = N-1) - corr.vec[3:(N+1)]) + sum.vec)*gamma(1-2*d)
  denom <- (taus*gamma(1-d))^2
  return(numerator/denom)
}

t1 <- tavar_ARFIMA(N = N/2, d = 0.49, sig.2.a = 1)
t2 <- tavar_ARFIMA(N = N/2, d = 0.25, sig.2.a = 1)

plot(t1,type = "l", ylim = c(0,2))
lines(t2,lty = 2, col = "red")

#checking out to what it should look like for ARFIMA(0,0.49,0) and ARFIMA(0,0.25,0), good

####Multiple Simulation Experiment ARFIMA(0,0.25,0), ARFIMA(0,0.49,0)

d <- 0.49
N = 4096
avar_saved_ARFIMA_49 <- matrix(NA, nrow = 100, ncol = 214)

for(i in 1:100){
  set.seed(i)
  y <- arfima.sim(N,model = list(dfrac = d))
  
  avar_saved_ARFIMA_49[i,] <- getAvars(N,y)$avarRes$avars
}


plot(log10(seq(2,N/2-1, by = 1)), log10(t1[2:length(t1)]),ylim = c(-2,0.2), type = "l", lwd = 2,ylab = "log10(AVAR)", xlab = "log10(tau)", main = "Theoretical AVAR vs. 100 Simulations for \n ARFIMA(0,0.49,0)")

for(i in 1:100){
  points(log10(AVAR.AR1$avarRes$taus), log10(avar_saved_ARFIMA_49[i,]),col = "red")
}
lines(log10(seq(2,N/2-1, by = 1)), log10(t1[2:length(t1)]),ylim = c(-6,1), type = "l", lwd = 2,ylab = "log10(AVAR)", xlab = "log10(tau)", main = "Theoretical AVAR vs. 100 Simulations for ARFIMA(0,0.25,0)")

#Show averages for each tau
avar_ARFIMA_means_49 <- apply(log10(avar_saved_ARFIMA_49), MARGIN = 2, FUN = "mean")

points(log10(AVAR.AR1$avarRes$taus), avar_ARFIMA_means_49, pch = 19, col = "blue")

legend(x = 0, y = -2, legend = c("Theoretical", "Estimated", "Mean"), col = c("black", "red", "blue"), lty = c(1,NA,NA),pch = c(NA,1,19), lwd = c(2,1,1))






#####AVAR computation functions #####

#### Regular Allan Variance #####
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


###Overlapping AVAR ####

########overlapping

y=1:10*1#rnorm(10)
m=3
overlapping_avar_fn=function(y,m){
  M=length(y)
  
  numberOfGroups=M-m
  groupmeans = numeric(numberOfGroups)
  
  for(i in 1:(M-m)){
    print(i:(i+m))
    groupmeans[i]=mean(y[i:(i+m)])
  }

  
  out=0
  for(i in 1:(M-2*m)){
    out=out+(groupmeans[i+m]-groupmeans[i])^2
  }
  out=out/(2*(M-2*m+1))
  
  return(out)
}


##### From: Frequency Stability Analysis Using R ######

#####MODIFIED ALLAN VARIANCE #####

# Function to calculate Modified Allan deviation
# MVAR for phase data
# Argument tau is basic data sampling interval
# Each analysis tau is tau*m
# where argument m is averaging factor 1 to N/3
pmavar<-function(x, tau=1, m=1){
  N=length(x)
  mvar=0
  # Outer loop
  for(j in 1:(N-3*m+1))
  {
    s=0
    # Inner loop
    for(i in j:(j+m-1))
    {
      s=s+(x[i+(2*m)]-2*x[i+m]+x[i])
    }
    mvar=mvar+s^2
  }
  # Scaling
  mvar=mvar/(2*m^4*tau^2*(N-3*m+1))
  return (mvar)
}


##Example of modified avar: white noise###

#frac. freq. deviates
N = 4096
y = rnorm(N)
AVAR.WN <- getAvars(N,y)

#phase data

x <- diffinv(y)

#calculate pmdev
max.tau <- floor(N/3)
pmdev.WN <- rep(NA, times = max.tau)
ms <- c(1:10,15,20,seq(25,floor(N/3), by = 20))
for(m in ms){
pmdev.WN[m] <- pmdev(x,tau = 1, m = m)
}


plot(log10(ms), na.omit(log10(pmdev.WN^2)), ylim = c(-4,0.5), ylab = "log(avar)", xlab = "log(tau)")

points(log10(AVAR.WN$avarRes$taus),log10(AVAR.WN$avarRes$avars), col = "red")


##Example of modified avar: AR(1) process###

#frac. freq. deviates
N = 4096
y = arima.sim(n = N, list(ar = 0.5))
AVAR.AR1 <- getAvars(N,y)

#phase data

x <- diffinv(y)

#calculate pmdev
pmdev.AR1 <- rep(NA, times = N/2)
for(m in 1:(N/2)){
  pmdev.AR1[m] <- pmdev(x,m)
}


plot(log10(1:(N/2)), log10(pmdev.AR1))

points(log10(AVAR.WN$avarRes$taus),log10(AVAR.WN$avarRes$avars), col = "red")



### fractional frequency deviates ####

y = rnorm(2048)

##this takes forever using fractional frequency data
modified_avar_fun <- function(y, m = 1){
  M = length(y)
  mvar=0
  # Outer loop
  for(j in 1:(M-3*m+2)){
    
    #print(paste("j = ", j))
    s = 0
    # First Inner loop
    for(i in j:(j+m-1)){
      
      #print(paste("i = ", i))
      #Second Inner loop
      for(k in i:(i + m - 1)){
        
        #print(paste("k = ", k))
        s = s + y[k + m] - y[k]
      }
      
    }
    mvar=mvar+s^2
  }
  
  # Scaling
  mvar=mvar/(2*m^4*(M-3*m+2))
  return (mvar)
}

test.modifiedAVAR <- rep(NA, times = 1000)

for(m in c(1:10,15,20,seq(25,500, by = 20))){
  test.modifiedAVAR[m] <-  modified_avar_fun(y,m = m)
    }

plot(log10(1:1000), log10(test.modifiedAVAR))
points(log10(getAvars(1000, y)$avarRes$avars), col = "red")





library(dplyr)

data_list = list()

for(i in 1:100){data_list[[i]] = data.frame("x" = rnorm(1), "y" = runif(1))}

df = bind_rows(data_list)



