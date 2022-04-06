# some code here: http://faculty.washington.edu/dbp/sauts/R-code/chapter-08.R

library(readr)
library(ggplot2)

primerData1 <- read_csv("Data/primerData1.txt", col_names = FALSE, skip = 10)
colnames(primerData1)="X1tilde"
primerData1$t = 1:dim(primerData1)[1]

ggplot(primerData1,aes(t,X1tilde))+
  geom_point()

test=spectrum(primerData1$X1tilde)

res=data.frame(f=test$freq,spec=test$spec)

ggplot(res,aes(f,spec))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()


arSpec=spectrum(primerData1$X1tilde,method = "ar")
arSpec

arRes=data.frame(f=arSpec$freq,spec=arSpec$spec)

ggplot(arRes,aes(f,spec))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()



################ second dataset

primerData2 <- read_csv("Data/primerData2.txt", col_names = FALSE, skip = 10)
colnames(primerData2)="X2tilde"
primerData2$t = 1:dim(primerData2)[1]

ggplot(primerData2,aes(t,X2tilde))+
  geom_point()

test=spectrum(primerData2$X2tilde)

res=data.frame(f=test$freq,spec=test$spec)

ggplot(res,aes(f,spec))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()


arSpec=spectrum(primerData2$X2tilde,method = "ar")
arSpec

arRes=data.frame(f=arSpec$freq,spec=arSpec$spec)

ggplot(arRes,aes(f,spec))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()




######## Periodogram ###########
#center the time series
X.tilde <- X - mean(X)
length(X.tilde)
plot(X.tilde, type = "l")
abline(h = 0)

##zero pad X.tilde
X.tilde.prime <- c(X.tilde,rep(0,times = 96))
length(X.tilde.prime)

N <-length(X)
N.prime <- length(X.tilde.prime)
freqs <- 2*pi*(0:(N.prime-1))/N.prime #f_j's

I <- cbind(abs(fft(X.tilde.prime))^2/N.prime) #periodogram

par(mfrow=c(1,1))
plot(I~freqs, ylab = "",xlab = "", type = "l")

plot(log10(I[1:(N.prime/2+1)])~log10(freqs[1:(N.prime/2+1)]), ylab = "",xlab = "", type = "l")



###### Multitaper Analysis ########
f.nyquist <- 1/2
f.j.prime <- seq(0,f.nyquist, length.out = N.prime/2 + 1)

#function to create taper h'_k,N,t
h.k.prime <- function(k,N,N.prime){
  t <- 0:(N-1)
  h.k.N <- (2/(N+1))^(1/2)*sin((k + 1)*pi*(t + 1)/(N + 1)) #value of h' taper for t = 0, .., N-1
  if(N.prime > N){
    return(c(h.k.N, rep(0, times = N.prime-N)))#add in the 0's to pad the end to get up to 4096
  }
  else{
    return(h.k.N)
  }
}

plot(h.k.prime(k = 0, N = 4000, N.prime = 4096), type = "l") #look at taper for k = 0 with zero padding
plot(h.k.prime(k = 1, N = 4000, N.prime = 4096), type = "l") #look at taper for k = 1 with zero padding

#Recreating Figure 4
par(mfrow = c(1,1), mar = c(3,3,3,3))
for(k in 0:5){
  h.k <- h.k.prime(k = k, N = 4000, N.prime = 4000)
  plot(h.k, type = "l", ylab = paste("k = ",k))
  plot(h.k*X.tilde, type = "l", ylab = "X.tilde*h.k")
}


#building multitaper S_x(f) estimator
t <- 0:(N.prime-1)
S.x.hat <- rep(NA, times = N.prime/2 + 1) #where S_x.hat(f'_j) will be saved

for(j in 0:N.prime/2){
  k.vec <- rep(NA,times = 6)
  for(k in 0:5){
    W.t <- h.k.prime(k = k, N=N, N.prime = N.prime)*X.tilde.prime
    inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*t*j/N.prime))
    k.vec[k + 1] <- abs(inner.sum)^2
  }
  S.x.hat[j] <- (1/length(k.vec))*sum(k.vec)
}

plot(log10(f.j.prime), log10(S.x.hat), type = "l") #plot of multitaper spectral estimator


