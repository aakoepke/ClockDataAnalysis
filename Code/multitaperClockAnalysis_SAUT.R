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



########################## MTSE function #################################
## inputs: time series, possibly with gaps, W = bandwidth, K = number of tapers, delta.t = time step
## outputs: tapers, eigenvalues, spectral leakage?, spectral estimate



multitaper_est <- function(X.t, NW, K){
  X.t <- X.t - mean(X.t, na.rm = TRUE) #demean
  N.long <- length(X.t)
  t.n <- 1:N.long
  missing.indices <- which(is.na(X.t))
  t.n[which(is.na(X.t))] <- NA


  W <- NW/length(na.omit(t.n))
  dist.mat <- rdist(na.omit(t.n))
  
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  eigdec <- eigs_sym(A.prime, k = K, which = "LM")
  
  
  eig_vecs <- eigdec$vectors #get only the vectors
  
  #some sign maintenance
  for(i in seq(1,K,by = 2)){
    if (mean(Re(eig_vecs[,i]))<0){
      eig_vecs[,i] <- -eig_vecs[,i]
    }
  }
  
  for(i in seq(2,K-1,by = 2)){
    if (Re(eig_vecs[2,i] - eig_vecs[1,i])<0){
      eig_vecs[,i] <- -eig_vecs[,i]
    }
  }
  
  
  ##use tapers to generate spectral estimate
  N <- length(na.omit(t.n))
  S.x.hat_MD <- rep(NA, times = floor(N/2) + 1)
  
  for(j in 0:floor(N/2)){
    k.vec <- rep(NA,times = K)
    for(k in 0:(K-1)){
      W.t <- eig_vecs[,k+1]*na.exclude(X.t)
      inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*na.omit(t.n)*j/N), na.rm = TRUE)
      k.vec[k + 1] <- abs(inner.sum)^2
    }
    S.x.hat_MD[j+1] <- mean(k.vec)
  }
  
  
  return(list("tapers" = eig_vecs, "e.values" = eigdec$values, "spectrum" = S.x.hat_MD))
}


plot(S.x.hat_MD)

### Test: White Noise ####
N <- 2048
set.seed(130)
X.t <- X.t_missing <-  rnorm(N)
X.t_missing[c(20:35, 600:800,900:1010)] <- NA

#EDA
plot(X.t)
hist(X.t)
acf(na.omit(X.t))

#MTSE of white noise with missing data
MTSE_wn <- multitaper_est(X.t = X.t, NW = 2, K = 5)
MTSE_wn_missing <- multitaper_est(X.t = X.t_missing, NW = 2, K = 5)

par(mfrow = c(2,1))
plot(log10(seq(0,0.5,length.out = length(MTSE_wn$spectrum)))[-1], log10(MTSE_wn$spectrum)[-1], type = "l")
lines(log10(seq(0,0.5,length.out = length(MTSE_wn_missing$spectrum)))[-1], log10(MTSE_wn_missing$spectrum[-1]), col = "red")


#fit a line to the spectrum:
f <- log10(seq(0,0.5,length.out = length(MTSE_wn$spectrum)))[-1]
f.missing <- log10(seq(0,0.5,length.out = length(MTSE_wn_missing$spectrum)))[-1]
S.f <- log10(MTSE_wn$spectrum)[-1]
S.f_missing <- log10(MTSE_wn_missing$spectrum)[-1]

lin.fit <- lm(S.f ~ f)
lin.fit.missing <- lm(S.f_missing ~ f.missing)

plot(lin.fit)
plot(lin.fit.missing)

summary(lin.fit)
summary(lin.fit.missing)

plot(f, S.f, type = "l")
lines(f.missing,S.f_missing, col = "red")
abline(a = lin.fit$coefficients[1], b = lin.fit$coefficients[2])
abline(a = lin.fit.missing$coefficients[1], b = lin.fit.missing$coefficients[2], col = "red")
abline(h = 0)
##look at eig_vecs
MTSE_wn$e.values
MTSE_wn_missing$e.values

###### Test 2: ARFIMA(0,d,0) #######
N <- 2^14
set.seed(100)
d <- 0.4
X.t <- X.t_missing <- arfima.sim(N,model = list(dfrac = d))
X.t_missing[c(20:35, 600:800,900:1010)] <- NA

#MTSE of arfima with and without missing data
MTSE_fd <- multitaper_est(X.t = X.t, NW = 2, K = 5)
MTSE_fd_missing <- multitaper_est(X.t = X.t_missing, NW = 2, K = 5)

f <- log10(seq(0,0.5,length.out = length(MTSE_fd$spectrum)))[-1]
f.missing <- log10(seq(0,0.5,length.out = length(MTSE_fd_missing$spectrum)))[-1]

plot(f,log10(MTSE_fd$spectrum[-1]), type = "l")
lines(f.missing, log10(MTSE_fd_missing$spectrum[-1]), col = "red")
lines(log10(seq(0.001,0.5,length.out = 1025)),-2*d*log10(2*abs(sin(abs(pi*seq(0.005,0.5,length.out = 1025))))), col = "blue")


lines(log10(f.j), -2*d*log10(2*sin(pi*f.j)), col = "blue")


MTSE_fd$e.values
MTSE_fd_missing$e.values


