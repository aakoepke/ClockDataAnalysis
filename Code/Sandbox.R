##### Sandbox #####

#this is a script to try things out before they go in real scripts

########################## looking at distribution of spectrum at different freqeuncies ##########################

###look at distribution of spectral-based AVAR estimates
tmat <- readRDS(file = "Code/Paper1/Results/tmat050523_W4_K6_N2048_300sims_WhiteNoiseNoGaps.Rds")


##fit with ML

##we'll do it numerically
fun.to.optimize <- function(x,k){
  n <- length(x)
  log.like <- (k*n/2)*log(2) + n*log(gamma(k/2)) - (k/2 - 1)*sum(log(x)) + sum(x)/2
}

mle.chisq <- optim(par = 3, fn = fun.to.optimize, x = spec.mat[,1], lower = 1, upper = 50, method = "Brent")
mle.chisq$par

k.seq <- seq(1,30, by = 1)
x <- rchisq(300, df = 5)
plot(k.seq, fun.to.optimize(spec.mat[,1],k.seq))



############### spectrum of white noise, just using periodogram ########

spec.mat <- matrix(NA, ncol = 33, nrow = 1000)

for(i in 1:1000){
  print(i)
  set.seed(i)
  x.t <- rnorm(64) # arfima.sim(64,model = list(dfrac = 0.45))
  t.vec <- 1:64
  
  ##For Multitaper
  s <- multitaper_est(x.t, W = 4/64, K = 2)
  spec.mat[i,] <- s$spectrum
  
  
  ##For the periodogram
  #for(j in 1:(length(x.t)/2)){
  #spec.mat[i,j] <- abs((1/sqrt(length(x.t)))*sum(x.t*exp(-im*2*pi*(j/length(x.t))*t.vec)))^2
  #}
}

hist(spec.mat[,1], freq = FALSE, breaks = 50)

spec.arfima <- function(f, d){
  (4*sin(pi*f)^2)^(-d)
}


lines(seq(0,15, length.out = 200), dchisq(seq(0,15, length.out = 200), df = 2)/2)

length(s$freq)


probs <- seq(0.01,0.99, by = 0.01)
data.quantiles <- quantile(spec.mat[,4], probs = probs)

true.quantiles <- qchisq(probs, df = 4)

plot(data.quantiles,true.quantiles/4) #*(max(data.quantiles)/max(true.quantiles)))
abline(a = 0, b = 1)

#########################################################################################################

###### breaking up multitaper_est function

#input: X.t = time series of length N (possibly with NA values), W = analysis half bandwidth, K = number of tapers
#output: L x K matrix of tapers where L = length of time series without missing values, K = number of tapers
get_tapers <- function(X.t, W, K){
  N.long <- length(X.t)
  t.n <- 1:N.long
  missing.indices <- which(is.na(X.t))
  t.n[which(is.na(X.t))] <- NA
  
  
  dist.mat <- rdist(na.omit(t.n))
  
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  print("A matrix computed")
  
  eigdec <- eigs(A.prime, k = K, which = "LM")
  eig_vecs <- eigdec$vectors #get only the vectors
  print("tapers computed")
  
  if(K ==1){
    if (mean(Re(eig_vecs))<0){
      eig_vecs <- -eig_vecs
    }
  }
  
  if(K == 2  || K == 3){
    
    if (mean(Re(eig_vecs[,1]))<0){
      eig_vecs[,1] <- -eig_vecs[,1]
    }
    if (Re(eig_vecs[2,2] - eig_vecs[1,2])<0){
      eig_vecs[,2] <- -eig_vecs[,2]
    }
    
    if(K == 3){
      if (mean(Re(eig_vecs[,3]))<0){
        eig_vecs[,3] <- -eig_vecs[,3]
      }
    }
  }
  if(K >=4){
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
  }
  
  print("sign maintenance done")

  #"tapers" = eig_vecs, "e.values" = eigdec$values,
  
  return(eig_vecs)
} 

#input: X.t = time series of length N with any missing values and length L without, 
#       V.mat = L X K dimension taper matrix
#output: freqs = fourier frequencies
#       spectrum = spectral estimate

MT_spectralEstimate <- function(X.t, V.mat){
  im <- complex(real = 0, imaginary = 1)
  X.t <- X.t - mean(X.t, na.rm = TRUE) #demean
  N.long <- length(X.t)
  t.n <- 1:N.long
  missing.indices <- which(is.na(X.t))
  t.n[which(is.na(X.t))] <- NA
  t.n_m <- rdist(na.omit(t.n))[1,]
  
  ##use tapers to generate spectral estimate
  N <- length(na.omit(t.n))
  print(N)
  S.x.hat <- rep(NA, times = floor(N/2) + 1)
  freqs <- seq(0,0.5, length.out = floor(N/2) + 1)
  K <- dim(V.mat)[2]
  
  for(j in 1:length(freqs)){
    k.vec <- rep(NA,times = K)
    for(k in 1:K){
      inner.sum <- sum(V.mat[,k]*na.exclude(X.t)*exp(-im*2*pi*freqs[j]*t.n_m))
      k.vec[k] <- abs(inner.sum)^2
    }
    S.x.hat[j] <- mean(k.vec)
  }
  return(list("spectrum" = S.x.hat, "freqs" = freqs))
}




test.spec <- MT_spectralEstimate(X.t = X.t, V.mat = test)
plot(test.spec)





###################################

A.mat <- matrix(1:4, nrow = 2)
B.mat <- matrix(4:7, nrow = 2)

alpha.mat <- matrix(c(1,1,2,2), nrow = 2)
alpha.mat

first.try <- geigen::geigen(A.mat, B.mat)
after.scaling <- geigen::geigen(A.mat*alpha.mat, B.mat*alpha.mat)
first.try
after.scaling




















