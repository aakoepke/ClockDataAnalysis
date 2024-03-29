
#########################################
#######           Libraries         #####
#########################################

library(tidyverse)
library(RSpectra) #eigensolving
library(fields) #dist.mat
library(RobPer) #flicker noise
library(arfima) #acvf for arfima

###########################
#### things we'll need ####
###########################

## i
im <- complex(real = 0, imaginary = 1)

## MTSE for gappy data function 
multitaper_est <- function(X.t, W, K){
  X.t <- X.t - mean(X.t, na.rm = TRUE) #demean
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
  
  # if(K ==1){
  #   if (mean(Re(eig_vecs))<0){
  #     eig_vecs <- -eig_vecs
  #   }
  # }
  # 
  # if(K == 2  || K == 3){
  #   
  #   if (mean(Re(eig_vecs[,1]))<0){
  #     eig_vecs[,1] <- -eig_vecs[,1]
  #   }
  #   if (Re(eig_vecs[2,2] - eig_vecs[1,2])<0){
  #     eig_vecs[,2] <- -eig_vecs[,2]
  #   }
  #   
  #   if(K == 3){
  #     if (mean(Re(eig_vecs[,3]))<0){
  #       eig_vecs[,3] <- -eig_vecs[,3]
  #     }
  #   }
  # }
  # if(K >=4){
  #   #some sign maintenance
  #   for(i in seq(1,K,by = 2)){
  #     if (mean(Re(eig_vecs[,i]))<0){
  #       eig_vecs[,i] <- -eig_vecs[,i]
  #     }
  #   }
  #   
  #   for(i in seq(2,K-1,by = 2)){
  #     if (Re(eig_vecs[2,i] - eig_vecs[1,i])<0){
  #       eig_vecs[,i] <- -eig_vecs[,i]
  #     }
  #   }
  # }
  # 
  # print("sign maintenance done")
  
  ##use tapers to generate spectral estimate
  N <- length(na.omit(t.n))
  S.x.hat_MD <- rep(NA, times = floor(N/2) + 1)
  #freqs <- seq(0,0.5, length.out = floor(N/2) + 1)
  
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


######## the multitaper_est function broken into two parts #####
####### 1. get_tapers() ########
####### 2. calculate estimate ######

#input: t.n = time points of length N (possibly with NA values), W = analysis half bandwidth, K = number of tapers
#output: L x K matrix of tapers where L = length of time series without missing values, K = number of tapers
get_tapers <- function(t.n, W, K){
  dist.mat <- rdist(na.omit(t.n))
  
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  print("A matrix computed")
  
  eigdec <- eigs(A.prime, k = K, which = "LM")
  eig_vecs <- eigdec$vectors #get only the vectors
  print("tapers computed")
  
  # if(K ==1){
  #   if (mean(Re(eig_vecs))<0){
  #     eig_vecs <- -eig_vecs
  #   }
  # }
  # 
  # if(K == 2  || K == 3){
  #   
  #   if (mean(Re(eig_vecs[,1]))<0){
  #     eig_vecs[,1] <- -eig_vecs[,1]
  #   }
  #   if (Re(eig_vecs[2,2] - eig_vecs[1,2])<0){
  #     eig_vecs[,2] <- -eig_vecs[,2]
  #   }
  #   
  #   if(K == 3){
  #     if (mean(Re(eig_vecs[,3]))<0){
  #       eig_vecs[,3] <- -eig_vecs[,3]
  #     }
  #   }
  # }
  # if(K >=4){
  #   #some sign maintenance
  #   for(i in seq(1,K,by = 2)){
  #     if (mean(Re(eig_vecs[,i]))<0){
  #       eig_vecs[,i] <- -eig_vecs[,i]
  #     }
  #   }
  #   
  #   for(i in seq(2,K-1,by = 2)){
  #     if (Re(eig_vecs[2,i] - eig_vecs[1,i])<0){
  #       eig_vecs[,i] <- -eig_vecs[,i]
  #     }
  #   }
  # }
  # 
  # print("sign maintenance done")
  
  #"tapers" = eig_vecs, "e.values" = eigdec$values,
  
  return(list("tapers" = eig_vecs, "e.values" = eigdec$values))
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


MT_spectralEstimate_freqs <- function(X.t, freqs, V.mat){
  im <- complex(real = 0, imaginary = 1)
  X.t <- X.t - mean(X.t, na.rm = TRUE) #demean
  N.long <- length(X.t)
  t.n <- 1:N.long
  missing.indices <- which(is.na(X.t))
  t.n[which(is.na(X.t))] <- NA
  t.n_m <- rdist(na.omit(t.n))[1,]
  
  ##use tapers to generate spectral estimate
  N <- length(na.omit(t.n))
  S.x.hat <- rep(NA, times = length(freqs))
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


########function to calculate spectrum and covariance matrix 

spectralEstWithUnc=function(x.t,t.vec,numTapers){
  N <- length(t.vec)
  
  V.mat <- get_tapers(t.vec, W = 7/N, K = numTapers)
  MTSE_full <- MT_spectralEstimate(x.t, V.mat$tapers) #new function, calculates just spectrum
  N <- length(t.vec)
  N.fourier <- floor(N/2) + 1
  freq <- seq(0,0.5, length.out = N.fourier)
  
  delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
  
  ### calculate the covariance matrix 
  Cov.mat_chave <- matrix(NA, nrow = N.fourier, ncol = N.fourier)
  
  for(i in 1:N.fourier){
    j = 1
    while(j <= i){
      Cov.mat_chave[i,j] <- norm(Conj(t(V.mat$tapers*exp(-im*2*pi*freq[i]*t.vec)*(1/sqrt(numTapers))))%*%(V.mat$tapers*exp(-im*2*pi*freq[j]*t.vec)*(1/sqrt(numTapers))), type = "2") 
      j = j+1
    }
  }
  
  Cov.mat_chave[upper.tri(Cov.mat_chave)] <- t(Cov.mat_chave)[upper.tri(Cov.mat_chave)]
  
  
  return(list(freq=freq,
              spec.hat=MTSE_full$spectrum,Cov.mat=Cov.mat_chave))
  
}

########## calculate avar and unc from results of spectralEstWithUnc



#input: X.t = time series of length N with any missing values and length L without, 
#       V.mat = L X K dimension taper matrix
#output: freqs = fourier frequencies
#       spectrum = spectral estimate

MT_spectralEstimate_fft <- function(X.t, V.mat){
  
  ##use tapers to generate spectral estimate
  N <- length(na.exclude(X.t))
  N.fourier <- floor(N/2) + 1
  S.x.hat <- rep(NA, times = N.fourier)
  freqs <- seq(0,0.5, length.out = N.fourier)
  K <- dim(V.mat)[2]
  S.k.mat <- matrix(NA,nrow = K, ncol = N.fourier)
  
    for(k in 1:K){
      spec.vec <- fft(taperMatrix[,k]*na.exclude(X.t))[1:N.fourier]
      S.k.mat[k,] <- abs(spec.vec)^2
    }
  
    S.x.hat <- apply(S.k.mat, MARGIN = 2, FUN = mean)
  
  return(list("spectrum" = S.x.hat, "freqs" = freqs))
}


######################################################################################################
###### Calculate sigma.hat_tau using transfer function method (see "reappraisal...", page 70) ########
######################################################################################################

#G_tau(f) function
transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}

###########   Allan Variance Calculation    #############
######## (eq'n at bottom of p. 70 of reappraisal) #######
#input: spectral estimate (as a vector), taus (as a vector) where you want the AVAR calculated
#output: a vector of the AVAR estimates

AVAR_trfunc <- function(spectral_est, taus){
  out <- rep(NA, times = length(taus))
  f <- seq(0,0.5,length.out = length(spectral_est))
  
  for(i in 1:length(taus)){
    G.vec <- transfer.func(f, taus[i])
    G.vec[1] <- 0
    out[i] <- f[2]*sum(G.vec*spectral_est)
  }
  
  return(out)
}


AVAR_trfunc_withUnc <- function(spectral_est, taus,Cov.mat_chave){
  f <- seq(0,0.5,length.out = length(spectral_est))
  
  cov.mat=avar=numeric(length(taus))
  
  for(i in 1:length(taus)){
    G.vec <- transfer.func(f,tau = taus[i]) 
    G.vec[1] <- 0 
    
    avar[i]=f[2]*sum(G.vec*spectral_est)
    
    #calculate variance for the AVAR estimate at the given tau
    cov.mat[i] <- t(G.vec)%*%(Cov.mat_chave)%*%G.vec*(f[2])^2
    
  }
  return(list(avar=avar,avarVar=cov.mat))
}

########## AVAR/OVAR Functions #######


#calculates regular AVAR for a given tau and data (y)
avar_fn <- function(y,tau){
  n=length(y)
  
  div=seq(1,n,by = tau)
  
  M=length(div)-1 #number of groups
  
  groupmeans = numeric(M)
  for(i in 1:M){
    groupmeans[i]=mean(y[div[i]:(div[i+1]-1)])
  }
  
  1/(2*(M-1)) * sum(diff(groupmeans)^2)
}


#calculates OVAR for a given averaging factor (m) and data (y)
overlapping_avar_fn <- function(y,m){
  
  M=length(y)
  
  outer.sum = 0
  for(j in 1:(M-2*m+1)){
    sum = 0
    for(i in j:(j + m - 1)){
      sum = sum+ y[i+m] - y[i]
    }
    outer.sum = outer.sum + sum^2
  }
  
  out <- 1/(2*m^2*(M - 2*m + 1))*outer.sum
  
  return(out)
}


# #calculates avar uncertainty for a given averaging factor (m) and data (y)
# # m=tau (can be vector)
# # avar can be vector
# avar_unc <- function(N,m,avar){
#   
#   # N=length(y)
#   
#   edf=((3*(N-1)/(2*m))-((2*(N-2))/N))*((4*m^2)/(4*m^2+5))
# 
#   lower=avar*edf/qchisq(p = c(.975),df = edf) 
#   upper=avar*edf/qchisq(p = c(.025),df = edf)
#   ## hard coding these for 95% intervals, but scientists usually use 1-sigma intervals, and below is Dave's code 
#   ## (not sure where those numbers come from exactly...)
#   ## osigMax(j) = sqrt(edf/chi2inv(.1573,edf));
#   ## osigMin(j) = sqrt(edf/chi2inv(.8427,edf));
#   
#   UncDF=data.frame(tau=m,avar=avar,lower=lower, upper=upper)
#   
#   return(UncDF)
# }
# uncertainty for overavar

avar_CI <- function(CI.level,noise_type = "white noise", avar_type, avars, taus,N){
  
  a <- (1-CI.level)/2
  s.2=avars
  
  edf <- rep(NA, times = length(taus))
  i=1
  
  if(noise_type == "white noise"){
    for(m in taus){
      edf[i] <- ((3*(N-1)/(2*m)) - (2*(N-2)/N))*(4*m^2)/(4*m^2 + 5)
      i=i+1
    }
  }
  
  if(avar_type == "regular"){
    CI.limits <- bind_rows("lower" = s.2 - s.2/N, "upper" = s.2 +  s.2/N)
  }else{
    CI.limits <- bind_rows("lower" = s.2*edf/qchisq(1-a,edf),"upper" = s.2*edf/qchisq(a, edf) )
  }
  return(CI.limits)
}


#inputs: length of data (N), data (y), 
#        vector of tau values at which you want AVAR/OVAR calculated (taus)
#outputs:avarRes = data frame with taus, avar estimates, ovar estimates, slope/int fit
#        SEests = SE estimates from the data and the truth? For white noise anyway, I'm not sure about this part
getAvars <-  function(N, y,taus){
  
  avars=numeric(length(taus))
  overlapping_avars= numeric(length(taus))
  
  for (i in 1:length(taus)){
    avars[i]=avar_fn(y,taus[i])
    overlapping_avars[i] = overlapping_avar_fn(y,taus[i])
  }
  
  # m1=data.frame(taus=taus,avars=avars)  
  # fit=lm(log(sqrt(avars))~log(taus),data = m1)
  # slope=as.numeric(fit$coefficients[2])
  # int=as.numeric(fit$coefficients[1])
  
  avarRes=data.frame(taus=taus,avars=avars, overavars=overlapping_avars) #,N=N,slope=slope,int=int)  
  
  ##########################################
  # get SE
  ##########################################
  # m2=data.frame(taus=taus,oavars=overlapping_avars)  
  # fit2=lm(log(sqrt(oavars))~log(taus),data = m2)
  # slope2=as.numeric(fit2$coefficients[2])
  # int2=as.numeric(fit2$coefficients[1])
  # 
  # SEests=data.frame()
  # 
  # onew=data.frame(N=N,out=exp(int2+slope2*log(N)), type="OAD")
  # new=data.frame(N=N,out=exp(int+slope*log(N)), type="AD")
  # new2=data.frame(N=N,out=sd(y)/sqrt(N), type="SE")
  # new3=data.frame(N=N,out=1/sqrt(N), type="true")
  # SEests=bind_rows(new,new2,new3,onew)
  # 
  return(list(avarRes=avarRes)) #,SEests=SEests
}

##calculates the true allan variance for ARFIMA(0,d,0)
#inputs: N.tau - max tau value you are looking for
#        d - parameter
#        sig.2.a - variance of white noise process in arfima
#output: vector of length N.tau-1 with  

tavar_ARFIMA <- function(N.tau,d, sig.2.a){
  rho.vec <- tacvfARFIMA(phi = 0, theta = 0, dfrac = d, maxlag = 2*N.tau)
  corr.vec <- rho.vec/max(rho.vec) #normalize
  taus <- 2:N.tau
  
  sum.vec <- rep(NA, times = N.tau-1)
  for(k in 2:N.tau){
    #print(k)
    total <- 0
    for(i in 1:(k - 1)){
      #print(i) 
      total = total + i*(2*corr.vec[k-i + 1] - corr.vec[i + 1] - corr.vec[2*k-i + 1])
    }
    sum.vec[k-1] <- total
  }
  numerator <-(taus*(rep(corr.vec[1],times = N.tau-1) - corr.vec[3:(N.tau+1)]) + sum.vec)*gamma(1-2*d)
  denom <- (taus*gamma(1-d))^2
  return(numerator/denom)
}


##calculates the lomb-scargle periodogram for data x.t at frequencies f
#inputs: x.t - data (possibly with NA values)
#        f - frequencies at which we calculate the lomb-scargle periodogram
#output: spectral estimate

#lomb-scargle function

lomb_scargle <- function(x.t,f){
  
  ## calculates the Lomb-Scargle Periodogram for data x.t at frequencies f##

  N <- length(x.t)
  L <- length(f)
  t.vec <- 1:N
  t.vec[which(is.na(x.t))] <- NA
  x.missing <- na.exclude(x.t)
  t.missing <- na.exclude(t.vec)
  
  lsperio <- rep(NA, times = L)
  
  for(i in 1:L){
    x.centered <- x.missing - mean(x.missing)
    x.var <- var(x.missing)
    tau.value <- tau.shift(f[i], t.missing)
    c.vec <- cos(2*pi*(f[i]*(t.missing - tau.value)))
    s.vec <- sin(2*pi*(f[i]*(t.missing - tau.value)))
    lsperio[i] <- (1/(2*x.var))*((x.centered%*%c.vec)^2/sum(c.vec^2) + 
                                   (x.centered%*%s.vec)^2/sum(s.vec^2))
  }
  
  return(lsperio)
}

#tau function
tau.shift <- function(f,t){
  (1/(4*pi*f))*atan(sum(sin(4*pi*f*t))/sum(cos(4*pi*f*t)))
}


#### generating pink noise from Tara

generate_pink_noise <- function(length, fs) {
  pink_noise <- cumsum(rnorm(length))
  pink_noise <- pink_noise - mean(pink_noise)
  pink_noise <- pink_noise / max(abs(pink_noise))
  return(pink_noise)
}
