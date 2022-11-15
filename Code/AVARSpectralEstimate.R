##### Calculate sigma.hat_tau using the Spectral Estimate ########

im <- complex(real = 0, imaginary = 1)
#Calculate sigma.hat_tau using transfer function

transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}

#input: spectral estimate, taus where you want the AVAR calculated
#output: a vector of the AVAR estimates
AVAR_trfunc <- function(spectral_est, taus){
out <- rep(NA, times = length(taus))
f <- seq(0,0.5,length.out = length(spectral_est))

for(i in 1:length(taus)){
G.vec <- transfer.func(f, taus[i])
G.vec[1] <- 1
out[i] <- f[2]*sum(G.vec*spectral_est)
}

return(out)
}


#Calculate sigma.hat_tau using bandpass variance

#input: spectral estimate, taus where you want the AVAR calculated
#output: a vector of the AVAR estimates
AVAR_bpvar <- function(spectral_est, taus){
  out <- rep(NA, times = length(taus))
  f <- seq(0,0.5,length.out = length(spectral_est))
  delta.f <- f[2]
  
  for(i in 1:length(taus)){
  tau <- taus[i]
  if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
    f.min.index <- which(f == 1/(4*tau))
    f.max.index <- which(f == 1/(2*tau))
    out[i] <- 4*delta.f*(sum(spectral_est[f.min.index:(f.max.index - 1)]))
  }
  else{
    f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
    f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
    if(f.min.index == f.max.index){
      out[i] <- 4*(spectral_est[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + spectral_est[f.min.index]*(1/(2*tau) - f[f.min.index]))
    }
    else{
      out[i] <- 4*delta.f*(sum(spectral_est[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
    }
  }
}

  return(out)
}


#Variance of sigma.hat_tau using transfer function

var_AVAR_trfunc <- function(tapers, t.n, X.t, taus){
  #distance matrix of time points
  dist.mat <- rdist(t.n)
  #R.x matrix
  X.t <- X.t - mean(X.t)
  L <- length(X.t)
  R.x <- matrix(NA, nrow = L, ncol = L)
  for(i in 1:L){
    for(j in 1:L){
      R.x[i,j] <- cor(X.t[i],X.t[j])
    }
  }
  
  #f vector
  f <- seq(0,0.5, length.out = floor(length(t.n))/2 + 1)
 
  
  #Calculate Cov matrix C
  N <- dim(tapers)[1]
  
  for(i in 1:N){
    for(j in 1:N){
      W.star <- Conj(t(tapers)%*%exp(-im*2*pi*f[i]*dist.mat))
      W <- t(tapers)%*%exp(-im*2*pi*f[j]*dist.mat)
      C[i,j] <- norm(W.star%*%R.x[i,j]%*%W, type = "F")
    }
  }
    
for(k in 1:length(taus)){
#Calculate transfer function vector G_tau
  G_tau <- transfer.func(f,tau)

#G_tau^T * C * G_tau
  out[k] <- t(G_tau)%*%C%*%G_tau 
}
  return(out)
}

#Variances of sigma.hat_tau using bandpass variance

