##### Sandbox #####

#this is a script to try things out before they go in real scripts


#test out the functions in AVARSpectralEstimate file

X.t <- X.t_missing <- rnorm(2048)
X.t_missing[c(100:200, 700:950, 1500:1600)] <- NA

plot(X.t)

#get spectral estimate
spec.est <- spec.pgram(X.t, demean = TRUE)
test_AVAR <- AVAR_bpvar(spectral_est = spec.est$spec, taus = 2^(0:8))
test_AVAR

spec.est <- multitaper_est(X.t_missing,W=0.005, K = 5)


length(spec.est$spectrum)

t.vec <- 1:2048
test_var <- var_AVAR_trfunc(tapers = spec.est$tapers, t.n =  t.vec[which(!is.na(X.t_missing))], X.t = na.omit(X.t_missing),taus = 2^(0:8))

tapers = spec.est$tapers
t.n =  t.vec[which(!is.na(X.t_missing))]
X.t = na.omit(X.t_missing)
taus = 2^(0:8)



#var_AVAR_trfunc <- function(tapers, t.n, X.t, taus){


  #distance matrix of time points
  dist.mat <- rdist(t.n)
  #R.x matrix
  X.t <- X.t - mean(X.t)
  L <- length(X.t)
  R.x <- matrix(NA, nrow = L, ncol = L)
  for(i in 1:L){
    for(j in 1:L){
      R.x[i,j] <- saved$acf[abs(t.n[i] - t.n[j]) + 1]
    }
  }
  
  #f vector
  f <- seq(0,0.5, length.out = floor(length(t.n))/2 + 1)
  
  
  #Calculate Cov matrix C
  N <- length(f)
  Cov.mat <- matrix(NA, nrow = N, ncol = N)
for(i in 1:N){
  print(i)
   for(j in 1:N){
      W.star <- Conj(t(tapers)%*%exp(-im*2*pi*f[i]*dist.mat))
      W <- t(t(tapers)%*%exp(-im*2*pi*f[j]*dist.mat))
      Cov.mat[i,j] <- sum(abs(W.star%*%R.x%*%W)^2)
    }
  }
  
  out <- rep(NA, times = length(taus))
  #for(k in 1:length(taus)){
  k = 1
    #Calculate transfer function vector G_tau
    G_tau <- transfer.func(f,taus[k])
    
    #G_tau^T * C * G_tau
    out[k] <- t(G_tau)%*%Cov.mat%*%G_tau 
  #}
  return(out)
#}
      
      
#R.x matrix
      
saved <- acf(X.t_missing, na.action = na.pass, lag.max = 2048)

length(saved$acf)

saved$acf[dist.mat[i,j]]
length(t.n)
saved$acf[abs(t.n[i] - t.n[j])]
