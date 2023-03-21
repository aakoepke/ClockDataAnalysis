##### Sandbox #####

#this is a script to try things out before they go in real scripts


#test out the functions in AVARSpectralEstimate file

X.t <- X.t_missing <- rnorm(300)
#X.t_missing[c(100:200, 700:950, 1500:1600)] <- NA
X.t_missing[c(8:13,18,35,51:52)] <- NA
plot(X.t)

#get spectral estimate
spec.est <- spec.pgram(X.t, demean = TRUE)
test_AVAR <- AVAR_bpvar(spectral_est = spec.est$spec, taus = 2^(0:8))
test_AVAR

spec.est <- multitaper_est(X.t_missing,W=0.08, K = 5)

tapers = spec.est$tapers
t.n =  t.vec[which(!is.na(X.t_missing))]
X.t = na.omit(X.t_missing)
taus = 2^(0:6)



length(spec.est$spectrum)

t.vec <- 1:100


#var_AVAR_trfunc <- function(tapers, t.n, X.t, taus){

saved <- acf(X.t_missing, na.action = na.pass, lag.max = 100)
  #distance matrix of time points
  dist.mat <- rdist(t.n)
  #R.x matrix
  X.t <- X.t - mean(X.t)
  L <- length(X.t)
  #R.x <- matrix(saved$acf[dist.mat+1], nrow = L, ncol = L)
  #f vector
  f <- seq(0,0.5, length.out = floor(length(t.n))/2 + 1)
  
  
  #Calculate Cov matrix C
  N <- length(f)
  Cov.mat <- matrix(NA, nrow = N, ncol = N)
  W_matlist <- list()
  
  
  ##### Messing with Tapers #####
  for(i in 1:5){
    tapers[,i] <- tapers[,i]*exp(-im*2*pi*f[1]*t.n)
  }
  
  tapers_i <- tapers
  tapers_j <- spec.est$tapers
  
  for(i in 1:5){
    tapers_j[,i] <- tapers_j[,i]*exp(-im*2*pi*f[10]*t.n)
  }
  
  #######################
  
  for(i in 1:N){
    print(i)
    W_matlist[[i]] <- tapers*exp(-im*2*pi*f[i]*t.n)
  }
  
  for(i in 1:N){
  print(i)
  j = i
   while(j>=i & j <=N){
      W.star <- t(Conj(W_matlist[[i]]))
      W <- W_matlist[[j]]
      Cov.mat[i,j] <- (1/5)*sum(abs(W.star%*%W)^2) #frobenius norm 
      #Cov.mat[i,j] <- (1/5)*norm(W%*%W.star, type = "2") #2-norm
      j = j + 1
    }
}
  
  Cov.mat[lower.tri(Cov.mat)] <- t(Cov.mat)[lower.tri(Cov.mat)]
  
  out <- rep(NA, times = length(taus))
  #for(k in 1:length(taus)){
  k = 1
    #Calculate transfer function vector G_tau
    G_tau <- transfer.func(f,taus[k])
    G_tau[1] <- 1
    #G_tau^T * C * G_tau
    out[k] <- t(G_tau)%*%Cov.mat%*%G_tau 
  #}
  return(out)
#}
      






X.t <- TK95(N = 5224, alpha = 1)
X.t[c(100:500,1300:1400, 1700:1874, 3000:3450)] <- NA
X.t_sims_flk_gps[i,] <- X.t

#calculate S.hat
MTSE_full <- multitaper_est(X.t, W = 0.0007, K = 5)
