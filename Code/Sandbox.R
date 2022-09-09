##### Sandbox #####

#this is a script to try things out before they go in real scripts

###bandpass vs. transfer function experiment

trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat <- bmat <- matrix(NA, ncol = 300, nrow = 9)
 

for(k in 1:9){
  tau = 2^k
for(i in 1:300){
  print(i)
  set.seed(i)
  #generate X.t
  X.t <- rnorm(2048,mean = 0, sd = 1)
  
  #calculate S.hat
  MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
  
  #calculate bandpass variance
  bpvar.vec[i] <- 2*(MTSE_full$spectrum[1]*(f[2] - 1/(4*tau)) + MTSE_full$spectrum[2]*(2/(1*tau) - f[2]))
  
  #calculate transfer function AVAR
  G.vec <- transfer.func(f, tau)
  G.vec[1] <- 1
  trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)

}
 tmat[k,] <- trfunc.vec
 bmat[k,] <- bpvar.vec
   
}

matplot(tmat)
lines(1/c(2^(1:9)))
matplot(bmat)
