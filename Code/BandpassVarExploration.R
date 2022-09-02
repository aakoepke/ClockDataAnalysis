#### Band-Pass Variance vs. transfer function #########

## calculating BPV estimates for known processes ##

library(fda) #functional boxplot



###bandpass vs. transfer function experiment: White Noise

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

taus <- c(2^(1:9),floor(N/3))

matplot(tmat)

#transfer function

plot(log10(taus), log10(tmat[,1]), type = "p", pch = 19, col = "red", ylim = c(-4,1),main =  "transfer function")

for(i in 1:300){
  points(log10(taus), log10(tmat[,i]), col = "red")
}

lines(log10(taus), log10(1/taus), col = "blue", lwd = 2)


#bandpass variance

plot(log10(taus), log10(bmat[,1]), type = "p", pch = 19, col = "red", ylim = c(-4,1), main = "bandpass variance")

for(i in 1:300){
  points(log10(taus), log10(bmat[,i]), col = "red")
}

lines(log10(taus), log10(1/taus), col = "blue", lwd = 2)


fbplot(log10(tmat), x = log10(taus), xlim = c(0.25,3), main = "Functional Boxplot for Transfer Function")
lines(log10(taus), log10(1/taus), col = "green", lwd = 2)

fbplot(log10(bmat), x = log10(taus), xlim = c(0.25,3), main = "Function Boxplot for Bandpass Variance")
lines(log10(taus), log10(1/taus), col = "green", lwd = 2)



##try another process




##write out comparison



