##### Sandbox #####

#this is a script to try things out before they go in real scripts


###look at distribution of spectral-based AVAR estimates
tmat <- readRDS(file = "Code/Paper1/Results/tmat050523_W4_K6_N2048_300sims_WhiteNoiseNoGaps.Rds")


##fit with ML

##we'll do it numerically
fun.to.optimize <- function(x,k){
  n <- length(x)
  log.like <- -(k*n/2)*log(2) - n*log(gamma(k/2)) + (k/2 - 1)*sum(log(x)) - sum(x)/2
}

mle.chisq <- optim(par = 3, fn = fun.to.optimize, x = tmat[10,], lower = 1, upper = 50, method = "Brent")
mle.chisq$par

k.seq <- seq(1,30, by = 1)
x <- rchisq(300, df = 5)
plot(k.seq, fun.to.optimize(x,k.seq))



