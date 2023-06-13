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

####################################
X.t <- rnorm(256)

test <- multitaper_est(X.t, W = 4/256, K = 5)
test$e.values


plot(clock_df - mean(clock_df, na.rm = TRUE))
var(clock_df, na.rm = TRUE)

plot(na.omit(clock_df))

acf(clock_df, na.action = na.pass)

diffedx2_clockdata <- diff(diff(clock_df))


















