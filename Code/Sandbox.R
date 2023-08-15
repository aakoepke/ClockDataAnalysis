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


#####################
##what if we took the sign maintenance out of the estimate?

get_tapers <- function(t.n, W, K){
#t.n <- 1:100
  dist.mat <- rdist(na.omit(t.n))
#W = 0.025 
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  print("A matrix computed")
  
  eigdec <- eigs(A.prime, k = K, which = "LM")
  #eigdec <- eigen(A.prime)
  eig_vecs <- eigdec$vectors[,1:K] #get only the vectors
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
  
  return(eig_vecs)
} 


X.t <- rnorm(1524)
X.t[c(300:799)] <- NA
t.n <- 1:1524
t.n[c(300:799)] <- NA
V.eigs <- get_tapers(t.n, W = 0.00390625, K = 5)
V.eigen <- get_tapers(t.n, W = 0.00390625, K = 5)
plot(V.eigen[,3])

test_eigs <- MT_spectralEstimate(X.t, V.mat = V.eigs)
test_eigen <- MT_spectralEstimate(X.t, V.mat = V.eigen)
plot(test_eigs$spectrum)
lines(test_eigen$spectrum)
sum(abs(test_eigs$spectrum - test_eigen$spectrum))

t(V.eigs[,3])%*%V.eigen[,3]



###

parResult = readRDS("Code/BronezResults.Rds")







