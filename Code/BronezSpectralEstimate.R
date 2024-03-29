######################Bronez Spectral Estimate ###############

#parResults object from BronezTapers.R file is a matlist where
##each element has two objects: weights, eigenvalues

tapers2048 <- readRDS("Code/BronezResults2048.Rds")
tapers256 <- readRDS("Code/BronezResults.Rds")
#look at them
plot(tapers2048[[1]][,1])


##Here we generate the spectral estimate P(A) = sum_{i = 1}^{K}|w_k^*x|^2
taper.mat <- tapers256
N = length(taper.mat[[1]]$weights[,1])
N.fourier = floor(N/2) + 1
freq = seq(0,0.5, length.out = N.fourier)
set.seed(1990)
x.t <- rnorm(N)

spec.hat <- rep(NA, times = N.fourier)

for(i in 1:length(freq)){
  print(i)
  spec.hat[i] <- (1/3)*sum(abs(Conj(t(weights_scaled[[i]]))%*%x.t)^2)
}

##compare to Chave estimate
test_chave <- multitaper_est(X.t = x.t, W = 4/N, K = 3)

plot(freq, test_chave$spectrum, col = "blue", ylim = c(0,3))
points(freq, spec.hat) #/(2*(4/N)/3))
plot(freq, abs(test_chave$spectrum - spec.hat))
w.k <- tapers2048[[1]]$weights[,1]
w.k <- w.k*sqrt(0.001)
t(w.k)%*%w.k

weights_scaled = list()

for(i in 1:N.fourier){
  weights_scaled[[i]] <- tapers256[[i]]$weights/sqrt(0.01041667)
}

Cov.mat_chave <- matrix(NA, nrow = N.fourier, ncol = N.fourier)

for(i in 1:N.fourier){
  print(i)
  j = 1
  while(j <= i){
    Cov.mat_chave[i,j] <- norm(Conj(t(test_chave$tapers*exp(-im*2*pi*freq[i]*1:N)))%*%(test_chave$tapers*exp(-im*2*pi*freq[j]*1:N)), type = "2") #norm(V_star%*%exp(im*2*pi*freq[i]*dist.mat)%*%exp(-im*2*pi*freq[j]*dist.mat)%*%V, type = "2")*(A.size/numTapers)^2
    j = j+1
  }
}

Cov.mat_chave[upper.tri(Cov.mat_chave)] <- t(Cov.mat_chave)[upper.tri(Cov.mat_chave)]

Cov.mat_bronez[1:5,1:5]
Cov.mat_chave[1:5,1:5]
norm(Cov.mat_bronez - Cov.mat_chave, type = "2")

for(i in 1:length(taus)){
G_tau <- transfer.func(freq,tau = taus[i]) #change the tau value to get different vectors
G_tau[1] <- 1
#calculate variance for the AVAR estimate at the given tau
cov.mat.calc <- t(G_tau)%*%(Cov.mat_bronez)%*%G_tau*(0.5/256)^2
sample.var <- var(tmat[i,])
print(abs(cov.mat.calc-sample.var)/sample.var)
}

