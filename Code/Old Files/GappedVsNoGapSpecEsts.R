####comparing spectral estimate with gaps to no gaps ##########


####ar(p) process?
N = 2048
phi_vec <- c(-0.1,0.2, 0.3, -0.1)
sim <- arima.sim(model = list(ar = phi_vec, sd = 1), n = N)
plot(sim)
acf(sim)

N = 1024
freqs <- seq(0, 0.5, length.out = N/2 + 1)
s_f <- sapply(freqs, FUN =  arp.spectrum, phi_vec = phi_vec, sigma2_z = 1)
plot(freqs, s_f, type = "l")
tapers_ar1 <- get_tapers(t.n = 1:N, W = 4/N, K = 5)
plot(tapers_ar1$tapers[,4])


mtse_ar1 <- MT_spectralEstimate(X.t = sim, V.mat = tapers_ar1$tapers)
periodogram_ar1 <- spec.pgram(sim)
  
plot(freqs,log(mtse_ar1$spectrum))
lines(freqs, log(s_f), type = "l")
points(periodogram_ar1$freq, log(periodogram_ar1$spec), col = "blue")


###arfima(0,d,0)

arfima.spec <- function(f, sigma2.e, d){
  sigma2.e/((2*abs(sin(pi*f)))^(2*d))
}

N = 1024
freqs <- seq(0, 0.5, length.out = N/2 + 1)
plot(log(freqs[-1]), log(arfima.spec(freqs[-1], 1, d = 0.25)))


library(arfima)
x.t <- arfima.sim(n = N, model = list(dfrac = 0.25))
plot(x.t)

tapers_arfima <- get_tapers(t.n = 1:N, W = 4/N, K = 7)
mtse_arfima <- MT_spectralEstimate(X.t = x.t, V.mat = tapers_arfima$tapers)
periodogram_arfima <- spec.pgram(x.t, plot = FALSE)
plot(log(mtse_arfima$freqs[-1]), log(mtse_arfima$spectrum[-1]))
lines(log(freqs[-1]), log(arfima.spec(freqs[-1], 1, d = 0.25)))
points(log(periodogram_arfima$freq), log(periodogram_arfima$spec), col = "blue", pch = 19)

##add in gaps

N.long = 5000
x.t <- arfima.sim(n = N.long, model = list(dfrac = 0.49999))
x.t[c(1200:2000, 2300:2500, 3000:3500)] <- NA
t.vec <- 1:N.long
t.vec[c(1200:2000, 2300:2500, 3000:3500)] <- NA

N <- length(na.exclude(x.t))
tapers_arfima <- get_tapers(t.n = na.exclude(t.vec), W = 8/N, K = 7)
mtse_arfima <- MT_spectralEstimate(X.t = x.t, V.mat = tapers_arfima$tapers)
tapers_arfima$e.values
plot(tapers_arfima$tapers[,3])

tapers_concat_arfima <- get_tapers(t.n = 1:N, W = 8/N, K = 7)
tapers_concat_arfima$e.values
plot(tapers_concat_arfima$tapers[,1])
mtse_concat_arfima <- MT_spectralEstimate(X.t = na.exclude(x.t), V.mat = tapers_concat_arfima$tapers)

plot(log(mtse_arfima$freqs[-1]), log(mtse_arfima$spectrum[-1]))
lines(log(freqs[-1]), log(arfima.spec(freqs[-1], 1, d = 0.49)))
points(log(mtse_concat_arfima$freqs[-1]), log(mtse_concat_arfima$spectrum[-1]), col = "blue", pch = 19)



#############################################
############## Functions ####################
#############################################

ar1.spectrum <- function(phi, sigma2_z, f){
  sigma2_z/(1 - phi^2 - 2*phi*cos(2*pi*f))
}

arp.spectrum <- function(phi_vec, sigma2_z, f){
  im <- complex(real = 0, imaginary = 1)
  sigma2_z/abs(1 - sum(phi*exp(-im*2*pi*f*1:length(phi_vec))))^2
}

freqs <- seq(0,0.5,length.out = 1000)
s_f <- lapply(freqs, FUN =  arp.spectrum, phi_vec = phi_vec, sigma2_z = 1)
plot(freqs, s_f)



dist.mat <- rdist(1:10)
dist.mat
exp(im*2*pi*0.1*dist.mat)

sin(pi*dist.mat)

