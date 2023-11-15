############# Lomb-Scargle Periodogram #############

library(lomb)

#calculate lomb-scargle periodogram of unevenly spaced data

#generate a time series with gaps
N <- 2048
set.seed(1)
x.t <- rnorm(N)
t.n <- 1:N
x.t[c(300:500, 700:900,1500:1750)] <- NA
t.n[c(300:500, 700:900,1500:1750)] <- NA

lsperio <- lsp(x = x.t, times = t.n)
length(lsperio$power)
