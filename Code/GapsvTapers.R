################ Gaps vs. Tapers #####################
###  In this script we test the interplay between  ###
### the length of gaps and the tapers              ###
######################################################

numberOfSimulations = 500
N.long = 2048 + 1000
t.n_missing <- 1:N.long
#t.n_missing[c(100:500,1300:1400, 1700:1874, 2400:2722)] <- NA
t.n_missing[c(500:1499)] <- NA
N <- length(na.omit(t.n_missing))

setWnum = 12
setW = setWnum/N
setK = 13

##calculate tapers
V.mat <- get_tapers(t.n_missing, W = setW, K = setK)
V.mat$e.values

plot(t.n_missing, rep(0, times = N.long))

par(mfrow = c(4,1))
colors <- rainbow(setK)

for(i in 1:setK){
plot(V.mat$tapers[,i], col = colors[i], ylim = c(-0.06,0.06))
abline(h = 0)
}







