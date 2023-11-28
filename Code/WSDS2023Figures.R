############## WSDS Talk 2023 Figures #####################
###########################################################

###### Integrated Spectrum #########
N <- 1024
freqs <- seq(0,0.5, length.out = N/2 + 1)

true.spec <- function(f){
  y <- 1:length(f)
  y[f<0.1] <- 0.2
  y[f>=0.1 & f < 0.3] <- 0.5
  y[f>=0.3 & f < 0.4] <- 0.2
  y[f>= 0.4 & f <=0.45]<- 0.1
  y[f>0.45] <- 0.2
  return(y)
}

truth <- true.spec(freqs)

om = 0.005
integrated.spec <- 1:length(freqs)

for(i in 1:(length(freqs))){
  tmp <-  integrate(true.spec, lower = freqs[i] - om, upper = freqs[i] + om)
  integrated.spec[i] <- tmp$value
}

#plots
par(mar = c(5,5,2,2), cex.lab = 1.8)
plot(freqs, truth, pch = 19 , ylab = "Spectral Density", xlab = "Frequency", type = "l", lwd = 2)
plot(freqs, integrated.spec, main = expression(paste(omega, " = ", 0.005)), ylab = "Integrated Spectrum", xlab = "frequency", type = "l", lwd = 2)




###### How eigenvalues change with analysis bandwidth #############
t.n <- 1:2500
t.n[c(100:200, 1050:1350,1800:1849)] <- NA
N <- length(na.omit(t.n))

t.n <- 1:9000
t.n[c(100:200, 1050:1350,1800:1849, 3000:7000, 6010:6500)] <- NA
N <- length(na.omit(t.n))
N

all.fw <- seq(0.0001,0.005,length.out = 51)

top.eight <- matrix(NA, nrow = length(all.fw), ncol = 15)

for(i in 1:length(all.fw)){
  print(i)
  tmp <- get_tapers(t.n, W = all.fw[i], K = 15) 
  top.eight[i,] <- tmp$e.values
}


matplot(N*all.fw, top.eight, type = "l", lwd = 2, xlab = "NW", ylab = expression(lambda[k])) 
abline(h = 0.95)
abline(v = 1)
abline(v = N*f.w2[min(which(top.eight2[,1]>0.95))])

lambda.crosspoints <- rep(NA, times = 15)

for(i in 1:15){
  lambda.crosspoints[i] <- N*all.fw[min(which(top.eight[,i]>0.95))]
}

lambda.close.to.one <- rep(NA, times = 6)

for(i in 1:6){
  lambda.close.to.one[i] <- sum(lambda.crosspoints <= 2*i - 1)
}

plot(2*1:6-1,lambda.close.to.one, ylim = c(0,15), xlim = c(0,15), pch = 19, xlab = "2NW-1", 
     ylab = expression(paste("Number of ", lambda, " > 0.95")))
abline(a = 0, b = 1)





