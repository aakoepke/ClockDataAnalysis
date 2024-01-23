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

#undo the normalization
# where unnormalized P(f) = L-S * R, 
#where R = sum(x.t^2) 
plot(lsperio$power*sum(x.t^2, na.rm = T), type = "l")

spec.pgram(x.t)

##code myself using https://www.mathworks.com/help/signal/ref/plomb.html#lomb
N <- 2048
g <- rnorm(N)
missing.points <- c(20:400, 600:800)
g[missing.points] <- NA
t.vec <- 1:N
t.vec[missing.points] <- NA


lsp.check <- lsp(x = g, times = t.vec, plot = TRUE)
f <- seq(0,0.5, length.out = floor(length(na.omit(t.vec))/2) + 1)


lsp.mine <- lomb_scargle(g,f[-1])

plot(lsp.mine$freqs, lsp.mine$lsp, type = "l")
check <- spec.pgram(g)


########### Functions ################

#lomb-scargle function

lomb_scargle <- function(x.t,f){
  
  ## calculates the Lomb-Scargle Periodogram for data x.t at frequencies f##
  
  N <- length(x.t)
  L <- length(f)
  t.vec <- 1:N
  t.vec[which(is.na(x.t))] <- NA
  x.missing <- na.omit(x.t)
  t.missing <- na.omit(t.vec)
  
  lsperio <- rep(NA, times = L)
  
  for(i in 1:L){
    x.centered <- x.missing - mean(x.missing)
    x.var <- var(x.missing)
    tau.value <- tau.shift(f[i], t.missing)
    c.vec <- cos(2*pi*(f[i]*(t.missing - tau.value)))
    s.vec <- sin(2*pi*(f[i]*(t.missing - tau.value)))
    lsperio[i] <- (1/(2*g.var))*((x.centered%*%c.vec)^2/sum(c.vec^2) + 
                                   (x.centered%*%s.vec)^2/sum(s.vec^2))
  }
  
  return(lsperio)
}

#tau function
tau.shift <- function(f,t){
  (1/(4*pi*f))*atan(sum(sin(4*pi*f*t))/sum(cos(4*pi*f*t)))
}



