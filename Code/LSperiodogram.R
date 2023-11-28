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

##code myself using VanderPlas2018
N <- 2048
g <- rnorm(N)
missing.points <- c(20:40, 60:80)
g[missing.points] <- NA
t.vec <- 1:N
t.vec[missing.points] <- NA
g.missing <- g %>% as_tibble() %>% drop_na()
t.missing <- t.vec %>% as_tibble() %>% drop_na()

lsp.check <- lsp(x = x.t, times = t.vec, plot = FALSE)
f <- lsp.check$scanned

lomb_scargle <- function(x.t,f){
  
  ## calculates the Lomb-Scargle Periodogram for data x.t at frequencies f##
  
  N <- length(x.t)
  t.vec <- 1:N
  t.vec[which(is.na(x.t))] <- NA
  x.missing <- na.omit(x.t)
  t.missing <- na.omit(t.vec)
  
  lsperio <- rep(NA, times = length(f))
  
  for(i in 1:length(f)){
    x.centered <- x.missing - mean(x.missing)
    x.var <- var(x.missing)
    tau.value <- tau.shift(f[i], t.missing)
    c.vec <- cos(2*pi*(f[i]*(t.missing - tau.value)))
    s.vec <- sin(2*pi*(f[i]*(t.missing - tau.value)))
    lsperio[i] <- (1/(2*g.var))*((g.centered%*%c.vec)^2/sum(c.vec^2) + 
                                   (g.centered%*%s.vec)^2/sum(s.vec^2))
  }
  
  return(list("lsp" = lsperio, "freqs" = f))
}


plot(lsperio, type = "l", ylim = c(0,1.5))



#tau function
tau.shift <- function(f,t){
  (1/(4*pi*f))*atan(sum(sin(4*pi*f*t))/sum(cos(4*pi*f*t)))
}



