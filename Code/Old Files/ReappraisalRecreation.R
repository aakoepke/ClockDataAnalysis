### Recreating Figures and Concepts covered in:  ###########
######### Reappraisal of Freq. Domain Techniques ###########
####### for Assessing Feq. Stability Measurements ##########
################### Don Percival, 1987 #####################

#########   Section 3    ####################
###Figure 2: Frequency Drift Problem
##Random Walk AVAR vs. Bandpass Variance with Random Walk

#1. Create some data
set.seed(1)
N <- 128
y.t <- cumsum(sample(c(-1, 1), N, TRUE))
#add in a linear drift
y.t.LD <- y.t + 2 + 0.1*(1:N)
plot(y.t, type = "l")
plot(y.t.LD, type = "l")

#2. Calculate AVAR(y.t) and AVAR(y.t.LD)

AVAR.RW_y.t <- getAvars(N=N,y=y.t)
AVAR.RW_y.t.LD <- getAvars(N,y.t.LD)

#2. Difference the Data

z.t <- diff(y.t.LD) + 0.1
plot(z.t, type = "l")

#3. Estimate S_Z (spectrum of z process) using multitaper:
f.nyquist <- 1/2
N.prime <- N
f.j.prime <- seq(0,f.nyquist, length.out = N.prime/2 + 1)
X.tilde.prime <- z.t

t <- 0:(N.prime-1)
S.x.hat <- rep(NA, times = N.prime/2 + 1) #where S_x.hat(f'_j) will be saved

for(j in 0:N.prime/2){
  k.vec <- rep(NA,times = 6)
  for(k in 0:5){
    W.t <- h.k.prime(k = k, N=N, N.prime = N.prime)*X.tilde.prime
    inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*t*j/N.prime))
    k.vec[k + 1] <- abs(inner.sum)^2
  }
  S.x.hat[j] <- (1/length(k.vec))*sum(k.vec)
}

#4. S_Y(f) = S_Z/4sin^2(pi*f)

spec.y <- spec.z$spec/(4*sin(pi*spec.z$freq)^2)

par(mfrow = c(2,1))
plot(log(spec.z$freq), log(spec.z$spec), type = "l")
plot(log(spec.z$freq), log(spec.y), type = "l")



############ Section 4 ##############################
##Estimation of Parameters of Power Law Process: Gaussian White Noise
##sdf vs. AVAR approach

#We already know for WN(0,1), S_y(f) = 1

#1. Create some data: 1000 realizations of N=128 WN(0,1)
N <- 128
Y <- matrix(NA, nrow = 1000,ncol = 128)

for(i in 1:1000){
set.seed(i)
Y[i,] <- rnorm(N,mean = 0, sd = 1)
}

#2. Calculate AVAR estimates
AVAR.Y <- matrix(NA, nrow = 1000, ncol = 6)
for(i in 1:1000){
  AVAR.Y[i,] <- getAvars(N,Y[i,])$avarRes$avars
}

#3. Calculate S_y(f_k) estimates
Spec.Y <- matrix(NA, nrow = 1000, ncol = 64)
for(i in 1:1000){
  Spec.Y[i,] <- spec.pgram(Y[i,], plot = FALSE)$spec
}

#4. fit lines in log-log space to these estimates
v_k <- (0:5)*log(2)
alpha.ests_AVAR <- alpha.ests_spec <- rep(NA, times = 1000)

for(i in 1:1000){
 alpha.ests_AVAR[i] <-  -(lm(log(AVAR.Y[i,]) ~ v_k)$coefficients[2] + 1)
}

log_fk <- log(1:64/128)

for(i in 1:1000){
  alpha.ests_spec[i] <-  lm(log(Spec.Y[i,]) ~  log_fk)$coefficients[2] 
}


mean(alpha.ests_AVAR)
mean(alpha.ests_spec)








###Allan Variance####
#### Functions from Amanda ####
avar_fn=function(y,tau){
  n=length(y)
  
  div=seq(1,n,by = tau)
  
  M=length(div)-1 #number of groups
  
  groupmeans = numeric(M)
  for(i in 1:M){
    groupmeans[i]=mean(y[div[i]:(div[i+1]-1)])
  }
  
  1/(2*(M-1)) * sum(diff(groupmeans)^2)
}

getAvars=function(N, y){
  taus = 2^(0:5)
  
  avars=numeric(length(taus))
  
  for (i in 1:length(taus)){
    avars[i]=avar_fn(y,taus[i])
  }
  
  m1=data.frame(taus=taus,avars=avars)  
  fit=lm(log(sqrt(avars))~log(taus),data = m1)
  slope=as.numeric(fit$coefficients[2])
  int=as.numeric(fit$coefficients[1])
  
  avarRes=data.frame(taus=taus,avars=avars,N=N,slope=slope,int=int)  
  
  ##########################################
  # get SE
  ##########################################
  SEests=data.frame()
  
  new=data.frame(N=N,out=exp(int+slope*log(N)), type="AD")
  new2=data.frame(N=N,out=sd(y)/sqrt(N), type="SE")
  new3=data.frame(N=N,out=1/sqrt(N), type="true")
  SEests=bind_rows(new,new2,new3)
  
  return(list(avarRes=avarRes,SEests=SEests))
}
