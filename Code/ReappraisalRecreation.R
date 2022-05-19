### Recreating Figures and Concepts covered in:  ###########
######### Reappraisal of Freq. Domain Techniques ###########
####### for Assessing Feq. Stability Measurements ##########
################### Don Percival, 1987 #####################


###Exammple 1: Random Walk

#1. Create some data
set.seed(1)
n <- 1024
y.t <- cumsum(sample(c(-1, 1), n, TRUE))
#add in a linear drift
y.t.LD <- y.t + 2 + 0.1*(1:1000)
plot(y.t.LD, type = "l")

#2. Difference the Data

z.t <- diff(y.t.LD) + 0.1
plot(z.t, type = "l")

#3. Estimate S_Z (spectrum of z process):
spec.z <- spec.pgram(z.t)

#4. S_Y(f) = S_Z/4sin^2(pi*f)

spec.y <- spec.z$spec/(4*sin(pi*spec.z$freq)^2)

par(mfrow = c(2,1))
plot(log(spec.z$freq), log(spec.z$spec), type = "l")
plot(log(spec.z$freq), log(spec.y), type = "l")


##Example 2: Gaussian White Noise

#We already know for WN(0,1), S_y(f) = 1

#1. Create some data
set.seed(1)
n <- 1024
y.t.WN <- rnorm(n,mean = 0, sd = 1)

plot(y.t.WN, type = "l")

#2. Difference the Data

z.t.WN <- diff(y.t.WN)
plot(z.t.WN, type = "l")

#3. Estimate S_Z (spectrum of z process):
spec.z.WN <- spec.pgram(z.t.WN)

#4. S_Y(f) = S_Z/4sin^2(pi*f)

spec.y.WN <- spec.z.WN$spec/(4*sin(pi*spec.z.WN$freq)^2)

par(mfrow = c(2,1))
plot(log(spec.z.WN$freq), log(spec.z.WN$spec), type = "l")
plot(spec.z.WN$freq, spec.y.WN, type = "l")
abline(h = 1)


###WN Example ctd.##

#look at log(S_hat_y) vs. log(f)
plot(log(spec.z.WN$freq), log(spec.y.WN))
lin.fit <- lm(log(spec.y.WN) ~ log(spec.z.WN$freq))

lin.fit
abline(a = lin.fit$coefficients[1], b = lin.fit$coefficients[2])
plot(lin.fit)

#I'm definitely missing something


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

getAvars=function(N){
  y=rnorm(N,0,1) 
  taus = c(1,2,3,5,10,20)
  taus = c(taus,seq(40,N,by=100)) ## this goes too far, gives lots of NAs from avar_fn, fix later 
  
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

try1 <- getAvars(1024)

plot(log(try1$avarRes$taus), log(try1$avarRes$avars))
