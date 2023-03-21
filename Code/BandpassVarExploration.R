#### Band-Pass Variance vs. transfer function #########

## calculating BPV estimates for known processes ##
library(tidyverse)
library(fda) #functional boxplot
library(fields)
library(RSpectra)
library(magrittr)
library(RobPer) #flicker noise

## Bandpass Variance or Transfer Function comparison to AVAR
delta.f <- f[2] #delta f, distance between each frequency in the frequency grid

#bandpass variance
tau <-  floor(N/3) #1/3 the length of the dataset

f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))

bp.var <- delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )


#transfer function
transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}
G.vec <- transfer.func(f, tau)
G.vec[1] <- 1

f[2]*sum(G.vec*MTSE_full$spectrum)




#################################
###Example 1: White Noise #######
#################################

N <- 14500
trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat <- bmat <- matrix(NA, ncol = 300, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]


r = 3
for(k in c(2^(0:9), floor(N/3))){
  r = r + 1
  tau = k
  print(paste("r = ", r))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- rnorm(2048,mean = 0, sd = 1)
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(MTSE_full$spectrum[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + MTSE_full$spectrum[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }
      else{
        bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat[r,] <- trfunc.vec
  bmat[r,] <- bpvar.vec
}





###also Calculate AVAR#####

amat <- oamat <- matrix(NA, nrow = 300, ncol = 11)

for(i in 1:300){
  set.seed(i)
  print(i)
  #generate X.t
  X.t <- rnorm(2048,mean = 0, sd = 1)
  
  avar.calc <- getAvars(N,X.t, taus = c(2^(0:9), floor(N/3)))
  amat[i,] <- avar.calc$avarRes$avars
  oamat[i,] <- avar.calc$avarRes$overavars
}

#for plotting: bmat, tmat, amat/oamat


###### Plots #####
N <- 2048
taus <- c(2^(0:9),floor(N/3))

##tidy the data
tmat %<>% t()
bmat %<>% t()
dim(bmat)
three.mat <- rbind(oamat, tmat)
method_labels <- rep(c("avar", "tr"), each = 300)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

#sumDat=dat %>% group_by(method,tau) %>%
 # summarise(median = median(measurement),
 #          lower25=quantile(measurement,prob=.25),
#            upper75=quantile(measurement,prob=.75),
 #           min=min(measurement),
  #          max=max(measurement))

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9,21 )*rep(10^(-10:10), each = 9)

linedat_w=data.frame(tau=1:1000)
linedat_w$truth=1/linedat_w$tau#+linedat$tau
linedat_w$method=NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat_w,aes(tau,truth), size = 1.2)+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("Allan","Spectral"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here 
  ### basically specifies the x and y position in terms of the plot size in unit scale. 
  theme_bw(base_size = 20)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  annotation_logticks()+

  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  ggtitle("White Noise Comparison, No Gaps")


save(dat, file = "TidyDataforBoxplots.RData")


##### WN with gaps #####
N <- 2048
trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat_gps2 <- bmat_gps2 <- matrix(NA, ncol = 300, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]



### create tapers
t.n <- 1:2926
t.n <- t.n[-c(100:500,1300:1400, 1700:1874, 2000:2200)]
dist.mat <- rdist(t.n)
W <-  9/2048

A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
A.prime[row(A.prime) == col(A.prime)] <- W*2
eigdec <- eigs_sym(A.prime, k = 5, which = "LM")

#eigdec <- eigen(A.prime, symmetric = TRUE)

##changing signs of the vectors
K = 5 #number of sequences we want
eig_vecs <- eigdec$vectors[,1:K] #get only those vectors
for(i in seq(1,K,by = 2)){
  if (mean(Re(eig_vecs[,i]))<0){
    eig_vecs[,i] <- -eig_vecs[,i]
  }
}

for(i in seq(2,K-1,by = 2)){
  if (Re(eig_vecs[2,i] - eig_vecs[1,i])<0){
    eig_vecs[,i] <- -eig_vecs[,i]
  }
}


r = 0
for(k in c(2^(0:9), floor(N/3))){
  r = r + 1
  tau = k
  print(paste("r = ", r))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- rnorm(2926,mean = 0, sd = 1)
    #put in gaps
    X.t[c(100:500,1300:1400, 1700:1874, 2000:2200)] <- NA
    
    #calculate S.hat
    #MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
    ##use tapers to generate spectral estimate

    S.x.hat_MD <- rep(NA, times = floor(N/2) + 1)
    
    for(j in 0:floor(N/2)){
      k.vec <- rep(NA,times = K)
      for(k in 0:(K-1)){
        W.t <- eig_vecs[,k+1]*na.exclude(X.t)
        inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*na.omit(t.n)*j/N), na.rm = TRUE)
        k.vec[k + 1] <- abs(inner.sum)^2
      }
      S.x.hat_MD[j+1] <- mean(k.vec)
    }
    
    
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(S.x.hat_MD[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(S.x.hat_MD[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + S.x.hat_MD[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }
      else{
        bpvar.vec[i] <- 4*delta.f*(sum(S.x.hat_MD[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*S.x.hat_MD)
    
  }
  tmat_gps2[r,] <- trfunc.vec
  bmat_gps2[r,] <- bpvar.vec
}


###also Calculate AVAR#####

amat_gps <- oamat_gps <- matrix(NA, nrow = 300, ncol = 11)

for(i in 1:300){
  set.seed(i)
  #generate X.t
  #generate X.t
  X.t <- rnorm(2926,mean = 0, sd = 1)
  #put in gaps
  y <- X.t[-c(100:500,1300:1400, 1700:1874, 2000:2200)]
  avar.calc <- getAvars(N,y, taus = c(2^(0:9), floor(N/3)))
  amat_gps[i,] <- avar.calc$avarRes$avars
  oamat_gps[i,] <- avar.calc$avarRes$overavars
}



### Plots #####
#for plotting: bmat_gps, tmat_gps, amat_gps/oamat_gps
N <- 2048
taus <- c(2^(0:9),floor(N/3))

##tidy the data
tmat_gps2 %<>% t()
bmat_gps2 %<>% t()
dim(bmat_gps2)
three.mat <- rbind(oamat_gps, tmat_gps2)
method_labels <- rep(c("avar", "tr"), each = 300)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)



ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat_w,aes(tau,truth), size = 1.2)+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("Allan", "Transfer"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here 
  ### basically specifies the x and y position in terms of the plot size in unit scale. 
  theme_bw(base_size = 20)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  annotation_logticks()+
  #theme_bw()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  ggtitle("White Noise Comparison, Gaps")





#####################################
###Example 3: Flicker Process #######
#####################################
N <- 2048
trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat_flk <- bmat_flk <- matrix(NA, ncol = 300, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]

X.t_sims_flk <- matrix(NA, nrow = 300, ncol = 2048)
g = 0
for(k in c(2^(0:9), floor(N/3))){
  g = g + 1
  tau = k
  print(paste("flk g = ", g))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- TK95(N = 2048, alpha = 1)
    X.t_sims_flk[i,] <- X.t
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(MTSE_full$spectrum[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + MTSE_full$spectrum[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }else{
        bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat_flk[g,] <- trfunc.vec
  bmat_flk[g,] <- bpvar.vec
}


##AVAR calculation

amat_flk <- oamat_flk <- matrix(NA, nrow = 300, ncol = 11)

for(i in 1:300){
  #get X.t from simulation matrix
  X.t <- X.t_sims_flk[i,]
  
  avar.calc <- getAvars(N,X.t, taus = c(2^(0:9), floor(N/3)))
  amat_flk[i,] <- avar.calc$avarRes$avars
  oamat_flk[i,] <- avar.calc$avarRes$overavars
}


###### Plots #####
#for plotting: bmat_flk, tmat_flk, amat_flk/oamat_flk

taus <- c(2^(0:9),floor(N/3))

##tidy the data
tmat_flk %<>% t()
bmat_flk %<>% t()
dim(bmat_flk)
three.mat <- rbind(oamat_flk, tmat_flk,bmat_flk)
method_labels <- rep(c("avar", "tr", "bp"), each = 300)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)


sumDat=dat %>% group_by(method,tau) %>%
  summarise(median = median(log10(measurement)),
            lower25=quantile(log10(measurement),prob=.25),
            upper75=quantile(log10(measurement),prob=.75),
            min=min(log10(measurement)),
            max=max(log10(measurement)))

ggplot(data = sumDat,aes(x=log10(tau),y=median,col=method,fill=method))+
  ### width controls how wide these are, can change all below if they are too wide to see tau linear relationship well with real data
  geom_errorbar(aes(ymin=min,ymax=max),width=.1,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),width=.1,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),color="black",width=.1,position = "dodge")+
  ####
  ylab(expression(log[10](sigma(tau))))+
  xlab(expression(log[10](tau)))+
  ggtitle("Flicker Noise, No Gaps")+
  ### add true line below
  geom_abline(slope = 0,intercept = 3,size=1)
  
#### flicker with gaps ####


N <- 2048
trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat_flk_gps <- bmat_flk_gps <- matrix(NA, ncol = 300, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]

X.t_sims_flk_gps <- matrix(NA, nrow = 300, ncol = 2926)
g = 0
for(k in c(2^(0:9), floor(N/3))){
  g = g + 1
  tau = k
  print(paste("flk g = ", g))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- TK95(N = 2926, alpha = 1)
    X.t[c(100:500,1300:1400, 1700:1874, 2000:2200)] <- NA
    X.t_sims_flk_gps[i,] <- X.t
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(MTSE_full$spectrum[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + MTSE_full$spectrum[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }else{
        bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat_flk_gps[g,] <- trfunc.vec
  bmat_flk_gps[g,] <- bpvar.vec
}


##AVAR calculation

amat_flk_gps <- oamat_flk_gps <- matrix(NA, nrow = 300, ncol = 11)

for(i in 1:300){
  #get X.t from simulation matrix
  X.t <- na.omit(X.t_sims_flk_gps[i,])
  
  avar.calc <- getAvars(N,X.t, taus = c(2^(0:9), floor(N/3)))
  amat_flk_gps[i,] <- avar.calc$avarRes$avars
  oamat_flk_gps[i,] <- avar.calc$avarRes$overavars
}


##tidy the data
tmat_flk_gps %<>% t()
bmat_flk_gps %<>% t()
dim(bmat_flk_gps)
three.mat <- rbind(oamat_flk_gps, tmat_flk_gps,bmat_flk_gps)
method_labels <- rep(c("avar", "tr", "bp"), each = 300)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

sumDat=dat %>% group_by(method,tau) %>%
  summarise(median = median(log10(measurement)),
            lower25=quantile(log10(measurement),prob=.25),
            upper75=quantile(log10(measurement),prob=.75),
            min=min(log10(measurement)),
            max=max(log10(measurement)))

ggplot(data = sumDat,aes(x=log10(tau),y=median,col=method,fill=method))+
  ### width controls how wide these are, can change all below if they are too wide to see tau linear relationship well with real data
  geom_errorbar(aes(ymin=min,ymax=max),width=.1,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),width=.1,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),color="black",width=.1,position = "dodge")+
  ####
  ylab(expression(log[10](sigma(tau))))+
  xlab(expression(log[10](tau)))+
  ggtitle("Flicker Noise, Gaps")+
  ### add true line below
  geom_abline(slope = 0,intercept = 3.1,size=1)




####Flicker Noise with ARFIMA ######
N <- 2048
X.t_sims_flk2 <- matrix(NA, nrow = 300, ncol = N)
d <- 0.4999999999
for(i in 1:300){
  print(i)
  set.seed(i)
  X.t_sims_flk2[i,] <- arfima.sim(N,model = list(dfrac = d))
  
  #avar_saved_ARFIMA_49[i,] <- getAvars(N,y, taus)$avarRes$avars
}

trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat_flk2 <- bmat_flk2 <- matrix(NA, ncol = 300, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]

g = 0
for(k in c(2^(0:9), floor(N/3))){
  g = g + 1
  tau = k
  print(paste("flk g = ", g))
  
  for(i in 1:300){
    print(i)
    
    #generate X.t
    X.t <- X.t_sims_flk2[i,]
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(MTSE_full$spectrum[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + MTSE_full$spectrum[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }else{
        bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat_flk2[g,] <- trfunc.vec
  bmat_flk2[g,] <- bpvar.vec
}


###AVAR Calculation ####

amat_flk2 <- oamat_flk2 <- matrix(NA, nrow = 300, ncol = 11)

for(i in 1:300){
  #get X.t from simulation matrix
  X.t <- X.t_sims_flk2[i,]
  
  avar.calc <- getAvars(N,X.t, taus = c(2^(0:9), floor(N/3)))
  amat_flk2[i,] <- avar.calc$avarRes$avars
  oamat_flk2[i,] <- avar.calc$avarRes$overavars
}

##tidy the data
tmat_flk2 %<>% t()
bmat_flk2 %<>% t()
dim(bmat_flk_gps)
three.mat <- rbind(oamat_flk2[,-1], bmat_flk2[,-1])
method_labels <- rep(c("avar", "tr"), each = 300)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus[-1])


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)



breaks <- 10^(-10:10)
minor_breaks <- rep(1:9,21 )*rep(10^(-10:10), each = 9)

linedat=data.frame(tau=2:682)
linedat$truth= tavar_ARFIMA(682, d, sig.2.a = 1)
linedat$method=NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+ #width=length(unique(dat$tau))
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat,aes(tau,truth), size = 1.2)+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("Allan","Spectral"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here 
  ### basically specifies the x and y position in terms of the plot size in unit scale. 
  theme_bw(base_size = 20)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  annotation_logticks()+
  
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  ggtitle("Flicker Noise Comparison, No Gaps")


#### Flicker Noise with ARFIMA: Gaps ####

N <- 2048
X.t_sims_flk2_gps <- matrix(NA, nrow = 500, ncol = 2926)
d <- 0.4999999999

for(i in 1:500){
  print(i)
  set.seed(i)
  X.t <- arfima.sim(2926,model = list(dfrac = d))
  X.t[c(100:500,1300:1400, 1700:1874, 2000:2200)] <- NA
  X.t_sims_flk2_gps[i,] <- X.t
}

trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat_flk3_gps <- bmat_flk3_gps <- matrix(NA, ncol = 300, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]

g = 0
for(k in c(2^(0:9), floor(N/3))){
  g = g + 1
  tau = k
  print(paste("flk g = ", g))
  
  for(i in 1:300){
    print(i)
    
    #generate X.t
    X.t <- X.t_sims_flk2_gps[i,]
    
    #calculate S.hat
    #MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
    
    S.x.hat_MD <- rep(NA, times = floor(N/2) + 1)
    
    for(j in 0:floor(N/2)){
      k.vec <- rep(NA,times = K)
      for(k in 0:(K-1)){
        W.t <- eig_vecs[,k+1]*na.exclude(X.t)
        inner.sum <- sum(W.t*exp(-complex(real = 0, imaginary = 1)*2*pi*na.omit(t.n)*j/N), na.rm = TRUE)
        k.vec[k + 1] <- abs(inner.sum)^2
      }
      S.x.hat_MD[j+1] <- mean(k.vec)
    }
    
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(S.x.hat_MD[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(S.x.hat_MD[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + S.x.hat_MD[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }else{
        bpvar.vec[i] <- 4*delta.f*(sum(S.x.hat_MD[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*S.x.hat_MD)
    
  }
  tmat_flk3_gps[g,] <- trfunc.vec
  bmat_flk3_gps[g,] <- bpvar.vec
}


###AVAR Calculation ####

amat_flk2_gps <- oamat_flk2_gps <- matrix(NA, nrow = 500, ncol = 11)

for(i in 1:500){
  print(i)
  #get X.t from simulation matrix
  X.t <- X.t_sims_flk2_gps[i,]
  y <- na.exclude(X.t)
  avar.calc <- getAvars(N,y, taus = c(2^(0:9), floor(N/3)))
  amat_flk2_gps[i,] <- avar.calc$avarRes$avars
  oamat_flk2_gps[i,] <- avar.calc$avarRes$overavars
}

##tidy the data
tmat_flk3_gps %<>% t()
bmat_flk3_gps %<>% t()
dim(bmat_flk3_gps)
three.mat <- rbind(oamat_flk2_gps[1:300,-1],bmat_flk3_gps[,-1])
method_labels <- rep(c("avar", "tr"), each = 300)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus[-1])


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

sumDat=dat %>% group_by(method,tau) %>%
  summarise(median = median(measurement),
            lower25=quantile(measurement,prob=.25),
            upper75=quantile(measurement,prob=.75),
            min=min(measurement),
            max=max(measurement))


breaks <- 10^(-10:10)
minor_breaks <- rep(1:9,21 )*rep(10^(-10:10), each = 9)

linedat=data.frame(tau=1:682)
linedat$truth= c(NA,tavar_ARFIMA(682, d, sig.2.a = 1))
linedat$method=NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+ #width=length(unique(dat$tau))
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat,aes(tau,truth), size = 1.2)+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("Allan", "Spectral"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here 
  ### basically specifies the x and y position in terms of the plot size in unit scale. 
  theme_bw(base_size = 20)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  annotation_logticks()+
  
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  ggtitle("Flicker Noise Comparison, Gaps")


##maybe try random walk?


###################################
####  Example 4: Clock Data #######
###################################

#starting index of non-missing
clock_df <- dat180403$SrYb
clock_df <- clock_df[min(which(!is.na(clock_df))):length(clock_df)]

###split function to get runs in data between missing values
t.n <- 1:length(clock_df)
t.n <- t.n[-c(which(is.na(clock_df)))]

splits <- split(t.n, cumsum(c(1,diff(t.n) != 1)))
length(splits)
splits[1:10]


count = 0
for(i in 1:length(splits)){
  if(length(splits[[i]])< 10){
    count = count + 1
  }
}
count

#get indices of elements that are part of runs < 10 in length
ind_to_rmv <- unlist(splits[which(lengths(splits)<10)], use.names = FALSE)

clock_df_omitted <- clock_df[-ind_to_rmv]
length(clock_df) - length(clock_df_omitted)
t.n <- 1:length(clock_df_omitted)
t.n <- t.n[-c(which(is.na(clock_df_omitted)))]



clock_MTSE_omitted <- multitaper_est(clock_df_omitted, NW = 2, K = 3)
f <- seq(0, 0.5, length.out = length(clock_MTSE_omitted$spectrum))
f.all <- seq(0, 0.5, length.out = length(clock_MTSE$spectrum))

plot(log10(f),log10(clock_MTSE_omitted$spectrum))
plot(log10(f.all),log10(clock_MTSE$spectrum), type = "l", ylab = "log10(S(f))", xlab = "log10(f)")
abline(v = -1.4, lwd = 2, col = "red")
line1.f <- seq(-1.4, log10(0.5), length.out = 100)
a1 = -1.59
b1 = -0.811
lines(line1.f, a1 + b1*line1.f, col = "blue", lwd = 2)
a2 = -0.04
b2 = 0.3
line2.f <- seq(-4.527321, -1.4, length.out = 100)
lines(line2.f, a2 + b2*line2.f, col = "blue", lwd = 2)

freqs <- log10(f.all[which(log10(f.all)>-1.4)])
spec.hat <- log10(clock_MTSE$spectrum[which(log10(f.all)> -1.4)])
lin.fit <- lm(spec.hat ~ freqs)

freqs1 <- log10(f.all[which(log10(f.all)< -1.4)])[-1]
spec.hat1 <- log10(clock_MTSE$spectrum[which(log10(f.all)< -1.4)])[-1]
lin.fit1 <- lm(spec.hat1 ~ freqs1)

summary(lin.fit1)
abline(a = -0.04, b = 0.3, col = "blue")
plot(lin.fit)

library(car)
durbinWatsonTest(lin.fit, method = "normal")




N <- (length(clock_MTSE$spectrum) - 1)*2
t.vec_clock_omitted <- bpvar.vec_clock_omitted <- rep(NA, times = 15)

r = 0
s.hat <- clock_MTSE_omitted$spectrum
delta.f <- f[2]
for(k in c(2^(0:13), floor(N/3))){
  r = r + 1
  tau = k
  
  #calculate bandpass variance
  
  if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
    f.min.index <- which(f == 1/(4*tau))
    f.max.index <- which(f == 1/(2*tau))
    bpvar.vec_clock[r] <- 4*delta.f*(sum(clock_MTSE$spectrum[f.min.index:(f.max.index - 1)]))
  } else{
    f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
    f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
    if(f.min.index == f.max.index){
      bpvar.vec_clock_omitted[r] <- 4*(s.hat[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + s.hat[f.min.index]*(1/(2*tau) - f[f.min.index]))
    }else{
      bpvar.vec_clock_omitted[r] <- 4*delta.f*(sum(s.hat[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
    }
  }
  
  #calculate transfer function AVAR
  G.vec <- transfer.func(f, tau)
  G.vec[1] <- 1
  t.vec_clock_omitted[r] <- f[2]*sum(G.vec*s.hat)
  
}

taus <-  c(2^(0:13), floor(N/3))
avar_clock <- getAvars(length(na.omit(clock_df)), na.omit(clock_df), taus)

par(mfrow = c(1,1), pch = 19, cex = 1.5)
plot(log10(taus), log10(t.vec_clock), ylim = c(-6,-1), main = "", ylab = expression(log10(sigma^2*(tau))), xlab = expression(log10(tau)))
#points( log10(taus), log10(bpvar.vec_clock), col = "red")
points( log10(taus), log10(avar_clock$avarRes$overavars), col = "red")
legend(x = 3, y = -1.3, legend = c("Spectral", "Allan"), col = c("black", "red"), pch = 19)

plot( log10(taus), log10(bpvar.vec_clock_omitted), main = "Clock Data Comparison (<10 Length Gaps Removed)", ylab = expression(log10(sigma(tau))), xlab = expression(log10(tau)), col = "green")
points( log10(taus), log10(t.vec_clock_omitted), col = "orange")
points( log10(taus), log10(avar_clock$avarRes$overavars), col = "blue")
legend(x = 3.5, y = -2, legend = c("bp", "tr", "oavar"), col = c("green", "orange", "blue"), pch = 19)


#### Confidence intervals ####

### AVAR:
N <- (length(clock_MTSE$spectrum) - 1)*2
taus <-  c(2^(0:13), floor(N/3))
edf <- rep(NA, times = length(taus))

for(i in 1:15){
  m = taus[i]
  edf[i] <- ((3*(N-1)/(2*m)) - (2*(N-2)/N))*(4*m^2)/(4*m^2 + 5)
}
length(edf)
pr = 0.025
avar.clock.CI <- data.frame(CI.upper = avar_clock$avarRes$overavars*edf/qchisq(1-pr, edf), 
                            CI.lower = avar_clock$avarRes$overavars*edf/qchisq(pr, edf))

points(log10(taus),log10(avar.clock.CI$CI.upper), pch = "-", col = "grey40")
points(log10(taus),log10(avar.clock.CI$CI.lower), pch = "-", col = "grey40")

plot(clock_df, type = "l", ylab = "Clock Noise Data", xlab = "Time [s]")











####Functions #####

getAvars=function(N, y,taus){
  
  avars=numeric(length(taus))
  overlapping_avars= numeric(length(taus))
  
  for (i in 1:length(taus)){
    avars[i]=avar_fn(y,taus[i])
    overlapping_avars[i] = overlapping_avar_fn(y,taus[i])
  }
  
  m1=data.frame(taus=taus,avars=avars)  
  fit=lm(log(sqrt(avars))~log(taus),data = m1)
  slope=as.numeric(fit$coefficients[2])
  int=as.numeric(fit$coefficients[1])
  
  avarRes=data.frame(taus=taus,avars=avars, overavars=overlapping_avars,N=N,slope=slope,int=int)  
  
  ##########################################
  # get SE
  ##########################################
  m2=data.frame(taus=taus,oavars=overlapping_avars)  
  fit2=lm(log(sqrt(oavars))~log(taus),data = m2)
  slope2=as.numeric(fit2$coefficients[2])
  int2=as.numeric(fit2$coefficients[1])
  
  SEests=data.frame()
  
  onew=data.frame(N=N,out=exp(int2+slope2*log(N)), type="OAD")
  new=data.frame(N=N,out=exp(int+slope*log(N)), type="AD")
  new2=data.frame(N=N,out=sd(y)/sqrt(N), type="SE")
  new3=data.frame(N=N,out=1/sqrt(N), type="true")
  SEests=bind_rows(new,new2,new3,onew)
  
  return(list(avarRes=avarRes,SEests=SEests))
}





############ scraps ###################

###################################
###Example 2: AR(1) Process #######
###################################

N <- 2048
trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat_AR <- bmat_AR <- matrix(NA, ncol = 300, nrow = 11)
X.t_sims_AR <- matrix(NA, nrow = 300, ncol = 2048)


f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]
phi <- 0.5
r <- 0

for(k in c(2^(0:9), floor(N/3))){
  r = r + 1
  tau = k
  print(paste("AR r = ", r))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- arima.sim(n = N, list(ar = c(phi)))
    X.t_sims_AR[i,] <- X.t
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(MTSE_full$spectrum[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + MTSE_full$spectrum[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }
      else{
        bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat_AR[r,] <- trfunc.vec
  bmat_AR[r,] <- bpvar.vec
}


###also Calculate AVAR#####

amat_AR <- oamat_AR <- matrix(NA, nrow = 300, ncol = 11)

for(i in 1:300){
  set.seed(i)
  #generate X.t
  X.t <- X.t_sims_AR[i,]
  
  avar.calc <- getAvars(N,X.t, taus = c(2^(0:9), floor(N/3)))
  amat_AR[i,] <- avar.calc$avarRes$avars
  oamat_AR[i,] <- avar.calc$avarRes$overavars
}




###### Plots #####
#for plotting: bmat_AR, tmat_AR, amat_AR/oamat_AR

taus <- c(2^(0:9),floor(N/3))

##tidy the data
tmat_AR %<>% t()
bmat_AR %<>% t()
dim(bmat_AR)
three.mat <- rbind(oamat_AR, tmat_AR,bmat_AR)
method_labels <- rep(c("avar", "tr", "bp"), each = 300)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

sumDat=dat %>% group_by(method,tau) %>%
  summarise(median = median(measurement),
            lower25=quantile(measurement,prob=.25),
            upper75=quantile(measurement,prob=.75),
            min=min(measurement),
            max=max(measurement))

ggplot(data = sumDat,aes(x=tau,y=median,col=method,fill=method))+
  ### width controls how wide these are, can change all below if they are too wide to see tau linear relationship well with real data
  geom_errorbar(aes(ymin=min,ymax=max),width=.1,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),width=.1,position = "dodge")+
  geom_crossbar(aes(ymin=lower25,ymax=upper75),color="black",width=.1,position = "dodge")+
  ####
  scale_y_log10()#+
#scale_x_log10()+
ylab(expression(log[10](sigma(tau))))+
  xlab(expression(log[10](tau)))+
  ggtitle("AR(1) Comparison, No Gaps")#+
### add true line below
#geom_abline(slope = -1,intercept = 0,size=1)





##### longer white noise #####
N <- 2^12
trfunc.vec <- bpvar.vec <- rep(NA, times = 500)

tmat_lng <- bmat_lng <- matrix(NA, ncol = 500, nrow = 12)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]


r = 0
for(k in c(2^(0:10), floor(N/3))){
  r = r + 1
  tau = k
  print(paste("r = ", r))
  
  for(i in 1:500){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- rnorm(N,mean = 0, sd = 1)
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, NW = 3, K = 5)
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(MTSE_full$spectrum[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + MTSE_full$spectrum[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }
      else{
        bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat_lng[r,] <- trfunc.vec
  bmat_lng[r,] <- bpvar.vec
}




###also Calculate AVAR LONG#####

amat_lng <- oamat_lng <- matrix(NA, nrow = 500, ncol = 12)

for(i in 1:300){
  set.seed(i)
  #generate X.t
  X.t <- rnorm(N,mean = 0, sd = 1)
  
  avar.calc <- getAvars(N,X.t, taus = c(2^(0:10), floor(N/3)))
  amat_lng[i,] <- avar.calc$avarRes$avars
  oamat_lng[i,] <- avar.calc$avarRes$overavars
}


###### White Noise with Gaps #####
##### With gaps
N <- 2^12
trfunc.vec <- bpvar.vec <- rep(NA, times = 500)

tmat_gps_lng <- bmat_gps_lng <- matrix(NA, ncol = 500, nrow = 12)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]


r = 0
for(k in c(2^(0:10), floor(N/3))){
  r = r + 1
  tau = k
  print(paste("r = ", r))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- rnorm(5851,mean = 0, sd = 1)
    #put in gaps
    X.t[c(100:500,1300:1800, 1700:1874, 2000:2200,3140:3600)] <- NA
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, NW = 2, K = 3)
    
    #calculate bandpass variance
    
    if(sum(f-1/(4*tau) == 0) & sum(f-1/(2*tau) == 0)){
      f.min.index <- which(f == 1/(4*tau))
      f.max.index <- which(f == 1/(2*tau))
      bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:(f.max.index - 1)]))
    }
    else{
      f.min.index <- min(which(f>1/(4*tau) & f<1/(2*tau)))
      f.max.index <- max(which(f>1/(4*tau) & f<1/(2*tau)))
      if(f.min.index == f.max.index){
        bpvar.vec[i] <- 4*(MTSE_full$spectrum[f.min.index-1]*(f[f.min.index] - 1/(4*tau)) + MTSE_full$spectrum[f.min.index]*(1/(2*tau) - f[f.min.index]))
      }
      else{
        bpvar.vec[i] <- 4*delta.f*(sum(MTSE_full$spectrum[f.min.index:f.max.index]) + (f[f.min.index] - 1/(4*tau)) + (1/(2*tau) - f[f.max.index]) )
      }
    }
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 1
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full$spectrum)
    
  }
  tmat_gps_lng[r,] <- trfunc.vec
  bmat_gps_lng[r,] <- bpvar.vec
}


###also Calculate AVAR#####

amat_gps_lng <- oamat_gps_lng <- matrix(NA, nrow = 500, ncol = 12)

for(i in 1:300){
  set.seed(i)
  #generate X.t
  #generate X.t
  X.t <- rnorm(5851,mean = 0, sd = 1)
  #put in gaps
  y <- X.t[-c(100:500,1300:1400, 1700:1874, 2000:2200)]
  avar.calc <- getAvars(N,y, taus = c(2^(0:9), floor(N/3)))
  amat_gps[i,] <- avar.calc$avarRes$avars
  oamat_gps[i,] <- avar.calc$avarRes$overavars
}



taus <- c(2^(0:10), floor(N/3))
##tidy the data
tmat_lng %<>% t()
bmat_lng %<>% t()
dim(bmat_lng)
three.mat <- rbind(oamat_lng, tmat_lng,bmat_lng)
method_labels <- rep(c("avar", "tr", "bp"), each = 500)
df.messy <- as.data.frame(cbind(method_labels, three.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

sumDat=dat %>% group_by(method,tau) %>%
  summarise(median = median(measurement),
            lower25=quantile(measurement,prob=.25),
            upper75=quantile(measurement,prob=.75),
            min=min(measurement),
            max=max(measurement))

linedat=data.frame(tau=1:1000)
linedat$truth=1/linedat$tau#+linedat$tau
linedat$method=NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat,aes(tau,truth))+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("Allan","Bandpass","Transfer Function"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here 
  ### basically specifies the x and y position in terms of the plot size in unit scale. 
  theme_bw()+
  theme(legend.position = c(.15, .2))+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  annotation_logticks()+
  #theme_bw()+
  ylab(expression(sigma(tau)))+
  xlab(expression(tau))+
  ggtitle("White Noise Comparison, Gaps")

