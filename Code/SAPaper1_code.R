############ SA Paper 1 #################
## This is code for simulation studies ##
## and clock data analysis for Paper 1 ##
#########################################

#########################################
#######           Libraries         #####
#########################################

library(tidyverse)
library(RSpectra) #eigensolving
library(fields) #dist.mat
library(RobPer) #flicker noise

###########################
#### things we'll need ####
###########################

## i
im <- complex(real = 0, imaginary = 1)

## MTSE for gappy data function 
multitaper_est <- function(X.t, W, K){
  X.t <- X.t - mean(X.t, na.rm = TRUE) #demean
  N.long <- length(X.t)
  t.n <- 1:N.long
  missing.indices <- which(is.na(X.t))
  t.n[which(is.na(X.t))] <- NA
  
  
  dist.mat <- rdist(na.omit(t.n))
  
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  print("A matrix computed")
  
  eigdec <- eigs(A.prime, k = K, which = "LM")
  eig_vecs <- eigdec$vectors #get only the vectors
  print("tapers computed")
  
  if(K ==1){
    if (mean(Re(eig_vecs))<0){
      eig_vecs <- -eig_vecs
    }
  }
  
  if(K == 2  || K == 3){
    
    if (mean(Re(eig_vecs[,1]))<0){
      eig_vecs[,1] <- -eig_vecs[,1]
    }
    if (Re(eig_vecs[2,2] - eig_vecs[1,2])<0){
      eig_vecs[,2] <- -eig_vecs[,2]
    }
    
    if(K == 3){
      if (mean(Re(eig_vecs[,3]))<0){
        eig_vecs[,3] <- -eig_vecs[,3]
      }
    }
  }
  if(K >=4){
    #some sign maintenance
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
  }
  
  print("sign maintenance done")
  
  ##use tapers to generate spectral estimate
  N <- length(na.omit(t.n))
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
  
  
  return(list("tapers" = eig_vecs, "e.values" = eigdec$values, "spectrum" = S.x.hat_MD))
}


######################################################################################################
###### Calculate sigma.hat_tau using transfer function method (see "reappraisal...", page 70) ########
######################################################################################################

#G_tau(f) function
transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}

###########   Allan Variance Calculation    #############
######## (eq'n at bottom of p. 70 of reappraisal) #######
#input: spectral estimate (as a vector), taus (as a vector) where you want the AVAR calculated
#output: a vector of the AVAR estimates

AVAR_trfunc <- function(spectral_est, taus){
  out <- rep(NA, times = length(taus))
  f <- seq(0,0.5,length.out = length(spectral_est))
  
  for(i in 1:length(taus)){
    G.vec <- transfer.func(f, taus[i])
    G.vec[1] <- 1
    out[i] <- f[2]*sum(G.vec*spectral_est)
  }
  
  return(out)
}




######## Simulation Studies ##########

######################################
###### Study 1: White Noise ##########
######   WN(0,1), no gaps    #########
######################################
N <- 7200
trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat <- bmat <- matrix(NA, ncol = 300, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]


r = 0
for(k in c(2^(0:9), floor(N/3))){
  r = r + 1
  tau = k
  print(paste("r = ", r))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- rnorm(N,mean = 0, sd = 1)
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, W = 0.00097, K = 5)
    
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


boxplot(tmat[6,], bmat[6,], oamat[,6])


 ###also Calculate AVAR###

amat <- oamat <- matrix(NA, nrow = 300, ncol = 11)

for(i in 1:300){
  set.seed(i)
  print(i)
  #generate X.t
  X.t <- rnorm(N,mean = 0, sd = 1)
  
  avar.calc <- getAvars(N,X.t, taus = c(2^(0:9), floor(N/3)))
  amat[i,] <- avar.calc$avarRes$avars
  oamat[i,] <- avar.calc$avarRes$overavars
}

#for plotting: bmat, tmat, amat/oamat


### Plots ###
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



######################################
###### Study 2: Flicker Noise ########
######      with  gaps     ###########
######################################

N <- 4096
trfunc.vec <- bpvar.vec <- rep(NA, times = 300)

tmat_flk_gps <- bmat_flk_gps <- matrix(NA, ncol = 300, nrow = 11)

f <- seq(0,0.5,length.out = N/2 + 1) #grid of frequencies
delta.f <- f[2]

X.t_sims_flk_gps <- matrix(NA, nrow = 300, ncol = (4096 + 1128)) #add in extras for the NAs we will be creating
g = 0
for(k in c(2^(0:9), floor(N/3))){
  g = g + 1
  tau = k
  print(paste("flk g = ", g))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    #generate X.t
    X.t <- TK95(N = 5224, alpha = 1)
    X.t[c(100:500,1300:1400, 1700:1874, 3000:3450)] <- NA
    X.t_sims_flk_gps[i,] <- X.t
    
    #calculate S.hat
    MTSE_full <- multitaper_est(X.t, W = 0.005, K = 5)
    
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




MTSE_saved <- list()

for(i in 1:300){
  print(i)
  set.seed(i)
  #generate X.t
  X.t <- TK95(N = 5224, alpha = 1)
  X.t[c(100:500,1300:1400, 1700:1874, 3000:3450)] <- NA
  X.t_sims_flk_gps[i,] <- X.t
  
  #calculate S.hat
  MTSE_full <- multitaper_est(X.t, W = 0.009, K = 5)
  
  MTSE_saved[[i]] <- MTSE_full
}


