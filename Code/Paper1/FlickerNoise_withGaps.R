### read in the file with functions
# source("") 

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


