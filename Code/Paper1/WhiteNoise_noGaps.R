### read in the file with functions
# source("")

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