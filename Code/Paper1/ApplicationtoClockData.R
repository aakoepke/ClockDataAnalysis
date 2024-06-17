################## Spectral-Based Method for AVAR estimate ################
source("Code/SA_ImportantFunctions.R")

### get clock data ###
dat180403 <- read_delim("//cfs2w.nist.gov/unix$/776unix/cmb15/ClockDataAnalysis/Data/180403 optical analysis TiS_fixed.dat",
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE, skip = 2)
colnames(dat180403)=c("MJD","AlYb","SrYb","AlSr")

plot(dat180306$AlYb[min(which(!is.na(dat180306$AlYb))):length(dat180306$AlYb)])

df <- dat180306$AlYb[min(which(!is.na(dat180306$AlYb))):length(dat180306$AlYb)]





####EDA on missing points
plot(df)
sum(is.na(df))

missing.indices <- which(is.na(df))
missing.indices[1:20]
differenced <- diff(missing.indices)
differenced[1:20]
sum(differenced > 1)
jumps <- which(differenced>1)
differenced[1220:1230]

missing.indices[1:10]
differenced[1:10]

new.df <- df
new.df[is.na(df)] <- Inf
runs <- rle(new.df)
runs
NA_period_length <- runs$lengths[runs$values == Inf]

NA_period_length


remove.indices <- c()
start <- 0
for(i in 1:length(NA_period_length)){
  if(NA_period_length[i]< 100){
    remove.indices <- append(remove.indices,missing.indices[(1 + start):sum(NA_period_length[1:i])])
    start <- sum(NA_period_length[1:i])
  }
  if(NA_period_length[i]>100){
    start <- sum(NA_period_length[1:i])
  }
  
}

#take out small gaps
new.df <- df[-c(remove.indices)]

#look at gaps now
new.df2 <- new.df
new.df2[is.na(new.df)] <- Inf
runs2 <- rle(new.df2)
runs2
NA_period_length2 <- runs2$lengths[runs2$values == Inf]

NA_period_length2

#only long gaps left

###EDA to check for white noise

hist(new.df, breaks = 600, 
     main = "", 
     xlab = "Fractional Frequency Deviates")

#sample acf
acf(new.df[1:2000], lag.max = 20)
acf(na.omit(new.df), lag.max = 20, main = "")



### calculate tapers ###
df.short <- new.df[22000:28000]
t.n <- 1:length(df.short)
t.n[which(is.na(df.short))] <- NA
N <- length(na.omit(t.n))
N
plot(df.short, pch = 19, ylab = "Clock Ratio Deviates", xlab = "t")

sum(is.na(t.n))
sum(is.na(df.short))
taperMatrix <- get_tapers(t.n, W = 12/N, K = 8)
taperMatrix$e.values
tapMat <- taperMatrix$tapers

### calculate MTSE ###
MTSE <- MT_spectralEstimate(df.short, tapMat)
plot(log(MTSE$freqs[-1]), log(MTSE$spectrum[-1]))
abline(v = -7)
abline(v = -6) 

### calculate AVAR ###
taus <- c(2^(0:10))  #for full data set: c(2^(0:12), 5000,6000, 9000, 10000,15000)
trfunc.vec <- rep(NA, times = length(taus))

#for(i in 1:length(taus)){
for(i in 1:length(taus)){ 
 tau = taus[i]

    #calculate bandpass variance
    # temp_bp <- integrate(approxfun(f, MTSE_mat[i,]), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 1000)
    # # temp_bp_LS <- integrate(approxfun(f[-1], lsperio), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 1000)
    # bpvar.vec[i] <- 4*temp_bp$value
    #bpvar.vec_LS[i] <- 4*temp_bp_LS$value
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(MTSE$freqs, tau)
    G.vec[1] <- 0
    trfunc.vec[i] <- MTSE$freqs[2]*sum(G.vec*MTSE$spectrum)
    #trfunc.vec_LS[i] <- f[2]*sum(G.vec[-1]*lsperio)
}



###also Calculate AVAR###
avar.calc <- getAvars(N,na.omit(df.short), taus = taus)
oavar <- avar.calc$avarRes$overavars

plot(log(taus), log(oavar), ylim = c(-10,-2), xlim = c(0,12), ylab = "AVAR Estimate")
points(log(taus), log(trfunc.vec), pch = 19, col = "blue")
legend(x = 1, y = -8, pch = c(1,19), legend = c("OAVAR", "Spectral"), col = c("black", "red"))
abline(a = 0, b = -1)


##save these to an .Rdata file so they don't get lost
save(oavar, taperMatrix, MTSE, Cov.mat_chave, file = "Application_Results_2.28.RData")
load("Application_Results_2.28.RData")

### error bars ###

  N.fourier <- floor(N/2) + 1
  freq <- seq(0,0.5, length.out = N.fourier)
  numTapers = 8
  delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
  t.vec <- na.omit(t.n)
  ### calculate the covariance matrix 
  Cov.mat_chave <- matrix(NA, nrow = N.fourier, ncol = N.fourier)
  
  for(i in 1:N.fourier){
    j = 1
    print(i)
    while(j <= i){
      Cov.mat_chave[i,j] <- norm(Conj(t(tapMat*exp(-im*2*pi*freq[i]*t.vec)*(1/sqrt(numTapers))))%*%(tapMat*exp(-im*2*pi*freq[j]*t.vec)*(1/sqrt(numTapers))), type = "2") 
      j = j+1
    }
  }
  
  Cov.mat_chave[upper.tri(Cov.mat_chave)] <- t(Cov.mat_chave)[upper.tri(Cov.mat_chave)]

  ours=AVAR_trfunc_withUnc(spectral_est = MTSE$spectrum,taus = taus,Cov.mat_chave = Cov.mat_chave)
  
######################################################################
### calculate the avar the old way
######################################################################

oldAvars=getAvars(N,y = na.omit(df.short),taus = taus)
oldAvarsUnc=avar_CI(CI.level = .68,noise_type = "white noise",avar_type = "NA",
                    avars = oldAvars$avarRes$overavars,
                    taus = taus,
                    N=N)

######################################################################
### make data frame of results
######################################################################

dfres=data.frame(tau=taus,
                 avar=ours$avar,
                 lower = ours$avar -sqrt(ours$avarVar),
                 upper=ours$avar +sqrt(ours$avarVar),
                 type="spectral")
olddfres=data.frame(tau=taus,
                    avar= oldAvars$avarRes$overavars,
                    lower=oldAvarsUnc$lower,
                    upper=oldAvarsUnc$upper,
                    type="OAVAR")
dfres=bind_rows(dfres,olddfres)

ggplot(dfres,aes(tau,avar,ymin=lower,ymax=upper,color=type))+
  geom_point()+
  geom_errorbar()+
  ### add true straight line below
  #geom_abline(slope = -1,intercept = 0,size=1)+
  scale_color_discrete(labels= c(expression(paste(hat(sigma)[AVAR]^2)), expression(paste(hat(sigma)[SPEC]^2))),name="Estimate")+
  theme_bw(base_size = 20)+
  theme(legend.text.align = 0, legend.key.height = unit(1.5, "cm"))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))






spectralEstWithUnc=function(x.t, tapers){
  N <- dim(tapers)[1]
  
  taper.mat = tapers*(1/dim(tapers)[2]) #with scaling
 
  
  N.fourier <- floor(N/2) + 1
  freq <- seq(0,0.5, length.out = N.fourier)

  delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
  
  ### calculate the covariance matrix 
  Cov.mat_chave <- matrix(NA, nrow = N.fourier, ncol = N.fourier)
  
  for(i in 1:N.fourier){
    j = 1
    print(i)
    while(j <= i){
      Cov.mat_chave[i,j] <- norm(Conj(t(taper.mat*exp(-1i*2*pi*freq[i]*t.vec))%*%(V.mat$tapers*exp(-1i*2*pi*freq[j]*t.vec)), type = "2")) 
      j = j+1
    }
  }
  #scale
  Cov.mat_chave = Cov.mat_chave*(1/8)
  Cov.mat_chave[upper.tri(Cov.mat_chave)] <- t(Cov.mat_chave)[upper.tri(Cov.mat_chave)]
  
  
  return(list(freq=freq,
              spec.hat=MTSE_full$spectrum,Cov.mat=Cov.mat_chave))
  
}

