rm(list=ls())
source("/home/aak3/NIST/ClockDataAnalysis/Code/SA_ImportantFunctions.R")
source("Code/SA_ImportantFunctions.R")

######################################################################
### simulate the data
### need x.t, t.vec, and N for different dataset
######################################################################
sizeOfData=1000
t.vec <- 1:sizeOfData #time vector

# omitted<-c(20:35,50:63, 100:130)
omitted<-c(200:350,500:580, 700:850)

t.vec[omitted] <- NA #take out values
t.vec <- na.omit(t.vec) #vector of times with data
x.t <- rnorm(sizeOfData) #data
x.t[omitted] <- NA #take out values
x.t=na.omit(x.t)

#########################

N=length(x.t)

######################################################################
######################################################################


######################################################################
### calculate the MT spectral est and our avar
######################################################################

test=spectralEstWithUnc(x.t = x.t,t.vec = t.vec,numTapers = 8)

taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]

ours=AVAR_trfunc_withUnc(spectral_est = test$spec.hat,taus = taus,Cov.mat_chave = test$Cov.mat)

######################################################################
### calculate the avar the old way
######################################################################

oldAvars=getAvars(N,y = x.t,taus = taus)
oldAvarsUnc=avar_CI(CI.level = .68,noise_type = "white noise",avar_type = "NA",
                    avars = oldAvars$avarRes$overavars,
                    taus = taus,
                    N=N)

######################################################################
### make data frame of results
######################################################################

dfres=data.frame(tau=taus,
                 avar=ours$avar,
                 lower=ours$avar-sqrt(ours$avarVar),
                 upper=ours$avar+sqrt(ours$avarVar),
                 type="spectral")
olddfres=data.frame(tau=taus,
                    avar=oldAvars$avarRes$overavars,
                    lower=oldAvarsUnc$lower,
                    upper=oldAvarsUnc$upper,
                    type="current")
dfres=bind_rows(dfres,olddfres)

ggplot(dfres,aes(tau,avar,ymin=lower,ymax=upper,color=type))+
  geom_point()+
  geom_errorbar()+
  ### add true straight line below
  geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))
