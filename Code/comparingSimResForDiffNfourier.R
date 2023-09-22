
resName="simDat_200_Sep22"
res200=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
res200$V.mat$e.values

resName="simDat_500_Sep22"
res500=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
res500$V.mat$e.values

resName="simDat_1000_Sep22"
res1000=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
res1000$V.mat$e.values

resName="simDat_1000_Sep22biggerW"
res1000biggerW=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
res1000biggerW$V.mat$e.values

resName="simDat_2000_Sep22"
res2000=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
res2000$V.mat$e.values

resName="simDat_2000_Sep22biggerW"
res2000biggerW=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
res2000biggerW$V.mat$e.values

dat1=data.frame(ratio="AlYb",
                n.fourier=length(res200$MTSE_full$freqs),
                W=res200$W,
                freq=res200$MTSE_full$freqs,
                spectrum=res200$MTSE_full$spectrum)
dat2=data.frame(ratio="AlYb",
                n.fourier=length(res500$MTSE_full$freqs),
                W=res500$W,
                freq=res500$MTSE_full$freqs,
                spectrum=res500$MTSE_full$spectrum)
dat3=data.frame(ratio="AlYb",
                n.fourier=length(res1000$MTSE_full$freqs),
                W=res1000$W,
                freq=res1000$MTSE_full$freqs,
                spectrum=res1000$MTSE_full$spectrum)
dat4=data.frame(ratio="AlYb",
                n.fourier=length(res1000biggerW$MTSE_full$freqs),
                W=res1000biggerW$W,
                freq=res1000biggerW$MTSE_full$freqs,
                spectrum=res1000biggerW$MTSE_full$spectrum)
dat5=data.frame(ratio="AlYb",
                n.fourier=length(res2000$MTSE_full$freqs),
                W=res2000$W,
                freq=res2000$MTSE_full$freqs,
                spectrum=res2000$MTSE_full$spectrum)
dat6=data.frame(ratio="AlYb",
                n.fourier=length(res2000biggerW$MTSE_full$freqs),
                W=res2000biggerW$W,
                freq=res2000biggerW$MTSE_full$freqs,
                spectrum=res2000biggerW$MTSE_full$spectrum)
allSpec=bind_rows(dat1,dat2,dat3,dat4,dat5,dat6)


# ggplot(allSpec,aes(freq,spectrum,col=factor(n.fourier)))+
#   geom_point(alpha = 0.5)
ggplot(allSpec,aes(freq,spectrum,col=factor(interaction(n.fourier,W))))+
  # geom_line()+
  # facet_wrap(~n.fourier)+
  geom_smooth()
  # scale_y_log10()+
  # scale_x_log10()+
  # geom_hline(yintercept = log10(1))
ggplot(allSpec,aes(freq,spectrum))+
  geom_line()+
  facet_wrap(~n.fourier+W)+
  geom_smooth()+
  scale_y_log10()+
  scale_x_log10()
  
dat1=data.frame(ratio="AlYb",
                n.fourier=length(res200$MTSE_full$freqs),
                W=res200$W,
                tau=res200$avarOut$tau,
                avar=res200$avarOut$avar,
                var=res200$avarOut$var)
dat2=data.frame(ratio="AlYb",
                n.fourier=length(res500$MTSE_full$freqs),
                W=res500$W,
                tau=res500$avarOut$tau,
                avar=res500$avarOut$avar,
                var=res500$avarOut$var)
dat3=data.frame(ratio="AlYb",
                n.fourier=length(res1000$MTSE_full$freqs),
                W=res1000$W,
                tau=res1000$avarOut$tau,
                avar=res1000$avarOut$avar,
                var=res1000$avarOut$var)
dat4=data.frame(ratio="AlYb",
                n.fourier=length(res2000$MTSE_full$freqs),
                W=res2000$W,
                tau=res2000$avarOut$tau,
                avar=res2000$avarOut$avar,
                var=res2000$avarOut$var)

dat5=data.frame(ratio="AlYb",
                n.fourier=length(res1000biggerW$MTSE_full$freqs),
                W=res1000biggerW$W,
                tau=res1000biggerW$avarOut$tau,
                avar=res1000biggerW$avarOut$avar,
                var=res1000biggerW$avarOut$var)
dat6=data.frame(ratio="AlYb",
                n.fourier=length(res2000biggerW$MTSE_full$freqs),
                W=res2000biggerW$W,
                tau=res2000biggerW$avarOut$tau,
                avar=res2000biggerW$avarOut$avar,
                var=res2000biggerW$avarOut$var)

allavar=bind_rows(dat1,dat2,dat3,dat4,dat5)




### plot avar
ggplot(allavar,aes(tau,avar,ymin=avar-var,ymax=avar+var,col=factor(n.fourier)))+
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

ggplot(allavar,aes(tau,avar,ymin=avar-var,ymax=avar+var,col=factor(interaction(n.fourier,W))))+
  geom_point()+
  geom_errorbar()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  facet_wrap(~tau,scales = "free")







### in this group, didn't account for bw
# resName="AlYb040318res"
# 
# resBIG=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# length(resBIG$freq)
# 
# resName="AlYb040318_100"
# res100=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# resName="AlYb040318_500"
# res500=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# 
# resName="AlYb040318_1000"
# res1000=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/Data/resultsFor",resName,".Rds",sep=""))
# dat1=data.frame(ratio="AlYb",
#                 n.fourier=length(resBIG$MTSE_full$freqs),
#                 freq=resBIG$MTSE_full$freqs,
#                 spectrum=resBIG$MTSE_full$spectrum)
# dat2=data.frame(ratio="AlYb",
#                 n.fourier=length(res100$MTSE_full$freqs),
#                 freq=res100$MTSE_full$freqs,
#                 spectrum=res100$MTSE_full$spectrum)
# dat3=data.frame(ratio="AlYb",
#                 n.fourier=length(res500$MTSE_full$freqs),
#                 freq=res500$MTSE_full$freqs,
#                 spectrum=res500$MTSE_full$spectrum)
# dat4=data.frame(ratio="AlYb",
#                 n.fourier=length(res1000$MTSE_full$freqs),
#                 freq=res1000$MTSE_full$freqs,
#                 spectrum=res1000$MTSE_full$spectrum)
# allSpec=bind_rows(dat1,dat2,dat3,dat4)
# dat1=data.frame(ratio="AlYb",
#                 n.fourier=length(resBIG$MTSE_full$freqs),
#                 tau=resBIG$avarOut$tau,
#                 avar=resBIG$avarOut$avar,
#                 var=resBIG$avarOut$var)
# dat2=data.frame(ratio="AlYb",
#                 n.fourier=length(res100$MTSE_full$freqs),
#                 tau=res100$avarOut$tau,
#                 avar=res100$avarOut$avar,
#                 var=res100$avarOut$var)
# 
# dat3=data.frame(ratio="AlYb",
#                 n.fourier=length(res500$MTSE_full$freqs),
#                 tau=res500$avarOut$tau,
#                 avar=res500$avarOut$avar,
#                 var=res500$avarOut$var)
#                 
# dat4=data.frame(ratio="AlYb",
#                 n.fourier=length(res1000$MTSE_full$freqs),
#                 tau=res1000$avarOut$tau,
#                 avar=res1000$avarOut$avar,
#                 var=res1000$avarOut$var)
# 
# allavar=bind_rows(dat1,dat2,dat3,dat4)