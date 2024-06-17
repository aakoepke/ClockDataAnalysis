library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 6+1
mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)


##################################################10/13 res


##########srYb res

resName="Data/resultsForSrYb1shortday_3391_Oct13.Rds"
res1=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res1$V.mat$e.values

resName="Data/resultsForSrYb1longday_22304_Oct13.Rds"
res2=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res2$V.mat$e.values

resName="Data/resultsForSrYb2days_30309_Oct13.Rds"
res3=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res3$V.mat$e.values

resName="Data/resultsForSrYb3days_33699_Oct13.Rds"
res4=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res4$V.mat$e.values



dat1=data.frame(ratio="SrYb_shortDay",
              n.fourier=length(res1$MTSE_full$freqs),
              W=res1$W,
              freq=res1$MTSE_full$freqs,
              spectrum=res1$MTSE_full$spectrum)
dat2=data.frame(ratio="SrYb_longDay",
              n.fourier=length(res2$MTSE_full$freqs),
              W=res2$W,
              freq=res2$MTSE_full$freqs,
              spectrum=res2$MTSE_full$spectrum)
dat3=data.frame(ratio="SrYb_2days",
              n.fourier=length(res3$MTSE_full$freqs),
              W=res3$W,
              freq=res3$MTSE_full$freqs,
              spectrum=res3$MTSE_full$spectrum)
dat4=data.frame(ratio="SrYb_3days",
              n.fourier=length(res4$MTSE_full$freqs),
              W=res4$W,
              freq=res4$MTSE_full$freqs,
              spectrum=res4$MTSE_full$spectrum)

allSpec=bind_rows(dat1,dat2,dat3,dat4)

ggplot(allSpec,aes(freq,spectrum,col=factor(interaction(n.fourier,W))))+
  geom_line()+
  scale_y_log10()+
  scale_x_log10()
  # geom_smooth()
ggplot(allSpec,aes(freq,spectrum,col=factor(interaction(n.fourier,W))))+
  # geom_line()+
  scale_y_log10()+
  scale_x_log10()+
  geom_smooth()
ggplot(allSpec,aes(freq,spectrum))+
  geom_line()+
  facet_wrap(~n.fourier+W)+
  geom_smooth()+
  scale_y_log10()+
  scale_x_log10()


###################look at tapers

# taperDF=data.frame(res4$V.mat$tapers)
# colnames(taperDF)=c("Taper 1",
#                     "Taper 2",
#                     "Taper 3",
#                     "Taper 4",
#                     "Taper 5",
#                     "Taper 6")
# taperDF$t=res4$t.vec
# 
# taperDFlong=melt(taperDF,id.vars = "t",variable.name = "tapers")
# taperDFlong$x=res4$x.t
# 
# ggplot(taperDFlong,aes(t,value,col=tapers))+
#   geom_line()
# 
# ggplot(filter(taperDFlong,tapers=="Taper 1"),aes(t,x,col=tapers))+
#   geom_line()
# 
# 
# ### 2 day tapers
# 
# taperDF=data.frame(res3$V.mat$tapers)
# colnames(taperDF)=c("Taper 1",
#                     "Taper 2",
#                     "Taper 3",
#                     "Taper 4",
#                     "Taper 5",
#                     "Taper 6")
# taperDF$t=res3$t.vec
# 
# taperDFlong=melt(taperDF,id.vars = "t",variable.name = "tapers")
# 
# ggplot(taperDFlong,aes(t,value,col=tapers))+
#   geom_line()


### 1 day tapers

taperDF=data.frame(res1$V.mat$tapers)
colnames(taperDF)=c("1","2","3","4","5","6")
taperDF$t=res1$t.vec

dataDF=data.frame(t=res1$t.vec,
                  value=res1$x.t)
###
alltData=data.frame(t=min(taperDF$t):max(taperDF$t))

taperDF=merge(alltData,taperDF,by="t",all.x = T)

dataDF=merge(alltData,dataDF,by="t",all.x = T)
dataDF$type="Data"
dataDF$Taper="Data"

taperDFlong=melt(taperDF,id.vars = "t",variable.name = "Taper")
taperDFlong$type="Tapers"

allTaperDat=bind_rows(taperDFlong,dataDF)

pdf("Plots/SrYbShortDayTapers.pdf",width = 6.5,height = 4)
ggplot(allTaperDat,aes(t,value,col=Taper))+
  geom_line()+
  facet_wrap(~type,nrow = 2,scales = "free_y",
             strip.position = "left", 
             labeller = as_labeller(c(Data = "Clock Ratio Data", Tapers = "Tapers") ) )+
  # scale_color_brewer(palette = "Accent")+
  scale_color_manual(values=mycolors)+
  # theme(axis.title.y=element_blank())+
  theme(legend.position = "none",legend.title=element_blank(),strip.background = element_blank(),
        strip.placement = "outside")+
  ylab(NULL)
dev.off()

#### 1 long day

taperDF=data.frame(res2$V.mat$tapers)
colnames(taperDF)=c("1","2","3","4","5","6")
taperDF$t=res2$t.vec

dataDF=data.frame(t=res2$t.vec,
                  value=res2$x.t)
###
alltData=data.frame(t=min(taperDF$t):max(taperDF$t))

taperDF=merge(alltData,taperDF,by="t",all.x = T)

dataDF=merge(alltData,dataDF,by="t",all.x = T)
dataDF$type="Data"
dataDF$Taper="Data"

taperDFlong=melt(taperDF,id.vars = "t",variable.name = "Taper")
taperDFlong$type="Tapers"

allTaperDat=bind_rows(taperDFlong,dataDF)

pdf("Plots/SrYbLongDayTapers.pdf",width = 6.5,height = 4)
ggplot(allTaperDat,aes(t,value,col=Taper))+
  geom_line()+
  facet_wrap(~type,nrow = 2,scales = "free_y",
             strip.position = "left", 
             labeller = as_labeller(c(Data = "Clock Ratio Data", Tapers = "Tapers") ) )+
  # scale_color_brewer(palette = "Accent")+
  scale_color_manual(values=mycolors)+
  # theme(axis.title.y=element_blank())+
  theme(legend.position = "none",legend.title=element_blank(),strip.background = element_blank(),
        strip.placement = "outside")+
  ylab(NULL)

dev.off()


### 3 day data

taperDF=data.frame(res4$V.mat$tapers)
colnames(taperDF)=c("1","2","3","4","5","6")
taperDF$t=res4$t.vec

dataDF=data.frame(t=res4$t.vec,
                  value=res4$x.t)
###
alltData=data.frame(t=min(taperDF$t):max(taperDF$t))

taperDF=merge(alltData,taperDF,by="t",all.x = T)

dataDF=merge(alltData,dataDF,by="t",all.x = T)
dataDF$type="Data"
dataDF$Taper="Data"

taperDFlong=melt(taperDF,id.vars = "t",variable.name = "Taper")
taperDFlong$type="Tapers"

allTaperDat=bind_rows(taperDFlong,dataDF)

pdf("Plots/SrYb3DaysTapers.pdf",width = 6.5,height = 4)
ggplot(allTaperDat,aes(t,value,col=Taper))+
  geom_line()+
  facet_wrap(~type,nrow = 2,scales = "free_y",
             strip.position = "left", 
             labeller = as_labeller(c(Data = "Clock Ratio Data", Tapers = "Tapers") ) )+
  # scale_color_brewer(palette = "Accent")+
  scale_color_manual(values=mycolors)+
  # theme(axis.title.y=element_blank())+
  theme(legend.position = "none",legend.title=element_blank(),strip.background = element_blank(),
        strip.placement = "outside")+
  ylab(NULL)

dev.off()

###################


dat1=res1$allAvarRes
dat1$Ratio="Short Day"

dat2=res2$allAvarRes
dat2$Ratio="Long Day"

dat3=filter(res3$allAvarRes,calculation=="spectral")
dat3$Ratio="2 days"

dat4=filter(res4$allAvarRes,calculation=="spectral")
dat4$Ratio="3 days"

allavar=bind_rows(dat1,dat2)#,dat3,dat4)
allavar$Ratio=factor(allavar$Ratio,levels = c("Short Day","Long Day", "2 days", "3 days"))
order(allavar$Ratio)

allavar$Method=NA
allavar$Method[allavar$calculation=="overlapping"]="Current"
allavar$Method[allavar$calculation=="spectral"]="Spectral"

allavar$Length=NA
allavar$Length[allavar$Ratio=="Short Day"]="Short Day"
allavar$Length[allavar$Ratio=="Long Day"]="Long Day"

pdf("Plots/SrYb1dayAVARcomparison.pdf",width = 6,height = 3)
ggplot(filter(allavar,calculation!="old"),aes(tau,avar,col=Method,shape=Length))+
  geom_point()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  # theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))
  # facet_wrap(~Ratio)
dev.off()

allavar=bind_rows(dat1,dat2,dat3,dat4)
allavar$Ratio=factor(allavar$Ratio,levels = c("short day","long day", "2 days", "3 days"),ordered = T)

pdf("Plots/SrYb3dayAVARcomparison_split.pdf",width = 6.5,height = 2.9)
ggplot(filter(allavar,calculation!="old"),aes(Ratio,avar,col=calculation))+
  geom_point()+
  # geom_jitter(width = .0001,height = 0)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  # theme(legend.position = c(.15, .2))+
  # scale_y_log10()+
  # scale_x_log10()+
  # annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  # xlab("")+
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1),legend.position="bottom",axis.title.x=element_blank())+
  # theme(
  #   axis.text.x = element_blank(),
  #   axis.text.y = element_blank(),
  #   axis.ticks = element_blank())+
  facet_wrap(~tau,scales = "free_y",nrow = 2)
  # theme(strip.text.x = element_text(angle = 45)) 
dev.off()

pdf("Plots/SrYb3dayAVARcomparison.pdf",width = 6,height = 2.3)
ggplot(filter(allavar,calculation!="old"),aes(tau,avar,col=calculation,shape=Ratio))+
  geom_point()+
  # geom_jitter(width = .0001,height = 0)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  # theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  theme(legend.title=element_blank())
# facet_wrap(~tau,scales = "free")
dev.off()



##### look at spectrum for 3 ratios, longest series

resName="Data/resultsForSrYb3days_33699_Oct13.Rds"
res1=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res1$V.mat$e.values
resName="Data/resultsForAlSr3days_16053_Oct13.Rds"
res2=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res2$V.mat$e.values
resName="Data/resultsForAlYb3days_17982_Oct13.Rds"
res3=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res3$V.mat$e.values


dat1=data.frame(ratio="SrYb",
                n.fourier=length(res1$MTSE_full$freqs),
                W=res1$W,
                freq=res1$MTSE_full$freqs,
                spectrum=res1$MTSE_full$spectrum)

dat2=data.frame(ratio="AlSr",
                n.fourier=length(res2$MTSE_full$freqs),
                W=res2$W,
                freq=res2$MTSE_full$freqs,
                spectrum=res2$MTSE_full$spectrum)
dat3=data.frame(ratio="AlYb",
                n.fourier=length(res3$MTSE_full$freqs),
                W=res3$W,
                freq=res3$MTSE_full$freqs,
                spectrum=res3$MTSE_full$spectrum)



allSpec=bind_rows(dat1,dat2,dat3)


##work on this plot!!!!!!!!!!!!!!
pdf("Plots/threeRatiosSpectrumComparison.pdf",width = 6,height = 4)
ggplot(allSpec,aes(freq,spectrum,col=ratio))+
  geom_line(alpha=.5)+
  scale_y_log10()+
  xlab("Frequency")+
  ylab("S")+
  # scale_x_log10()
  geom_smooth()
  # geom_hline(yintercept = .05)
dev.off()

ggplot(allSpec,aes(freq,spectrum,col=ratio))+
  geom_line(alpha=.4)+
  scale_y_log10()+
  scale_x_log10()+
  geom_smooth()+
  xlim(c(1e-02,.5))

ggplot(allSpec,aes(freq,spectrum))+
  geom_line()+
  facet_wrap(~ratio)+
  geom_smooth()+
  scale_y_log10()+
  scale_x_log10()


dat1=res1$allAvarRes
dat1$ratio="SrYb"

dat2=res2$allAvarRes
dat2$ratio="AlSr"

dat3=res3$allAvarRes#filter(res3$allAvarRes,calculation=="spectral")
dat3$ratio="AlYb"

allavar=bind_rows(dat1,dat2,dat3)


pdf("Plots/threeRatios3daysAVARcomparison.pdf",width = 6,height = 4)
ggplot(filter(allavar,calculation!="old"),aes(tau,avar,shape=calculation,col=ratio))+
  geom_point()+
  ### add true straight line below
  geom_abline(slope = -1,intercept = 0,size=1)+
  # theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))
# facet_wrap(~ratio)
dev.off()


###################################################################################
########################### looking at part of long day

resName="Data/resultsForSrYb1shortday_3391_Oct13.Rds"
res1=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res1$V.mat$e.values

resName="Data/resultsForSrYb1longday_22304_Oct13.Rds"
res2=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res2$V.mat$e.values

resName="Data/resultsForSrYb1longday_shortDayLength_3391_Oct15.Rds"
res3=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res3$V.mat$e.values




dat1=data.frame(ratio="SrYb_shortDay",
                n.fourier=length(res1$MTSE_full$freqs),
                W=res1$W,
                freq=res1$MTSE_full$freqs,
                spectrum=res1$MTSE_full$spectrum)
dat2=data.frame(ratio="SrYb_longDay",
                n.fourier=length(res2$MTSE_full$freqs),
                W=res2$W,
                freq=res2$MTSE_full$freqs,
                spectrum=res2$MTSE_full$spectrum)
dat3=data.frame(ratio="SrYb_longdayPart",
                n.fourier=length(res3$MTSE_full$freqs),
                W=res3$W,
                freq=res3$MTSE_full$freqs,
                spectrum=res3$MTSE_full$spectrum)
allSpec=bind_rows(dat1,dat2,dat3)

ggplot(allSpec,aes(freq,spectrum,col=ratio))+
  geom_line()+
  scale_y_log10()+
  scale_x_log10()
# geom_smooth()
ggplot(allSpec,aes(freq,spectrum))+
  geom_line()+
  facet_wrap(~ratio)+
  geom_smooth()+
  scale_y_log10()+
  scale_x_log10()


dat1=res1$allAvarRes
dat1$ratio="SrYb_shortDay"

dat2=res2$allAvarRes
dat2$ratio="SrYb_longDay"

dat3=res3$allAvarRes
dat3$ratio="SrYb_longDayPart"

allavar=bind_rows(dat1,dat2,dat3)


# pdf("Plots/SrYb1dayAVARcomparison.pdf",width = 6,height = 4)
ggplot(filter(allavar,calculation!="old"),aes(tau,avar,col=calculation,shape=ratio))+
  geom_point()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  # theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))
# facet_wrap(~ratio)
# dev.off()

allavar=bind_rows(dat1,dat2,dat3,dat4)

# pdf("Plots/SrYb3dayAVARcomparison.pdf",width = 6,height = 4)
ggplot(filter(allavar,calculation!="old"),aes(tau,avar,col=calculation,shape=ratio))+
  # geom_point()+
  geom_jitter(width = .0001,height = 0)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  # theme(legend.position = c(.15, .2))+
  # scale_y_log10()+
  # scale_x_log10()+
  # annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  facet_wrap(~tau,scales = "free")
# dev.off()











############################################################################################
############################################################################################
############################################################################################

resName="Data/resultsForSrYb040318_16839_Oct10.Rds"
res1=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res1$V.mat$e.values

resName="Data/resultsForAlYb040318_11443_Oct10.Rds"
res2=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res2$V.mat$e.values

resName="Data/resultsForAlSr040318_10384_Oct10.Rds"
res3=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res3$V.mat$e.values


res=res3
plot(log10(res$MTSE_full$freqs),log10(res$MTSE_full$spectrum),type="l")

ggplot(res$allAvarRes,aes(tau,avar,col=calculation))+
  geom_point()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))

plot(res$V.mat$tapers[,1])



##############2 days

resName="Data/resultsForSrYb2days_25694_Oct10.Rds"
res1=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res1$V.mat$e.values

resName="Data/resultsForAlYb2days_17982_Oct10.Rds"
res2=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res2$V.mat$e.values

resName="Data/resultsForAlSr2days_16053_Oct10.Rds"
res3=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res3$V.mat$e.values


res=res3
plot(log10(res$MTSE_full$freqs),log10(res$MTSE_full$spectrum),type="l")

ggplot(res$allAvarRes,aes(tau,avar,col=calculation))+
  geom_point()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))





##############3 days

resName="Data/resultsForSrYb3days_33699_Oct12.Rds"
res1=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res1$V.mat$e.values

resName="Data/resultsForAlYb3days_17982_Oct12.Rds"
res2=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res2$V.mat$e.values

resName="Data/resultsForAlSr3days_16053_Oct12.Rds"
res3=readRDS(file = paste("/home/aak3/NIST/ClockDataAnalysis/",resName,sep=""))
res3$V.mat$e.values


res=res3
plot(log10(res$MTSE_full$freqs),log10(res$MTSE_full$spectrum),type="l")

ggplot(res$allAvarRes,aes(tau,avar,col=calculation))+
  geom_point()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))
