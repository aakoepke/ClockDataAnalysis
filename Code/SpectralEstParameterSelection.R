########## Chave Parameter Selection ###############
### This is a script to look at the bias properties of the
### MTSE for missing data

source("Code/SA_ImportantFunctions.R")

##### length: 2048 without missing pieces
2048*1.2

t.n <- 1:2500
t.n[c(100:200, 1050:1350,1800:1849)] <- NA
N <- length(na.omit(t.n))

try1 <- get_tapers(t.n, W = 1/N, K = 3)
try2 <- get_tapers(t.n, W = 2/N, K = 3)
try3 <- get_tapers(t.n, W = 3/N, K = 3)
try4 <- get_tapers(t.n, W = 4/N, K = 3)
try5 <- get_tapers(t.n, W = 5/N, K = 3)
try6 <- get_tapers(t.n, W = 6/N, K = 3)
try7 <- get_tapers(t.n, W = 7/N, K = 3)

try1$e.values
try2$e.values
e.values.mat <- matrix(NA, ncol = 3, nrow = 12)

plot(try7$tapers[,3])

for(i in 1:8){
  tmp <- get_tapers(t.n, W = i/N, K = 3)
  e.values.mat[i,] <- tmp$e.values
}
e.values.mat <- e.values.mat[-c(9:12),]

bandwidths <- 1:12
numTapers <- 1:23
e.values.mat <- matrix(NA, ncol = length(numTapers), nrow = length(bandwidths))
e.vectors.matlist <- list()

#calculate eigenvalues/vectors up to 15
for(i in 1:length(bandwidths)){
  print(i)
  tmp <- get_tapers(t.n, W = i/N, K = length(numTapers))
  e.values.mat[i,] <- tmp$e.values
  e.vectors.matlist[[i]] <- tmp$tapers
}

#go along and calculate spectral leakage
spec.leakage.mat <- matrix(NA, nrow = length(bandwidths), ncol = length(numTapers))
for(i in 1:length(bandwidths)){
  for(j in 1:length(numTapers)){
    spec.leakage.mat[i,j] <- leakage_fun(e.values.mat[i,1:j])
  }
}


x <- rep(1:length(numTapers), each = length(bandwidths))
y <- rep(1:length(bandwidths), times = length(numTapers))
z <- matrix(spec.leakage.mat, ncol = 1)

library(fields)
bubblePlot(x,y,z, size = 1.5, xlab = "# of Tapers", ylab = "NW")
points(y = bandwidths, x = 2*bandwidths-1, pch = 1, cex = 3)
apply(spec.leakage.mat, FUN = which.min, MARGIN = 2)






#leakage
leakage_fun <- function(x){
  (1/length(x))*sum(1-x)
}

(1/3)*sum(1-try3$e.values)

#variance
1/3



########### Generate Spectral Estimates with Different Parameter Values #############

#generate series with known spectrum
# p <- noise(kind = "pink", samp.rate = 2500)
# p.vec <- p@left
# length(p.vec)
# plot(p.vec)
# p.vec_missing <- p.vec
# p.vec_missing[c(100:200, 1050:1350,1800:1849)] <- NA

library(arfima)
source("Code/SA_ImportantFunctions.R")

t.n <- t.n_missing <- 1:2500
t.n_missing[c(100:200, 1050:1350,1800:1849)] <- NA

N <- length(t.n)
x.t <- arfima.sim(n = N, model = list(dfrac = 0.25))
x.t[c(100:200, 1050:1350,1800:1849)] <- NA
plot(x.t)

N <- length(na.omit(t.n_missing))
par(mfrow = c(3,3))
for(i in 1:15){
tapers_arfima <- get_tapers(t.n = t.n_missing, W = 2/N, K = i)
mtse_arfima <- MT_spectralEstimate(X.t = x.t, V.mat = tapers_arfima$tapers)
plot(log(mtse_arfima$freqs[-1]), log(mtse_arfima$spectrum[-1]))
lines(log(freqs[-1]), log(arfima.spec(freqs[-1], 1, d = 0.25)), col = "blue")
}


#generate 300 simulated arfima processes
X.t_sims_ARFIMA <- matrix(NA, nrow = 300, ncol = 2500)

for(i in 1:300){
set.seed(i)
  X.t_sims_ARFIMA[i,] <- arfima.sim(n = 2500, model = list(dfrac = 0.25))
}

X.t_sims_ARFIMA[,c(100:200, 1050:1350,1800:1849)] <- NA

#calculate their spectrum with various choices of parameters
ARFIMA_specs4 <- matrix(NA, nrow = 300, ncol = 1025)
tapers_arfima4 <- get_tapers(t.n = t.n_missing, W = 4/N, K = 15)
spectrum_saved4 <- list()

for(j in 2:15){
  print(j)
for(i in 1:300){
  x.t <- X.t_sims_ARFIMA[i,]
  tmp <- MT_spectralEstimate(X.t = x.t, V.mat = tapers_arfima4$tapers[,1:j])
  ARFIMA_specs4[i,] <- tmp$spectrum
  }
  spectrum_saved4[[j]] <- ARFIMA_specs4
}


#plot their distributions

#tidy the data
A <- matrix(1:80, nrow = 8, ncol = 10)
A.vec <- matrix(spectrum_saved[[2]], ncol = 1)
freqs <- seq(0,0.5, length.out = 1025)
freqs.rep <- rep(freqs, each = 300)
spec.df <- cbind(A.vec,freqs.rep)
spec.df %<>% as.data.frame()
spec.df$freqs.rep <- as.factor(spec.df$freqs.rep)
##plot as violins

plot.violins <- ggplot(data = spec.df, aes(x = freqs.rep, y = spectrum)) + 
  geom_violin()
plot.violins
#not super useful

##calculate AVAR using these
numberOfSimulations = 300
trfunc.vec <- bpvar.vec <- rep(NA, times = numberOfSimulations)
taus <- c(2^(0:9), floor(N/3))
tmat <- bmat <- matrix(NA, ncol = numberOfSimulations, nrow = length(taus))

f <- seq(0,0.5,length.out = N/2 + 1)
avar_calcs_tmat <- avar_calcs_bmat <- list()

for(j in 2:15){
  r = 0
for(k in taus){
  r = r + 1
  tau = k
  print(paste("r = ", r))
  
  for(i in 1:300){
    print(i)
    set.seed(i)
    
    #calculate S.hat
    MTSE_full <- spectrum_saved4[[j]][i,]
    
    #calculate bandpass variance
    temp_bp <- integrate(approxfun(f, MTSE_full), lower = 1/(4*tau), upper = 1/(2*tau), subdivisions = 1000)
    bpvar.vec[i] <- 4*temp_bp$value
    
    #calculate transfer function AVAR
    G.vec <- transfer.func(f, tau)
    G.vec[1] <- 0
    trfunc.vec[i] <- f[2]*sum(G.vec*MTSE_full)
    
  }
  tmat[r,] <- trfunc.vec
  bmat[r,] <- bpvar.vec
}

avar_calcs_tmat[[j]] <- tmat
avar_calcs_bmat[[j]] <- bmat

}

####plots

##tidy the data
tmat2 <- avar_calcs_tmat[[2]]
tmat10 <- avar_calcs_tmat[[10]]
tmat15 <- avar_calcs_tmat[[15]]
tmat2 %<>% t()
tmat10 %<>% t()
tmat15 %<>% t()
dim(bmat)
three.mat <- rbind(tmat2, tmat10,tmat15)
method_labels <- rep(c("a", "b", "c"), each = numberOfSimulations)
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

# breaks <- 10^(-10:10)
# minor_breaks <- rep(1:9,21 )*rep(10^(-10:10), each = 9)
# 
# linedat = data.frame(tau=2:floor(N/3))
# linedat$truth = tavar_ARFIMA(floor(N/3), d = 0.25, sig.2.a = 1)
# linedat$method = NA

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat,aes(tau,truth), linewidth = 1.2)+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("2 tapers","10 tapers", "15 tapers"),name="Method")+
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
  xlab(expression(tau)) #+
#ggtitle("White Noise Comparison, No Gaps")





################# With Data Spacing ####################
## Repeat analysis above with spacing similar to that ##
## seen in real clock data                            ##
########################################################
t.n <- 1:9000
t.n[c(100:200, 1050:1350,1800:1849, 3000:7000, 6010:6500)] <- NA
N <- length(na.omit(t.n))
N

bandwidths <- 1:12
numTapers <- 1:23
e.values.mat <- matrix(NA, ncol = length(numTapers), nrow = length(bandwidths))
e.vectors.matlist <- list()

#calculate eigenvalues/vectors up to 15
for(i in 1:length(bandwidths)){
  print(i)
  tmp <- get_tapers(t.n, W = i/N, K = length(numTapers))
  e.values.mat[i,] <- tmp$e.values
  e.vectors.matlist[[i]] <- tmp$tapers
}

#go along and calculate spectral leakage
spec.leakage.mat <- matrix(NA, nrow = length(bandwidths), ncol = length(numTapers))
for(i in 1:length(bandwidths)){
  for(j in 1:length(numTapers)){
    spec.leakage.mat[i,j] <- leakage_fun(e.values.mat[i,1:j])
  }
}


x <- rep(1:length(numTapers), each = length(bandwidths))
y <- rep(1:length(bandwidths), times = length(numTapers))
z <- matrix(spec.leakage.mat, ncol = 1)

library(fields)
bubblePlot(x,y,z, size = 1.5, xlab = "# of Tapers", ylab = "NW")
points(y = bandwidths, x = 2*bandwidths-1, pch = 1, cex = 3)
points(y = bandwidths, x = 2*bandwidths-1, pch = 2, cex = 3)
apply(spec.leakage.mat, FUN = which.min, MARGIN = 1)

which(e.values.mat > 0.9)


library(arfima)
N.full <- 9000
t.n <- t.n_missing <- 1:N.full
t.n_missing[c(100:200, 1050:1350,1800:1849, 3000:7000, 6010:6500)] <- NA

N <- length(t.n)
x.t <- arfima.sim(n = N, model = list(dfrac = 0.25))
x.t[c(100:200, 1050:1350,1800:1849, 3000:7000, 6010:6500)] <- NA
plot(x.t)
tapers_arfima <- get_tapers(t.n = t.n_missing, W = 2/N, K = 3)
plot(tapers_arfima$tapers[,1])
tapers_arfima$e.values
N <- length(na.omit(t.n_missing))
par(mfrow = c(3,3))
for(i in 1:15){
  tapers_arfima <- get_tapers(t.n = t.n_missing, W = 2/N, K = i)
  mtse_arfima <- MT_spectralEstimate(X.t = x.t, V.mat = tapers_arfima$tapers)
  plot(log(mtse_arfima$freqs[-1]), log(mtse_arfima$spectrum[-1]))
  lines(log(freqs[-1]), log(arfima.spec(freqs[-1], 1, d = 0.25)), col = "blue")
}


#generate 300 simulated arfima processes
X.t_sims_ARFIMA2 <- matrix(NA, nrow = 300, ncol = N.full)

for(i in 1:300){
  set.seed(i)
  X.t_sims_ARFIMA2[i,] <- arfima.sim(n = N.full, model = list(dfrac = 0.25))
}

X.t_sims_ARFIMA2[,c(100:200, 1050:1350,1800:1849, 3000:5000, 6010:6500)] <- NA

#calculate their spectrum with various choices of parameters
ARFIMA_specs2 <- matrix(NA, nrow = 300, ncol = floor(N/2) + 1)
tapers_arfima2 <- get_tapers(t.n = t.n_missing, W = 2/N, K = 15)
spectrum_saved2 <- list()

for(j in 2:15){
  print(j)
  for(i in 1:300){
    print(i)
    x.t <- X.t_sims_ARFIMA2[i,]
    tmp <- MT_spectralEstimate(X.t = x.t, V.mat = tapers_arfima2$tapers[,1:j])
    ARFIMA_specs2[i,] <- tmp$spectrum
  }
  spectrum_saved2[[j]] <- ARFIMA_specs2
}



