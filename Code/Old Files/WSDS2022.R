###### WSDS 2022 Talk Results ######

#### Multitaper Estimator Example ######

getwd()
clk <- read_csv(file = "Data/clockDataExample.csv", col_names = FALSE)
plot(clk$X1)


#introduce gaps
gappy_dat <- clk$X1
gappy_dat[c(100:500, 1300:2000, 3000:3400)] <- NA
plot(gappy_dat)

gappy_dat_MTSE <- multitaper_est(gappy_dat, NW = 3, K = 5)

dim(gappy_dat_MTSE$tapers)

taper1 <- c(gappy_dat_MTSE$tapers[1:99,1], rep(NA, times = 401),
            gappy_dat_MTSE$tapers[100:898,1], rep(NA, times = 2000-1300 + 1),
            gappy_dat_MTSE$tapers[899:1897,1], rep(NA, times = 3400-3000 + 1),
            gappy_dat_MTSE$tapers[1898:2497,1])
k <- 5
taper5 <- c(gappy_dat_MTSE$tapers[1:99,k], rep(NA, times = 401),
            gappy_dat_MTSE$tapers[100:898,k], rep(NA, times = 2000-1300 + 1),
            gappy_dat_MTSE$tapers[899:1897,k], rep(NA, times = 3400-3000 + 1),
            gappy_dat_MTSE$tapers[1898:2497,k])


par(mfrow = c(1,1))
plot(taper1, type = "l", lwd = 1, ylab = "k = 1", xlab = "")
plot(taper2, type = "l", lwd = 1, ylab = "k = 2", xlab = "")
plot(taper3, type = "l", lwd = 1, ylab = "k = 3", xlab = "")
plot(taper4, type = "l", lwd = 1, ylab = "k = 4", xlab = "")
plot(taper5, type = "l", lwd = 1, ylab = "k = 5", xlab = "")

tapers <- cbind(taper1,taper2,taper3,taper4,taper5)
plot(taper1, type = "l", lwd = 1, ylab = "kth taper", xlab = "", col = colors[1], ylim = c(-0.05,0.05))

for(i in 2:4){
  lines(tapers[,i], col = colors[i])
}


#clock noise

plot(clock_df[1:400], type = "l")

#clock noise with gaps

plot(clock_df[1:4000], type = "l", ylab = "Fractional Frequency Deviates", xlab = "Time [s]", cex.lab = 1.5)

#white noise
set.seed(124)
y.wn <- rnorm(200)
plot(y.wn, type = "l", ylab = "X(t)", xlab = "t", lwd = 2, cex.lab = 1.5)

plot(rep(1,times = 200), type = "l", ylab = "S(f)", xlab = "f", xaxt = "n", lwd = 2, cex.lab = 1.5)
axis(side = 1, labels = seq(0,0.5, length.out = 5), at = seq(0,200, length.out = 5))


#sines + cosines
t <- seq(0,10, length.out = 1000)
f1 <- 0.5
f2 <- 3
f3 <- 5
sin.wave.1 <- sin(2*pi*f1*t)
sin.wave.2 <- sin(2*pi*f3*t)
cos.wave.1 <- cos(2*pi*f2*t)
par(mfrow = c(1,1), lwd = 2, cex.lab = 1.8, mar = c(4,5,4,4))
plot(t, sin.wave.1, type = "l", col = "red", ylab = expression(paste("sin(2", pi, f[1],"t)")))
plot(t, sin.wave.2, type = "l", col = "red", ylab = expression(paste("sin(2", pi, f[3],"t)")))
plot(t, cos.wave.1, type = "l", col = "blue", ylab = expression(paste("cos(2", pi, f[2],"t)")))
plot(t, cos.wave.1 + sin.wave.1 + sin.wave.2, type = "l", col = "purple", ylab = expression(X[t]))



### periodogram/tapering plots

##see SAUTS_FigureRecreation.R File


###Multitaper (Chave)

t.n <- 1:14500
t.n <- t.n[-c(4745:5447, 8378:9545,12823:13051)]
W <- 0.00097
dist.mat <- rdist(t.n)

A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
A.prime[row(A.prime) == col(A.prime)] <- W*2
eigdec <- eigs_sym(A.prime, k = 15, which = "LM")

#eigdec <- eigen(A.prime, symmetric = TRUE)

##changing signs of the vectors
K = 15 #number of sequences we want
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




##from matlab:vk.mat, uk.mat
uk.mat <- readMat("C:/Users/cmb15/OneDrive - UCB-O365/NIST/ClockDataAnalysis/Data/uk.mat")
dim(uk.mat$u)


### Fig 1
par(mfrow = c(2,1), mar = c(4,5,2,4))
x.t_example <- rnorm(14500)
x.t_example[c(4745:5447, 8378:9545,12823:13051)] <- NA
plot(x.t_example, type = "l", ylab = expression(x[t]), xlab = "t")
colors = c("blue", "red", "green", "magenta", "cyan")
k = 0
s = 1
mdss1 <- eig_vecs[,s]#uk.mat$u[,s]#
mdss1_long <- c(mdss1[1:4744],rep(0,times = 5447-4745 + 1), mdss1[4745:(4745+8378-5447-1)],rep(0,times = 9545-8378 + 1), mdss1[7676:(7676 + 12823-9545-1)],rep(0,times = 13051-12823+1),mdss1[10954:length(mdss1)])
#p <- sqrt(sum(mdss1))
plot(mdss1_long,type = "l",ylim = c(-0.03,0.03),col = colors[1], xlab = "t", ylab = "Tapers") #sign(sum(mdss1))*

for(i in 2:5){
  mdss1 <- eig_vecs[,i] #uk.mat$u[,i]#
  mdss1_long <- c(mdss1[1:4744],rep(0,times = 5447-4745 + 1), mdss1[4745:(4745+8378-5447-1)],rep(0,times = 9545-8378 + 1), mdss1[7676:(7676 + 12823-9545-1)],rep(0,times = 13051-12823+1),mdss1[10954:length(mdss1)])
  lines(mdss1_long,type = "l",col = colors[i-5*k])
  
}





##### White/Flicker Noise data examples

set.seed(1)
white.noise <- rnorm(2048)
wn.gaps <- rnorm(2926)
wn.gaps[c(100:500,1300:1400, 1700:1874, 2000:2200)] <- NA

flicker.noise <- X.t_sims_flk2[1,]
FN.gaps <- X.t_sims_flk2_gps[1,]

par(mfrow = c(1,2))
plot(white.noise, ylab = "", main  = "White Noise", xlab = "t", type = "l")
plot(wn.gaps, ylab = "", main = "White Noise with Gaps", xlab = "t", type = "l")

plot(flicker.noise, ylab = "", main  = "Flicker Noise", xlab = "t", type = "l")
plot(FN.gaps, ylab = "", main = "Flicker Noise with Gaps", xlab = "t", type = "l")






