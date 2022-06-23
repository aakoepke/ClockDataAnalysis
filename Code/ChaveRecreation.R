##### Chave 2019 ####


#time sequence with gaps
t.n <- 1:14500
t.n <- t.n[-c(4745:5447, 8378:9545,12823:13051)]
W <- 0.00097
dist.mat <- rdist(t.n)

A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
A.prime[row(A.prime) == col(A.prime)] <- W*2
eigdec <- eigs_sym(A.prime, k = 15, which = "LM")

eigdec <- eigen(A.prime, symmetric = TRUE)

##from matlab:vk.mat
uk.mat <- readMat("C:/Users/cmb15/OneDrive - UCB-O365/NIST/ClockDataAnalysis/Data/uk.mat")
dim(uk.mat$u)

colors = c("blue", "red", "green", "magenta", "cyan")
k = 1
s = 6
mdss1 <- uk.mat$u[,s]#eigdec$vectors[,s]
mdss1_long <- c(mdss1[1:4744],rep(0,times = 5447-4745 + 1), mdss1[4745:(4745+8378-5447-1)],rep(0,times = 9545-8378 + 1), mdss1[7676:(7676 + 12823-9545-1)],rep(0,times = 13051-12823+1),mdss1[10954:length(mdss1)])
#p <- sqrt(sum(mdss1))
plot(mdss1_long,type = "l",ylim = c(-0.03,0.03),col = colors[1]) #sign(sum(mdss1))*

for(i in 7:10){
  mdss1 <- uk.mat$u[,i]#eigdec$vectors[,i]
  mdss1_long <- c(mdss1[1:4744],rep(0,times = 5447-4745 + 1), mdss1[4745:(4745+8378-5447-1)],rep(0,times = 9545-8378 + 1), mdss1[7676:(7676 + 12823-9545-1)],rep(0,times = 13051-12823+1),mdss1[10954:length(mdss1)])
  lines(mdss1_long,type = "l",col = colors[i-5*k])
  
}


mdss1 <- eigdec$vectors[,6]
mdss1_long <- c(mdss1[1:4744],rep(0,times = 5447-4745 + 1), mdss1[4745:(4745+8378-5447-1)],rep(0,times = 9545-8378 + 1), mdss1[7676:(7676 + 12823-9545-1)],rep(0,times = 13051-12823+1),mdss1[10954:length(mdss1)])
plot(mdss1_long,type = "l", ylim = c(-0.2,0.2))
abline(v = 4745)
abline(v = 5447)