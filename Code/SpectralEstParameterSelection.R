########## Chave Parameter Selection ###############
### This is a script to look at the bias properties of the
### MTSE for missing data



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

try2$e.values
try7$e.values
e.values.mat <- matrix(NA, ncol = 3, nrow = 12)

plot(try7$tapers[,3])

for(i in 1:8){
  tmp <- get_tapers(t.n, W = i/N, K = 3)
  e.values.mat[i,] <- tmp$e.values
}
e.values.mat <- e.values.mat[-c(9:12),]

bandwidths <- 1:8
numTapers <- 1:10
BV.mat <- matrix(NA, ncol = length(numTapers), nrow = length(bandwidths))

for(i in bandwidths){
  for(j in numTapers){
    print(i,j)
    tmp <- get_tapers(t.n, W = i/N, K = j)
    BV.mat[i,j] <- (1/j)*sum(1-tmp$e.values) 
  }
}


A <- matrix(BV.mat,ncol = 1)
A <- matrix(1:9, ncol = 3)
A



x <- rep(1:10, each = 8)
y <- rep(1:8, times = 10)
z <- matrix(BV.mat, ncol = 1)

quilt.plot(x,y,z)


#leakage
leakage_fun <- function(x){
  sum(1-x)
}

(1/3)*sum(1-try3$e.values)

#variance
1/3
