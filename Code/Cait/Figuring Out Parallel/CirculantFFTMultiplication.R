######## testing out fft mutliplication #########


##### Toy Example Matrices ######
circulant <- function(vec) {
  n <- length(vec)
  mat <- matrix(0, n, n)
  for (i in 1:n) {
    mat[i, ] <- vec
    vec <- c(tail(vec, 1), head(vec, -1))
  }
  return(mat)
}

C <- circulant(c(1:10000))
dim(C)

C.mult <- C[1:(dim(C)[1]/2),1:(dim(C)[1]/2)]
dim(C.mult)
B <- matrix(rnorm(dim(C)[1]*5), nrow = dim(C)[1]/2, ncol = 5)
dim(B)

C%*%B

####### Now we do with FFT of C #####

# C = F^{-1}\Lambda F
#where \Lambda = diag(Fc) where c = C[,1]

c_vec <- C[,1]
#for a vector c_vec
Fc <- (1/length(c_vec))*fft(c_vec)
#
B.zero <- rbind(B, matrix(0, nrow = length(c_vec) - dim(B)[1], ncol = dim(B)[2]))
FB <- apply(B.zero, MARGIN = 2, FUN = fft)

CB <- apply(FB, MARGIN = 2, FUN = function(x) {Re(fft(Fc*x, inverse = TRUE))})
CB

### time test ###
### regular
startTime = Sys.time()
for(i in 1:10){
  C.mult%*%B
}
print(Sys.time() - startTime)


###FFT
mult_fun <- function(x){
  Re(fft(Fc*x, inverse = TRUE))
}

startTime = Sys.time()
for(i in 1:10){
  apply(FB, MARGIN = 2, FUN = mult_fun)
}
print(Sys.time() - startTime)


### try this on a real example to make sure they match up

N=5000
t.vec <- 1:N  # Time vector

## tapers
setW = 8/N
setK = 5
taperObject <- get_tapers(t.vec, W = setW, K = setK)
taperMat <- taperObject$tapers


tN = 10
R_mat <- toeplitz(c(seq(1,0.1, length.out = tN), rep(0, times = N-tN))) #to start
c_vec <- c(R_mat[,1], 0, rev(R_mat[1,2:N]))
length(c_vec)

freq <- seq(0, .5, length.out = N)  # Frequency vector

V_star_mat <- t(taperMat*exp(1i*2*pi*freq[3]*t.vec))
V_mat <- taperMat*exp(-1i*2*pi*freq[505]*t.vec)

B = V_mat
##using toeplitz
Fc <- (1/length(c_vec))*fft(c_vec)
B.zero <- rbind(B, matrix(0, nrow = length(c_vec) - dim(B)[1], ncol = dim(B)[2]))
FB <- apply(B.zero, MARGIN = 2, FUN = fft)

CB <- apply(FB, MARGIN = 2, FUN = function(x) {fft(Fc*x, inverse = TRUE)})
dim(CB)
CB <- CB[1:N,]
dim(CB)

RV <- R_mat%*%V_mat
dim(RV)
RV[1:10,1]
CB[1:10,1]

sum(abs(V_star_mat%*%RV-V_star_mat%*%CB))

norm(V_star_mat%*%RV/K, type = "2")
norm(V_star_mat%*%CB/K, type = "2")

### time test ###
### regular
startTime = Sys.time()
for(i in 1:10){
  R_mat%*%V_mat
}
print(Sys.time() - startTime)


###FFT
mult_fun <- function(x){
  fft(Fc*x, inverse = TRUE)
}

startTime = Sys.time()
for(i in 1:10){
  apply(FB, MARGIN = 2, FUN = mult_fun)
}
print(Sys.time() - startTime)

