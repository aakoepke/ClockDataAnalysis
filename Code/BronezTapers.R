################ Bronez tapers ################

get_k_eigs=function(Ra,thek,f.w){
  evs=eigs(Ra,k = thek)
  
  return(list("weights" = evs$vectors*sqrt(2*f.w), "eigenvalues" = evs$values))
}

get_k_geigen=function(Ra,Rb,thek,f.w){
  evs=geigen(Ra,Rb, symmetric = TRUE)
  
  N=dim(evs$vectors)[1]
  return(list("weights" = evs$vectors[,(N-thek + 1):N]*sqrt(2*f.w), "eigenvalues" = sort(evs$values, decreasing = TRUE)[1:thek]))
}

get.weights_bronez <- function(t.n = 1:50, K = 1, f.c = 0, f.w){
  
  dist.mat <- #rdist(t.n)
    outer(t.n,t.n,"-")
  
  # B = [-pi,pi] for omega or [-1/2,1/2] for f
  R.b <- 1/(pi*(dist.mat))*(sin(dist.mat*pi))
  R.b[row(R.b) == col(R.b)] <- 1
  # B = identity
  
  j = complex(real = 0, imaginary = 1)
  R.a <- exp(j*2*pi*f.c*dist.mat)*(sin(2*pi*f.w*(dist.mat))/(pi*dist.mat))
  R.a[row(R.a) == col(R.a)] <- f.w*2
  
  #Solve the eigenvalue problem
  out <- tryCatch(get_k_eigs(R.a,K,f.w),error=function(err){get_k_geigen(R.a,R.b,K,f.w)})
  
  return(out)
}

dim(evs$vectors)
length(evs$values)
N <- 256
N.fourier <- floor(N/2) + 1
freq <- seq(0,0.5, length.out = N.fourier)

tapers_matlist <- list()
e.vals_list <- list()

for(i in 1:length(freq)){
  print(i)
  temp <- get.weights_bronez(t.n = 1:256, K = 3, f.c = freq[i], f.w = 4/256)
  #tapers_matlist[[i]] <- temp$weights
  #e.vals_list[[i]] <- temp$eigenvalues
}

