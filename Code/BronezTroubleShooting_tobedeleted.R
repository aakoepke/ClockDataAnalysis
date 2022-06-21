###what's going on with the changing of the f center?

top.eight <- matrix(NA, nrow = length(f.c), ncol = 8)

for(i in 1:length(f.c)){
  R.a <- 1/(pi*(dist.mat))*(sin(dist.mat*f.w))
  
  R.a[row(R.a) == col(R.a)] <- f.w/(pi)
  
  #Solve the generalized eigenvalue problem
  evs <- geigen(R.a,R.b, symmetric = TRUE)
  lambdas <- sort(evs$values, decreasing = TRUE)
  
  top.eight[i,] <- lambdas[1:8]
  
}

matplot(f.c/(2*pi),top.eight) 
top.eight



####### Chave 2019 Paper #######

f.seq <- seq(-0.45,0.45,by = 0.05)
#R.a.prime 


R.a.prime <- (1/(pi*dist.mat))*sin(2*pi*dist.mat*0.05)

R.b.prime <- exp(-j*2*pi*f.seq[i]*dist.mat)*(1/(pi*dist.mat))*sin(2*pi*dist.mat)



