###What's going on with the weighting sequences plots?

i = 21 #0.05 width, center 0, arithmetic sampling

R.a <- 1/(pi*(dist.mat))*(sin(dist.mat*f.w[i]))

R.a[row(R.a) == col(R.a)] <- f.w[i]/pi

#Solve the generalized eigenvalue problem
evs <- geigen(R.a,R.b, symmetric = TRUE)
lambdas <- sort(evs$values, decreasing = TRUE)


Conj(t(evs$vectors[,50]))%*%R.b%*%evs$vectors[,50] #=1
#so multiply the vector by sqrt(0.05/(pi))
evec.test <- sqrt(0.05/pi)*evs$vectors[,50]
plot(1:50,abs(evec.test))


