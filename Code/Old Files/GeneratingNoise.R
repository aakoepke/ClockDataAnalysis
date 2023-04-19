######Noise Generation####
library(RobPer)

#White Noise
y <- TK95(N = 2000, alpha = 0)
t <- seq(along = y)
plot(t,y,type = "l", main = "white noise", ylab = "y", xlab = "t")

#Flicker (Pink) Noise
y <- TK95(N = 2000, alpha = 1)
t <- seq(along = y)
plot(t,y,type = "l", main = "flicker noise", ylab = "y", xlab = "t")


#Random Walk
y<- TK95(N = 2000, alpha = 2)
t <- seq(along = y)
plot(t,y,type = "l", main = "random walk noise", ylab = "y", xlab = "t")
