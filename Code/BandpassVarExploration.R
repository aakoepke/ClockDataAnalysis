#### Band-Pass Variance #########

## calculating BPV estimates for known processes ##

#generate a process
N = 2048
y <- rnorm(N) #white noise

#calculate spectral estimate
periodogram.y <- spec.pgram(y)

#use to estimate BPV