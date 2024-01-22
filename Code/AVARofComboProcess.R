############# Combination Noise Process ##################


##generate noise process
N <- 10000
set.seed(2)

#white + pink noise
white_noise <- rnorm(N) 

generate_pink_noise <- function(length, fs) {
  pink_noise <- cumsum(rnorm(length))
  pink_noise <- pink_noise - mean(pink_noise)
  pink_noise <- pink_noise / max(abs(pink_noise))
  return(pink_noise)
}


pink_noise <- generate_pink_noise(length = N, fs = 1)
combo_process <- 3e-17*white_noise + 1e-17*pink_noise
acf(combo_process)
hist(combo_process)
acf(dat)

par(mfrow = c(3,1))

plot(white_noise)
plot(pink_noise)
plot(combo_process)


##calculate Allan Deviation

tmp <- getAvars(N,combined_noise_data$Frequency, taus = taus)
oavar_combo <- tmp$avarRes$overavars
length(oavar_combo)
length(taus)
plot(log10(taus), log10(sqrt(oavar_combo)))



#### Using Tara's data #####
plot(combined_noise_data)
combined_noise_data <- read_csv("Data/combined_noise_data.csv")

taus <-  c(2^(0:15))
taus
N <- length(combined_noise_data$Frequency)

tmp <- getAvars(N,combined_noise_data$Frequency, taus = taus)
oavar_combo <- tmp$avarRes$overavars
length(oavar_combo)
length(taus)
plot(log10(taus), log10(sqrt(oavar_combo)))

N <- 10000
white_noise <- rnorm(N) 
pink_noise <- TK95(N = N, alpha = 1)
combo_process <- white_noise + 0.01*pink_noise
taus <-  c(2^(0:8), seq(300,500, by = 20))
taus
tmp_wn <- getAvars(N,white_noise, taus = taus)
tmp_pink <- getAvars(N,0.01*pink_noise, taus = taus)
tmp_combo <- getAvars(N,combo_process, taus = taus)
  
par(mfrow = c(3,1))

plot(log10(taus), log10(sqrt(tmp_wn$avarRes$overavars)), type = "l")
lines(log10(taus), log10(sqrt(tmp_pink$avarRes$overavars)))
lines(log10(taus), log10(sqrt(tmp_pink$avarRes$overavars)) + log10(sqrt(tmp_wn$avarRes$overavars)), col = "blue")
