######### Multitaper Spectral Estimate Troubleshooting #########
library(waveslim) #dpss tapers

############# Recreate example from SAUTS ########

### R CODE FOR REPRODUCING CONTENT OF FIGURES AND TABLES IN CHAPTER 8 ...

ar2_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_1.txt")
ar4_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_1.txt")
ocean_wave <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_wave.txt")
ac_time_differences <- scan("http://faculty.washington.edu/dbp/sauts/Data/maser_deglitched.txt")
ac_fractional_frequencies <- diff(ac_time_differences)*100/6
ts_100 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ts_100.txt")
wavelength_spec <- read.table("http://faculty.washington.edu/dbp/sauts/Data/wavelength_spec.txt")

### functions used to compute content of figures in Chapter 8 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/amt_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H_bar.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/circular_shift.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/compute_slepian_eigenvalue.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/cosine_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/create_tapered_series.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dft.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/direct_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_lwe.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_mt.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_wosa.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_hanning_50_wosa.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_it_innov_var.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_OLS_pgram_mt.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/d_to_alpha.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/hanning_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/is_even.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/lag_windows.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/mt_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sim_ar_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sine_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/slepian_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/spec_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/thomson_chave_miller_jackknife.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/trapezoidal_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ts_to_lag_window_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/var_biased.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/var_power_law_alpha_mt.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/var_power_law_alpha_pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/width_e_R_mt.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/wosa_sdf_est.R")

### NOTE: to install the sapa library, uncomment the following three statements
###       and execute them:
###
### install.packages("devtools")
### devtools::install_github("wconstan/ifultools")
### devtools::install_github("wconstan/sapa")

library(sapa)

###

fig_taper_or_tapered_series <- function(the_taper,tag,ts=NULL,x_text=1024,main=" ")
{
  N <- length(the_taper)
  xs <- 0:(N-1)
  ys <- if(is.null(ts)) the_taper else the_taper*ts
  left_p <- is.null(ts)
  y_lim <- if(left_p) c(-0.1,0.1) else c(-6,6)*sqrt(0.002)
  plot(xs,ys,
       xlim=c(0,N),xaxs="i",xlab=expression(italic(t)),
       ylim=y_lim,yaxs="i",ylab=" ",
       typ="l",lwd=if(left_p) 1 else 0.5,axes=FALSE,
       main=main)
  if(left_p) abline(h=0,lwd=0.5)
  axis(1,at=seq(0,1024,1024))
  axis(1,at=seq(0,1024,256),label=FALSE,tcl=-0.25)
  axis(2,at=c(y_lim[1],0,y_lim[2]),label=FALSE)
  text(x_text,0.075,tag,pos=2)
  box(bty="l")
}

###

fig_spectral_window <- function(the_tapers,W=4,v_solid=W/N,main=" ")
{
  the_tapers <- as.matrix(the_tapers)
  N <- nrow(the_tapers)
  K <- ncol(the_tapers)
  the_sws <- vector("list",K)
  for(k in 1:K)
    the_sws[[k]] <- spec_window(the_tapers[,k],pad=2^8,fix_nulls_p=TRUE)
  the_sw <- rowMeans(matrix(unlist(lapply(the_sws,function(x) x$sw)),ncol=K))
  plot(the_sws[[1]]$freqs,dB(the_sw),
       xlim=c(0,0.01),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-80,40),yaxs="i",ylab="dB",
       typ="l",axes=FALSE,
       main=main)
  abline(v=c(v_solid,B_H_bar(the_tapers)/2),lwd=0.5,lty=c("solid","dashed"))
  axis(1,at=seq(0,0.01,0.01))
  axis(2,at=seq(-80,40,20),las=2)
  axis(2,at=seq(-80,40,10),label=FALSE,tcl=-0.25)
  box(bty="l")
}

###

fig_true_and_est_sdfs <- function(the_tapers,tag,ts=ar4_1,coeffs=ar4_coeffs,innov_var=ar4_innov_var,cc_center=c(0.1,-50),main=" ")
{
  mtse <- mt_sdf_est(ts,the_tapers,center=FALSE)
  plot(mtse$freqs,dB(mtse$sdfe),
       xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-80,40),yaxs="i",ylab="dB",
       typ="l",lwd=0.5,axes=FALSE,
       main=main)
  true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
  lines(true_sdf$freqs,dB(true_sdf$sdf),lwd=1.0)
  lines(rep(cc_center[1],2),cc_center[2]+c(mtse$cc$up,-mtse$cc$down),lwd=0.5)
  lines(cc_center[1]+c(-mtse$cc$width/2,mtse$cc$width/2),rep(cc_center[2],2),lwd=0.5)
  axis(1,at=seq(0,0.5,0.5))
  axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-80,40,20),las=2)
  axis(2,at=seq(-80,40,10),label=FALSE,tcl=-0.25)
  text(0.49,25,tag,pos=2)
  box(bty="l")
}

###

ar2_innov_var <- 1
ar2_coeffs    <- c(0.75,-0.5)

ar4_innov_var <- 0.002
ar4_coeffs    <- c(2.7607, -3.8106, 2.6535, -0.9238)

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES


### Figure 364, left-hand column

fig_spectral_window(tapers_all[,1],main="Figure 364")
fig_spectral_window(tapers_all[,1:2],main="Figure 364")
fig_spectral_window(tapers_all[,1:3],main="Figure 364")
fig_spectral_window(tapers_all[,1:4],main="Figure 364")

### Figure 364, right-hand column

fig_true_and_est_sdfs(tapers_all[,1],expression(italic(K == 1)),main="Figure 364")
fig_true_and_est_sdfs(tapers_all[,1:2],expression(italic(K == 2)),main="Figure 364")
fig_true_and_est_sdfs(tapers_all[,1:3],expression(italic(K == 3)),main="Figure 364")
fig_true_and_est_sdfs(tapers_all[,1:4],expression(italic(K == 4)),main="Figure 364")

### Figure 365, left-hand column

fig_spectral_window(tapers_all[,1:5],main="Figure 365")
fig_spectral_window(tapers_all[,1:6],main="Figure 365")
fig_spectral_window(tapers_all[,1:7],main="Figure 365")
fig_spectral_window(tapers_all[,1:8],main="Figure 365")

### Figure 365, right-hand column

fig_true_and_est_sdfs(tapers_all[,1:5],expression(italic(K == 5)),main="Figure 365")
fig_true_and_est_sdfs(tapers_all[,1:6],expression(italic(K == 6)),main="Figure 365")
fig_true_and_est_sdfs(tapers_all[,1:7],expression(italic(K == 7)),main="Figure 365")
fig_true_and_est_sdfs(tapers_all[,1:8],expression(italic(K == 8)),main="Figure 365")



##### our function


mtse_ar4 <- multitaper_est(X.t = ar4_1, W = 0.00390625, K = 5)

mtse_sauts <- mt_sdf_est(ar4_1,tapers_all[,1:5],center=FALSE)

plot(seq(0,0.5, length.out = length(mtse_ar4$spectrum)), 10*log10(mtse_ar4$spectrum))
lines(mtse_sauts$freqs, dB(mtse_sauts$sdfe), col = "red")


##they match! woohoo!


##let's try on one of the simulations of white noise
wn <- rnorm(2048)
mtse_wn <- multitaper_est(X.t = wn, W = 0.001953125, K = 5)
tapers_4096 <- dpss.taper(n = 4096, k = 5, nw = 4, nmax = 4096)
mtse_sauts_wn <- mt_sdf_est(wn,tapers_2048[,1:5],center=FALSE)
periodogram_wn <- spec.pgram(wn)

plot(seq(0,0.5, length.out = length(mtse_wn$spectrum)), 10*log10(mtse_wn$spectrum))
points(mtse_sauts_wn$freqs, dB(mtse_sauts_wn$sdfe), col = "red")


boxplot(10*log10(mtse_wn$spectrum), dB(mtse_sauts_wn$sdfe), dB(periodogram_wn$spec))

ours <- dons <- perio <- matrix(data = NA, nrow = 500, ncol = 2049)

for(i in 101:500){
  set.seed(i)
  print(i)
  wn <- rnorm(4096)
  mtse_wn <- multitaper_est(X.t = wn, W = 0.0009765625, K = 5)
  mtse_sauts_wn <- mt_sdf_est(wn,tapers_4096[,1:5],center=FALSE)
  periodogram_wn <- spec.pgram(wn)
  
  ours[i,] <- mtse_wn$spectrum
  dons[i,] <- mtse_sauts_wn$sdfe
  perio[i,] <- c(NA,periodogram_wn$spec)
}


oursdB <- dB(ours)
donsdB <- dB(dons)
perioddB <- dB(perio)

hist(ours[,1], breaks = 25, freq = FALSE)
lines(seq(0,3, length.out = 100), dchisq(seq(0,3, length.out = 100), df = 2)/2)
hist(dons[,1], breaks = 25)
hist(perio[,2], breaks = 25)
