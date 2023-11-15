######## SAUTS Don Percival Code #########

wind_speed <- scan("http://faculty.washington.edu/dbp/sauts/Data/wind_speed_128.txt")
ar2_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_1.txt")
ar2_2 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_2.txt")
ar2_3 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_3.txt")
ar2_4 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_4.txt")
ar4_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_1.txt")
ar4_2 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_2.txt")
ar4_3 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_3.txt")
ar4_4 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_4.txt")
earth_20 <- scan("http://faculty.washington.edu/dbp/sauts/Data/earth_20.txt")
ocean_wave <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_wave.txt")
chaotic_beam <- scan("http://faculty.washington.edu/dbp/sauts/Data/chaotic_beam.txt")
ocean_noise <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_noise_128.txt")

### functions used to compute content of figures in Chapter 6 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/circular_shift.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/cosine_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/create_tapered_series.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dft.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/direct_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_DCTII.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_lag_window_sdf_estimator.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_shp.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_shp_squared.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/fejer_kernel.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/hanning_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/is_even.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/next_power_of_2.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/rectangular_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sim_ar_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/slepian_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/spec_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")

ar2_innov_var <- 1
ar2_coeffs    <- c(0.75,-0.5)

ar4_innov_var <- 0.002
ar4_coeffs    <- c(2.7607, -3.8106, 2.6535, -0.9238)


fig_173 <- function(ts,coeffs,innov_var,y_ats,tag)
{
  the_pgram <- pgram(ts,center=FALSE)
  plot(the_pgram$freqs,dB(the_pgram$sdfe),
       xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-60,20),yaxs="i",ylab=paste("AR(",length(coeffs),") spectra  (dB)",sep=""),
       typ="l",lwd=0.25,col="gray40",axes=FALSE)
  the_ar_spec <- ar_coeffs_to_sdf(coeffs, innov_var, N_pad=1024)
  lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
  if(length(coeffs) == 4)
  {
    N <- length(ts)
    temp <- ev_lag_window_sdf_estimator(ar_coeffs_to_acvs(coeffs,N-1,innov_var,FALSE))
    lines(temp$freqs, dB(temp$sdf_ev), lwd=0.5)
  }
  axis(1,at=seq(0,0.5,0.5))
  axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-60,20,20),las=2)
  axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
  text(x=0.5,y=10,tag,pos=2)
  box(bty="l")
}



fig_173(ar4_2,ar4_coeffs,ar4_innov_var,seq(0,150,50),"")


## functions used to compute content of figures in Chapter 2 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/sim_ar_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sim_arma_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")

###

generate_double_exponential <- function() -log(1-runif(1))/(if(runif(1) <= 0.5) sqrt(2) else -sqrt(2))

generate_discrete_rv <- function()
{
  rn <- runif(1)
  if(rn <= 0.02) -5 else if(rn <= 0.98) 0 else 5
}

generate_mismash_deviate <- function()
{
  rn <- runif(1)
  if(rn <= 0.25) rnorm(1) else if(rn <= 0.5) runif(1,-sqrt(3),sqrt(3)) else if(rn <= 0.75) generate_double_exponential() else generate_discrete_rv()
}
fig_34 <- function(ts,tag="(a)",two_or_four="2")
{
  xs <- 0:1023
  plot(xs,ts,
       xlim=c(0,1024),xlab=expression(italic(t)),
       ylim=c(-5,5),ylab=paste("AR(4) series",sep=""),
       typ="l",lwd=0.5,axes=FALSE)
  axis(1,at=seq(0,1024,512))
  axis(1,at=seq(0,1024,256),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-5,5,5),las=2)
  axis(2,at=seq(-5,5,1),label=FALSE,tcl=-0.25)
  text(1000,5,tag,pos=4)
  box(bty="l")
}

par(mfrow = c(1,2), mar = c(4,5,1,4))

fig_34(ar4_1,"","")
fig_173(ar4_1,ar4_coeffs,ar4_innov_var,seq(0,150,50),"")



fig_185 <- function(ys,big_y_ats=seq(-5,5,5),little_y_ats=NULL,y_lab="AR(4) series")
{
  N <- length(ys)
  plot(0:(N-1),ys,
       xlim=c(0,N),xlab=expression(italic(t)),
       ylim=c(big_y_ats[1],big_y_ats[length(big_y_ats)]),ylab=y_lab,
       typ="l",lwd=0.25,axes=FALSE)
  axis(1,at=seq(0,1024,512))
  axis(1,at=seq(0,1024,256),label=FALSE,tcl=-0.25)
  axis(2,at=big_y_ats,las=2)
  axis(2,at=little_y_ats,label=FALSE,tcl=-0.25)
  box(bty="l")
}

the_taper <- hanning_taper(1024)

par(mfrow = c(3,1))
fig_185(ar4_1,little=seq(-5,5,1), y_lab = expression(X[t]))
fig_185(the_taper,big=seq(0,0.06,0.02),y_lab= expression(h[t]))
fig_185(the_taper*ar4_1,big=seq(-0.2,0.2,0.1),y_lab=expression(W[t]))



#Chapter 8

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
install.packages("devtools")
devtools::install_github("wconstan/ifultools")
devtools::install_github("wconstan/sapa")

library(sapa)

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
  text(x_text,0.075,tag,pos=2, cex = 1.8)
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
  text(0.49,25,tag,pos=2, cex = 2)
  box(bty="l")
}

###

ar2_innov_var <- 1
ar2_coeffs    <- c(0.75,-0.5)

ar4_innov_var <- 0.002
ar4_coeffs    <- c(2.7607, -3.8106, 2.6535, -0.9238)


tapers_all <- t(as.matrix(taper("dpss",1024,n.taper=8,bandwidth=4)))


par(mfrow = c(4,2))
### Figure 360, left-hand column

fig_taper_or_tapered_series(tapers_all[,1],expression(italic(k == 0)),main="")
fig_taper_or_tapered_series(tapers_all[,1]," ",ar4_1,main = "")
fig_taper_or_tapered_series(tapers_all[,2],expression(italic(k == 1)),main="")
fig_taper_or_tapered_series(tapers_all[,2]," ",ar4_1,main="")
fig_taper_or_tapered_series(tapers_all[,3],expression(italic(k == 2)),main="")
fig_taper_or_tapered_series(tapers_all[,3]," ",ar4_1,main="")
fig_taper_or_tapered_series(tapers_all[,4],expression(italic(k == 3)),main="")
fig_taper_or_tapered_series(tapers_all[,4]," ",ar4_1,main="")

par(mfrow = c(1,1))
fig_true_and_est_sdfs(tapers_all[,1:4],expression(italic(K == 4)),main="")


##Chapter 1: time series vs. sdf

source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sim_ar_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")



fig_129a <- function(the_sdf,tag)
{
  freqs <- seq(0,0.5,length=length(the_sdf)) 
  plot(freqs,the_sdf,
       xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
       ylim=c(0,10),yaxs="i",ylab=expression(S(italic(f))),
       typ="l",axes=FALSE)
  axis(1,at=seq(0,0.5,0.5))
  axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(0,10,5),las=2)
  axis(2,at=seq(0,10,1),label=FALSE,tcl=-0.25)
  #text(0.5,9,tag,pos=2)
  box(bty="l")
}

ar1_innov_var <- 0.36
ar1_coeffs    <- 0.8
ar1_sdf <- ar_coeffs_to_sdf(ar1_coeffs,ar1_innov_var,N_pad=1024)$sdf

ar2_innov_var <- 0.48
ar2_coeffs    <- c(-0.8,-0.6)
ar2_sdf <- ar_coeffs_to_sdf(ar2_coeffs,ar2_innov_var,N_pad=1024)$sdf

wn_1_sdf <- rep(3,2)

wn_2_sdf <- rep(8,2)

### Figure 129a, top row

fig_129a(ar1_sdf,"(a)")

fig_129a(ar2_sdf,"(b)")

### Figure 129a, bottom row

fig_129a(wn_1_sdf,"(c)")

fig_129a(wn_2_sdf,"(d)")



## realization of AR(2) X_t = 0.8X_{t-1} + e_t, sigma_e^2 = 0.36

ar1_sim <- arima.sim(n = 100, list(ar = 0.8), sd = 0.6)
plot(ar1_sim, ylim = c(-10,10), xlab = "t", ylab = expression(X[t]))





####################################################################################
####################################################################################
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
install.packages("devtools")
devtools::install_github("wconstan/ifultools")
devtools::install_github("wconstan/sapa")
library(devtools)
library(ifultools)
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

### Figures 360, 361, 362, 363, 364, 365

tapers_all <- t(as.matrix(taper("dpss",1024,n.taper=8,bandwidth=4)))

### Figure 360, left-hand column

fig_taper_or_tapered_series(tapers_all[,1],expression(italic(k == 0)),main="Figure 360")
fig_taper_or_tapered_series(tapers_all[,2],expression(italic(k == 1)),main="Figure 360")
fig_taper_or_tapered_series(tapers_all[,3],expression(italic(k == 2)),main="Figure 360")
fig_taper_or_tapered_series(tapers_all[,4],expression(italic(k == 3)),main="Figure 360")

### Figure 360, right-hand column

fig_taper_or_tapered_series(tapers_all[,1]," ",ar4_1,main="Figure 360")
fig_taper_or_tapered_series(tapers_all[,2]," ",ar4_1,main="Figure 360")
fig_taper_or_tapered_series(tapers_all[,3]," ",ar4_1,main="Figure 360")
fig_taper_or_tapered_series(tapers_all[,4]," ",ar4_1,main="Figure 360")

### Figure 361, left-hand column

fig_spectral_window(tapers_all[,1],main="Figure 361")
fig_spectral_window(tapers_all[,2],main="Figure 361")
fig_spectral_window(tapers_all[,3],main="Figure 361")
fig_spectral_window(tapers_all[,4],main="Figure 361")

### Figure 361, right-hand column

fig_true_and_est_sdfs(tapers_all[,1],expression(italic(k == 0)),main="Figure 361")
fig_true_and_est_sdfs(tapers_all[,2],expression(italic(k == 1)),main="Figure 361")
fig_true_and_est_sdfs(tapers_all[,3],expression(italic(k == 2)),main="Figure 361")
fig_true_and_est_sdfs(tapers_all[,4],expression(italic(k == 3)),main="Figure 361")

### Figure 362, left-hand column

fig_taper_or_tapered_series(tapers_all[,5],expression(italic(k == 4)),main="Figure 362")
fig_taper_or_tapered_series(tapers_all[,6],expression(italic(k == 5)),main="Figure 362")
fig_taper_or_tapered_series(tapers_all[,7],expression(italic(k == 6)),main="Figure 362")
fig_taper_or_tapered_series(tapers_all[,8],expression(italic(k == 7)),main="Figure 362")

### Figure 362, right-hand column

fig_taper_or_tapered_series(tapers_all[,5]," ",ar4_1,main="Figure 362")
fig_taper_or_tapered_series(tapers_all[,6]," ",ar4_1,main="Figure 362")
fig_taper_or_tapered_series(tapers_all[,7]," ",ar4_1,main="Figure 362")
fig_taper_or_tapered_series(tapers_all[,8]," ",ar4_1,main="Figure 362")

### Figure 363, left-hand column

fig_spectral_window(tapers_all[,5],main="Figure 363")
fig_spectral_window(tapers_all[,6],main="Figure 363")
fig_spectral_window(tapers_all[,7],main="Figure 363")
fig_spectral_window(tapers_all[,8],main="Figure 363")

### Figure 363, right-hand column

fig_true_and_est_sdfs(tapers_all[,5],expression(italic(k == 4)),main="Figure 363")
fig_true_and_est_sdfs(tapers_all[,6],expression(italic(k == 5)),main="Figure 363")
fig_true_and_est_sdfs(tapers_all[,7],expression(italic(k == 6)),main="Figure 363")
fig_true_and_est_sdfs(tapers_all[,8],expression(italic(k == 7)),main="Figure 363")

### Figure 364, left-hand column

fig_spectral_window(tapers_all[,1],main="Figure 364")
fig_spectral_window(tapers_all[,1:2],main="Figure 364")
fig_spectral_window(tapers_all[,1:3],main="Figure 364")
fig_spectral_window(tapers_all[,1:4],main="Figure 364")

### Figure 364, right-hand column

par(mfrow = c(2,4))
fig_true_and_est_sdfs(tapers_all[,1],expression(italic(K == 1)),main="")
fig_true_and_est_sdfs(tapers_all[,1:2],expression(italic(K == 2)),main="")
fig_true_and_est_sdfs(tapers_all[,1:3],expression(italic(K == 3)),main="")
fig_true_and_est_sdfs(tapers_all[,1:4],expression(italic(K == 4)),main="")

### Figure 365, left-hand column

# fig_spectral_window(tapers_all[,1:5],main="Figure 365")
# fig_spectral_window(tapers_all[,1:6],main="Figure 365")
# fig_spectral_window(tapers_all[,1:7],main="Figure 365")
# fig_spectral_window(tapers_all[,1:8],main="Figure 365")

### Figure 365, right-hand column

fig_true_and_est_sdfs(tapers_all[,1:5],expression(italic(K == 5)),main="")
fig_true_and_est_sdfs(tapers_all[,1:6],expression(italic(K == 6)),main="")
fig_true_and_est_sdfs(tapers_all[,1:7],expression(italic(K == 7)),main="")
fig_true_and_est_sdfs(tapers_all[,1:8],expression(italic(K == 8)),main="")

### Figures 367 and 368 ###

fig_367 <- function(the_sq_tapers,tag,main="Figure 367")
{
  ys <- rowSums(as.matrix(the_sq_tapers))/0.0146
  N <- length(ys)
  xs <- 0:(N-1)
  plot(xs,ys,
       xlim=c(0,N),xaxs="i",xlab=expression(italic(t)),
       ylim=c(0,1.01),yaxs="i",ylab=" ",
       typ="l",axes=FALSE,
       main=main)
  axis(1,at=seq(0,1024,1024))
  axis(1,at=seq(0,1024,256),label=FALSE,tcl=-0.25)
  axis(2,at=c(0,1),label=FALSE)
  axis(2,at=seq(0,1,0.1),label=FALSE,tcl=-0.25)
  text(512,0.9,tag,pos=1)
  box(bty="l")
}

tapers_sq_all <- t(as.matrix(taper("dpss",1024,n.taper=9,bandwidth=4))^2)

### Figure 367, left-hand column

fig_367(tapers_sq_all[,1],expression(italic(K == 1)))
fig_367(tapers_sq_all[,1:2],expression(italic(K == 2)))
fig_367(tapers_sq_all[,1:3],expression(italic(K == 3)))
fig_367(tapers_sq_all[,1:4],expression(italic(K == 4)))

### Figure 367, right-hand column

fig_367(tapers_sq_all[,1:5],expression(italic(K == 5)))
fig_367(tapers_sq_all[,1:6],expression(italic(K == 6)))
fig_367(tapers_sq_all[,1:7],expression(italic(K == 7)))
fig_367(tapers_sq_all[,1:8],expression(italic(K == 8)))

### Figure 368, top row of plots

tapers_all <- t(as.matrix(taper("dpss",1024,n.taper=9,bandwidth=4)))

fig_taper_or_tapered_series(tapers_all[,9],expression(italic(k == 8)),x_text=550,main="Figure 368")
fig_taper_or_tapered_series(tapers_all[,9]," ",ar4_1,main="Figure 368")

### Figure 368, middle row of plots

fig_spectral_window(tapers_all[,9],main="Figure 368")
fig_true_and_est_sdfs(tapers_all[,9]," ",main="Figure 368")

### Figure 368, bottom plot

fig_367(tapers_sq_all[,1:9],expression(italic(K == 9)),main="Figure 368")

### Table 369 ###

tab_369 <- function(NW,N=1024)
{
  tapers_all <- t(as.matrix(taper("dpss",N,n.taper=2*NW,bandwidth=NW)))
  return(round(sapply(1:(2*NW),function(K) B_H_bar(tapers_all[,1:K])/width_e_R_mt(tapers_all[,1:K])),2))
}

### Table 369, rows 1 to 4

tab_369(1)
tab_369(2)
tab_369(4)
tab_369(8)

### Figure 370 ###

fig_370 <- function(ts,coeffs,innov_var,jackknife_p=FALSE)
{
  ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
  plot(ar_spec$freqs,dB(ar_spec$sdf),
       xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-67,24),yaxs="i",ylab="dB",
       typ="l",axes=FALSE,
       main="Figure 370")
  mtse <- mt_sdf_est(ts,t(as.matrix(taper("dpss",length(ts),n.taper=5,bandwidth=4))),center=FALSE)
  if(jackknife_p)
  {
    temp <- thomson_chave_miller_jackknife(mtse)
    lines(mtse$freqs, dB(temp$tc_lci), lwd=0.5)
    lines(mtse$freqs, dB(temp$tc_uci), lwd=0.5)
  }
  else
  {
    lines(mtse$freqs,dB(mtse$sdfe)-mtse$cc$down,lwd=0.5)
    lines(mtse$freqs,dB(mtse$sdfe)+mtse$cc$up,lwd=0.5)
    cc_center <- c(0.05,-40)
    lines(rep(cc_center[1],2),cc_center[2]+c(mtse$cc$up,-mtse$cc$down),lwd=0.5)
    lines(cc_center[1]+c(-mtse$cc$width/2,mtse$cc$width/2),rep(cc_center[2],2),lwd=0.5)
    abline(v=c(mtse$cc$width/2,0.5-mtse$cc$width/2),lty="dashed")
  }
  axis(1,at=seq(0,0.5,0.5))
  axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-60,20,20),las=2)
  axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
  box(bty="l")
}

### Figure 370, top plot

fig_370(ar4_1,ar4_coeffs,ar4_innov_var)

### Figure 370, bottom plot

fig_370(ar4_1,ar4_coeffs,ar4_innov_var,jackknife_p=TRUE)

### Figures 373a, 373b and 403

fig_373a <- function(which_tapers="dpss",bandwidth=4,N=64,main=" ")
{
  xs <- 1:N
  tapers_all <- t(as.matrix(taper(which_tapers,N,n.taper=N,bandwidth=bandwidth)))
  covs_N_by_N <- matrix(nrow=N,ncol=N)
  for(i in 1:N) covs_N_by_N[i,i:N] <- covs_N_by_N[i:N,i] <- sapply(i:N,function(j){sum(tapers_all[,i]*tapers_all[,j]*(-1)^(1:N))^2})
  diag(covs_N_by_N) <- diag(covs_N_by_N) + 1
  ys <- dB(sapply(1:N,function(k) sum(covs_N_by_N[1:k,1:k])/k^2))
  plot(xs,ys,
       xlim=c(0,N),xaxs="i",xlab=expression(italic(K)),
       ylim=c(-20,0),ylab="variance",
       typ="l",lwd=2,axes=FALSE,
       main=main)
  lines(xs,dB(1/(1:N)), lwd=1)
  abline(h=dB(2/N),lty="dotted")
  if(which_tapers == "dpss") abline(v=2*bandwidth,lty="dotted")
  axis(1,at=seq(0,N,16))
  axis(1,at=seq(0,N,4),label=FALSE,tcl=-0.25)
  axis(2, at=seq(-20,0,10),labels=c(0.01,0.1,1),las=2)
  axis(2, at=dB(c((1:10)/1000,(1:10)/100,(1:10)/10,1:10)), labels=FALSE, tcl=-0.25)
  box(bty="l")
}

### Figure 373a

fig_373a(main="Figure 373a")

### Figure 373b

fig_373a(bandwidth=16,main="Figure 373b")

### Figure 403

fig_373a("sine",main="Figure 403")

### Figure 393 ###

fig_393 <- function(N,k,tag,y_text=y_lim[1]+diff(y_lim)*0.9)
{
  xs <- 0:(N-1)
  y_lim <- c(if(k==0) 0 else -0.25,0.25)
  plot(xs,sine_taper(N,k),
       xlim=c(0,N-1),xlab=expression(italic(t)),
       ylim=y_lim,ylab=" ",
       typ="l",lty="dashed",col="red",axes=FALSE,
       main="Figure 393")
  gen_Sigma_pl <- function(N)
  {
    pl_acvs <- sapply(1:(N-1),function(tau) if(tau==0) 1/12 else (-1)^tau/(2*pi^2*tau^2))
    Sigma_pl <- diag(rep(1/12,N))
    for(tau in 1:(N-1))
    {
      for(off in 0:(N-1-tau))
      {
        Sigma_pl[1+off,tau+1+off] <- pl_acvs[tau]
        Sigma_pl[tau+1+off,1+off] <- pl_acvs[tau]
      }
    }
    return(Sigma_pl)
  }
  Sigma_pl <- gen_Sigma_pl(N)
  min_bias_taper <- eigen(Sigma_pl,TRUE)$vectors[,N-k]
  if(min_bias_taper[1] < 0) min_bias_taper <- -min_bias_taper
  lines(xs,min_bias_taper)
  abline(h=0,lwd=0.5)
  axis(1,at=range(xs))
  axis(2,at=c(round(y_lim,1),0),las=2)
  axis(2,at=seq(-0.2,0.2,0.1),label=FALSE,tcl=-0.25)
  text(N,y_text,tag,pos=2)
  box(bty="l")
}

### Figure 393, left-hand column (but with dashed lines colored red)

fig_393(32,0,expression(italic(k == 0)))
fig_393(32,1,expression(italic(k == 1)))
fig_393(32,2,expression(italic(k == 2)),y_text=-0.15)
fig_393(32,3,expression(italic(k == 3)))

### Figure 393, right-hand column (but with dashed lines colored red)

fig_393(99,0,expression(italic(k == 0)))
fig_393(99,1,expression(italic(k == 1)))
fig_393(99,2,expression(italic(k == 2)))
fig_393(99,3,expression(italic(k == 3)))

### Figures 394, 395, 396, 397, 398 and 399

tapers_all <- t(as.matrix(taper("sine",1024,n.taper=8)))

### Figure 394, left-hand column

fig_taper_or_tapered_series(tapers_all[,1],expression(italic(k == 0)),main="Figure 394")
fig_taper_or_tapered_series(tapers_all[,2],expression(italic(k == 1)),main="Figure 394")
fig_taper_or_tapered_series(tapers_all[,3],expression(italic(k == 2)),main="Figure 394")
fig_taper_or_tapered_series(tapers_all[,4],expression(italic(k == 3)),main="Figure 394")

### Figure 394, right-hand column

fig_taper_or_tapered_series(tapers_all[,1]," ",ar4_1,main="Figure 394")
fig_taper_or_tapered_series(tapers_all[,2]," ",ar4_1,main="Figure 394")
fig_taper_or_tapered_series(tapers_all[,3]," ",ar4_1,main="Figure 394")
fig_taper_or_tapered_series(tapers_all[,4]," ",ar4_1,main="Figure 394")

### Figure 395, left-hand column

fig_spectral_window(tapers_all[,1],v_solid=NA,main="Figure 395")
fig_spectral_window(tapers_all[,2],v_solid=NA,main="Figure 395")
fig_spectral_window(tapers_all[,3],v_solid=NA,main="Figure 395")
fig_spectral_window(tapers_all[,4],v_solid=NA,main="Figure 395")

### Figure 395, right-hand column

fig_true_and_est_sdfs(tapers_all[,1],expression(italic(k == 0)),main="Figure 395")
fig_true_and_est_sdfs(tapers_all[,2],expression(italic(k == 1)),main="Figure 395")
fig_true_and_est_sdfs(tapers_all[,3],expression(italic(k == 2)),main="Figure 395")
fig_true_and_est_sdfs(tapers_all[,4],expression(italic(k == 3)),main="Figure 395")

### Figure 396, left-hand column

fig_taper_or_tapered_series(tapers_all[,5],expression(italic(k == 4)),main="Figure 396")
fig_taper_or_tapered_series(tapers_all[,6],expression(italic(k == 5)),main="Figure 396")
fig_taper_or_tapered_series(tapers_all[,7],expression(italic(k == 6)),main="Figure 396")
fig_taper_or_tapered_series(tapers_all[,8],expression(italic(k == 7)),main="Figure 396")

### Figure 396, right-hand column

fig_taper_or_tapered_series(tapers_all[,5]," ",ar4_1,main="Figure 396")
fig_taper_or_tapered_series(tapers_all[,6]," ",ar4_1,main="Figure 396")
fig_taper_or_tapered_series(tapers_all[,7]," ",ar4_1,main="Figure 396")
fig_taper_or_tapered_series(tapers_all[,8]," ",ar4_1,main="Figure 396")

### Figure 397, left-hand column

fig_spectral_window(tapers_all[,5],v_solid=NA,main="Figure 397")
fig_spectral_window(tapers_all[,6],v_solid=NA,main="Figure 397")
fig_spectral_window(tapers_all[,7],v_solid=NA,main="Figure 397")
fig_spectral_window(tapers_all[,8],v_solid=NA,main="Figure 397")

### Figure 397, right-hand column

fig_true_and_est_sdfs(tapers_all[,5],expression(italic(k == 4)),main="Figure 397")
fig_true_and_est_sdfs(tapers_all[,6],expression(italic(k == 5)),main="Figure 397")
fig_true_and_est_sdfs(tapers_all[,7],expression(italic(k == 6)),main="Figure 397")
fig_true_and_est_sdfs(tapers_all[,8],expression(italic(k == 7)),main="Figure 397")

### Figure 398, left-hand column

fig_spectral_window(tapers_all[,1],v_solid=2/(2*1025),main="Figure 398")
fig_spectral_window(tapers_all[,1:2],v_solid=3/(2*1025),main="Figure 398")
fig_spectral_window(tapers_all[,1:3],v_solid=4/(2*1025),main="Figure 398")
fig_spectral_window(tapers_all[,1:4],v_solid=5/(2*1025),main="Figure 398")

### Figure 398, right-hand column

fig_true_and_est_sdfs(tapers_all[,1],expression(italic(K == 1)),main="Figure 398")
fig_true_and_est_sdfs(tapers_all[,1:2],expression(italic(K == 2)),main="Figure 398")
fig_true_and_est_sdfs(tapers_all[,1:3],expression(italic(K == 3)),main="Figure 398")
fig_true_and_est_sdfs(tapers_all[,1:4],expression(italic(K == 4)),main="Figure 398")

### Figure 399, left-hand column

fig_spectral_window(tapers_all[,1:5],v_solid=6/(2*1025),main="Figure 399")
fig_spectral_window(tapers_all[,1:6],v_solid=7/(2*1025),main="Figure 399")
fig_spectral_window(tapers_all[,1:7],v_solid=8/(2*1025),main="Figure 399")
fig_spectral_window(tapers_all[,1:8],v_solid=9/(2*1025),main="Figure 399")

### Figure 399, right-hand column

fig_true_and_est_sdfs(tapers_all[,1:5],expression(italic(K == 5)),main="Figure 399")
fig_true_and_est_sdfs(tapers_all[,1:6],expression(italic(K == 6)),main="Figure 399")
fig_true_and_est_sdfs(tapers_all[,1:7],expression(italic(K == 7)),main="Figure 399")
fig_true_and_est_sdfs(tapers_all[,1:8],expression(italic(K == 8)),main="Figure 399")

### Figure 401 ###

fig_401_upper_plots <- function(k,tag,N=1024,pot=18)
{
  a_taper <- sine_taper(N,k)
  n_freqs <- ceiling(0.01*2^pot) + 1
  freqs <- (0:(n_freqs-1))/2^pot
  freqs_two_sided <- c(-rev(freqs),freqs[-1])
  sw <- (abs(dft(c(a_taper,rep(0,2^pot-N))))^2)[1:n_freqs]
  patch_me <- function(x)
  {
    x[x[1:(N-2)] > x[2:(N-1)] & x[2:(N-1)] < x[3:N]] <- -100
    return(x)
  }
  sw_dB <- patch_me(careful_dB(sw))
  sw_two_sided_dB <- c(rev(sw_dB),sw_dB[-1])
  plot(freqs_two_sided,sw_two_sided_dB,
       xlim=c(-0.01,0.01),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-20,40),yaxs="i",ylab="dB",
       typ="l",axes=FALSE,
       main="Figure 401")
  c_freq <- if(k==0) 0 else (k+1)/(2*(N+1))
  max_sw_dB <- max(sw_dB)
  offset <- 8.0
  diddle <- 1.5
  lines(rep(c_freq-1/(2*(N+1)),2),c(max_sw_dB+offset-diddle,max_sw_dB+offset+diddle))
  lines(rep(c_freq+1/(2*(N+1)),2),c(max_sw_dB+offset-diddle,max_sw_dB+offset+diddle))
  lines(c(c_freq-1/(2*(N+1)),c_freq+1/(2*(N+1))),rep(max_sw_dB+offset,2))
  if(k > 0)
  {
    lines(rep(-c_freq-1/(2*(N+1)),2),c(max_sw_dB+offset-diddle,max_sw_dB+offset+diddle))
    lines(rep(-c_freq+1/(2*(N+1)),2),c(max_sw_dB+offset-diddle,max_sw_dB+offset+diddle))
    lines(c(-c_freq-1/(2*(N+1)),-c_freq+1/(2*(N+1))),rep(max_sw_dB+offset,2))
  }
  dB_2 <- dB(2)
  upper_freqs <- freqs > c_freq
  i_min <- which.min(abs(sw_dB[upper_freqs] - (max_sw_dB-dB_2)))
  f_right_3dB_down <- (freqs[upper_freqs])[i_min]
  if(k==0) f_left_3dB_down <- - f_right_3dB_down
  else
  {
    lower_freqs <- freqs < c_freq
    i_min <- which.min(abs(sw_dB[lower_freqs] - (max_sw_dB-dB_2)))
    f_left_3dB_down <- (freqs[lower_freqs])[i_min]
  }
  lines(c(f_left_3dB_down,f_right_3dB_down),rep(max_sw_dB-dB_2,2),lwd=0.5)
  if(k>0) lines(c(-f_left_3dB_down,-f_right_3dB_down),rep(max_sw_dB-dB_2,2),lwd=0.5)
  axis(1,at=seq(-0.01,0.01,0.01))
  axis(1,at=seq(-0.01,0.01,0.001),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-20,40,20),las=2)
  axis(2,at=seq(-20,40,10),label=FALSE,tcl=-0.25)
  text(0.0096,32.5,tag,pos=2)
  box(bty="l")
}

fig_401_bottom_plot <- function(K=8,N=1024,v_dB=34,delta_dB=7,diddle=1.5)
{
  plot(rep(-1/(2*(N+1)),2),c(v_dB-diddle,v_dB+diddle),
       xlim=c(-0.01,0.01),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-20,40),yaxs="i",ylab=" ",
       typ="l",axes=FALSE,
       main="Figure 401")
  lines(rep(1/(2*(N+1)),2),c(v_dB-diddle,v_dB+diddle))
  lines(c(-1/(2*(N+1)),1/(2*(N+1))),rep(v_dB,2))
  abline(v=((1+1)/(2*(N+1))) * c(-1,1),lwd=0.5,lty="solid")
  abline(v=((K+1)/(2*(N+1))) * c(-1,1),lwd=0.5,lty="solid")
  abline(v=((1+0.5)/(2*(N+1))) * c(-1,1),lwd=1.0,lty="dotted")
  abline(v=((K+0.5)/(2*(N+1))) * c(-1,1),lwd=1.0,lty="dotted")
  for(k in 1:(K-1))
  {
    v_dB <- v_dB - delta_dB
    c_freq <- (k+1)/(2*(N+1))
    lines(rep(c_freq-1/(2*(N+1)),2),c(v_dB-diddle,v_dB+diddle))
    lines(rep(c_freq+1/(2*(N+1)),2),c(v_dB-diddle,v_dB+diddle))
    lines(c(c_freq-1/(2*(N+1)),c_freq+1/(2*(N+1))),rep(v_dB,2))
    lines(rep(-c_freq-1/(2*(N+1)),2),c(v_dB-diddle,v_dB+diddle))
    lines(rep(-c_freq+1/(2*(N+1)),2),c(v_dB-diddle,v_dB+diddle))
    lines(c(-c_freq-1/(2*(N+1)),-c_freq+1/(2*(N+1))),rep(v_dB,2))
  }
  axis(1,at=seq(-0.01,0.01,0.01))
  axis(1,at=seq(-0.01,0.01,0.001),label=FALSE,tcl=-0.25)
  axis(2,at=c(-20,40),label=FALSE)
  box(bty="l")
}

### Figure 401, top plot

fig_401_upper_plots(0,expression(italic(k == 0)))

### Figure 401, 2nd plot

fig_401_upper_plots(1,expression(italic(k == 1)))

### Figure 401, 3rd plot

fig_401_upper_plots(7,expression(italic(k == 7)))

### Figure 401, bottom plot

fig_401_bottom_plot()

### GROUP 3: 21 to 30

### Table 402 ###

tab_402 <- function(K=16,N=1024)
{
  tapers_all <- t(as.matrix(taper("sine",N,n.taper=K)))
  return(round(sapply(1:K,function(K) B_H_bar(tapers_all[,1:K])/width_e_R_mt(tapers_all[,1:K])),2))
}

### Table 402

tab_402()

### Figure 406 ###

fig_406 <- function(dj,pn,hn,mt,innov_var,tag,L_K_max=25)
{
  xs <- 2:L_K_max
  ys_mt <- log10(sqrt(sapply(1:(L_K_max-1),function(i) mean((mt[,i] - innov_var)^2)))/innov_var)
  plot(xs,ys_mt,
       xlim=c(1,L_K_max+1),xaxs="i",xlab=expression(paste(italic(L)," or ",italic(K))),
       ylim=c(-1.35,1.2),ylab="normalized RMSEs",
       typ="p",pch=1,axes=FALSE,
       main="Figure 406")
  i <- which.min(ys_mt)
  points(xs[i],ys_mt[i],pch=16)
  ys_hn <- log10(sqrt(sapply(1:(L_K_max-1),function(i) mean((hn[,i] - innov_var)^2)))/innov_var)
  points(xs,ys_hn,pch=2)
  abline(h=log10(sqrt(mean((dj - innov_var)^2))/innov_var),lty="dashed")
  abline(h=log10(sqrt(mean((pn - innov_var)^2))/innov_var),lty="dotted")
  axis(1,at=seq(5,25,5))
  axis(1,at=seq(2,25,1),label=FALSE,tcl=-0.25)
  axis(2, at=(-4):2, labels=expression(10^{-4},10^{-3},10^{-2},10^{-1},10^0,10^1,10^2),las=2)
  minors <- log10(c((2:9)/10000,(2:9)/1000,(2:9)/100,(2:9)/10,2:9,(2:9)*10))
  axis(2, at=minors, labels=FALSE, tcl=-0.25)
  text(L_K_max-1,log10(14),tag,pos=2)
  box(bty="l")
}

### NOTE: evaluation of the following R code is time consuming
###       (for the AR(2) process, it took about 3 hours on a 2017-vintage
###       MacBook Pro to create iv_dj_ar2, iv_hn_ar2, iv_pn_ar2 and
###       iv_mt_ar2; a similar amount of time was needed to create
###       the corresponding settings to the AR(4) process)
###
###         N <- 1024
###         L_or_K_max <- 64
###         N_rep <- 10000
###         iv_dj_ar2  <- rep(NA,N_rep)
###         iv_hn_ar2 <- matrix(nrow=N_rep,ncol=L_or_K_max-1) 
###         iv_pn_ar2 <- matrix(nrow=N_rep,ncol=3)
###         iv_mt_ar2 <- matrix(nrow=N_rep,ncol=L_or_K_max-1) 
###         my_tapers <- t(as.matrix(taper("sine",N,n.taper=L_or_K_max)))
###         LD_stuff_ar2 <- step_down_LD_recursions(ar2_coeffs, var=ar2_innov_var, proc=FALSE)
###         set.seed(42)
###         for(i in 1:N_rep)
###         {
###             an_ar2_ts <- sim_ar_process(N,LD_stuff_ar2)
###             temp <- do_it_innov_var(an_ar2_ts,the_tapers=my_tapers,Ls=2:L_or_K_max,Ks=2:L_or_K_max)
###             iv_dj_ar2[i]  <- temp$dj
###             iv_hn_ar2[i,] <- temp$hn
###             iv_pn_ar2[i,] <- temp$pn
###             iv_mt_ar2[i,] <- temp$mt
###         }
###         
###         iv_dj_ar4  <- rep(NA,N_rep)
###         iv_hn_ar4 <- matrix(nrow=N_rep,ncol=L_or_K_max-1) 
###         iv_pn_ar4 <- matrix(nrow=N_rep,ncol=3)
###         iv_mt_ar4 <- matrix(nrow=N_rep,ncol=L_or_K_max-1) 
###         my_tapers <- t(as.matrix(taper("sine",N,n.taper=L_or_K_max))
###         LD_stuff_ar4 <- step_down_LD_recursions(ar4_coeffs, var=ar4_innov_var, proc=FALSE)
###         set.seed(42)
###         for(i in 1:N_rep)
###         {
###             an_ar4_ts <- sim_ar_process(N,LD_stuff_ar4)
###             temp <- do_it_innov_var(an_ar4_ts,the_tapers=my_tapers,Ls=2:L_or_K_max,Ks=2:L_or_K_max)
###             iv_dj_ar4[i]  <- temp$dj
###             iv_hn_ar4[i,] <- temp$hn
###             iv_pn_ar4[i,] <- temp$pn
###             iv_mt_ar4[i,] <- temp$mt
###         }
###
###       Evaluation of the following eight load forms alleviates
###       having to recreate iv_dj_ar2 etc.

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/iv_dj_ar2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/iv_hn_ar2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/iv_pn_ar2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/iv_mt_ar2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/iv_dj_ar4.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/iv_hn_ar4.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/iv_pn_ar4.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/iv_mt_ar4.Rdata"))

### Figure 406, left-hand plot

fig_406(iv_dj_ar2,iv_pn_ar2,iv_hn_ar2,iv_mt_ar2,ar2_innov_var,"AR(2)")

### Figure 406, right-hand plot

fig_406(iv_dj_ar4,iv_pn_ar4,iv_hn_ar4,iv_mt_ar4,ar4_innov_var,"AR(4)")

### Figure 408 ###

fig_408 <- function()
{
  fd_sdf <- function(f,delta) 1/(4*(sin(pi*f)^2))^delta
  freqs <- seq(0.001,0.5,0.001)
  xs <- log10(freqs)
  ys <- cbind(sapply(freqs,fd_sdf,0.25),
              sapply(freqs,fd_sdf,0.4),
              sapply(freqs,fd_sdf,0.45),
              sapply(freqs,fd_sdf,0.49))
  temp <- range(log10(range(ys)))
  plot(xs,log10(ys[,4]),
       xaxs="i",xlab=expression(italic(f)),
       ylim=c(floor(temp[1])+0.5,ceiling(temp[2])-0.5),yaxs="i",ylab=expression(italic(S(f))),
       typ="l",axes=FALSE,
       main="Figure 408")
  lines(xs,log10(ys[,3]),lty="dashed")
  lines(xs,log10(ys[,2]),lty="dotted")
  lines(xs,log10(ys[,1]),lty="longdash")
  abline(v=log10(1/6),lty="dotted")
  abline(v=log10(1/3),lty="dashed")
  axis(1, at=log10(c(1/100000,1/10000,1/1000,1/100,1/10,1)), labels=expression(10^-5,10^-4,10^-3,10^-2,10^-1,10^0))
  ttn <- 2:9
  axis(1, at=log10(c(ttn/100000,ttn/10000,ttn/1000,ttn/100,ttn/10)), labels=FALSE)
  axis(2, at=log10(c(1/10,1,10,100)), labels=expression(10^-1,10^0,10^1,10^2), las=2)
  axis(2, at=log10(c(ttn/10,ttn,ttn*10,ttn*100,ttn*1000)), labels=FALSE, tcl=-0.25)
  legend("bottomleft", legend=expression(delta==0.49,delta==0.45,delta==0.4,delta==0.25), lty=c("solid","dashed","dotted","longdash"),inset=0.02,cex=1/1.2)
  box(bty="l")
}

### Figure 408

fig_408()

### Figure 410 ###

fig_410 <- function(Ns=c(256,512,1024),
                    Ks=2:16,
                    the_labels=expression(italic(N==256),italic(N==512),italic(N==1024)),
                    pches=c(1,0,2),
                    bolds=c(16,15,17),
                    y_lim=c(0.8,1.4))
{
  N_curves <- length(Ns)
  ys <- sapply(Ks,function(a_K) var_power_law_alpha_mt(Ns[1],K=a_K)$var_alpha)/var_power_law_alpha_pgram(Ns[1])
  plot(Ks,ys,
       xlab=expression(italic(K)),
       ylim=y_lim,ylab=expression(paste("ratio of variances of ",alpha," estimators")),
       typ="b",pch=pches[1],cex=0.75,axes=FALSE,
       main="Figure 410")
  i_min <- which.min(ys)
  points(Ks[i_min],ys[i_min],pch=bolds[1],cex=0.75)
  text(13,ys[12],the_labels[1],pos=1)
  abline(h=1,lty="dashed")
  for(n in 2:N_curves) 
  {
    ys <- sapply(Ks,function(a_K) var_power_law_alpha_mt(Ns[n],K=a_K)$var_alpha)/var_power_law_alpha_pgram(Ns[n])
    points(Ks,ys,pch=pches[n],typ="b",cex=0.75)
    i_min <- which.min(ys)
    points(Ks[i_min],ys[i_min],pch=bolds[n],cex=0.75)
    text(13,ys[12],the_labels[n],pos=1)
  }
  axis(1,at=seq(min(Ks),max(Ks),2))
  axis(1,at=seq(min(Ks),max(Ks),1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(y_lim[1],y_lim[2],0.1),las=2)
  box(bty="l")
}

### Figure 410

fig_410()

### Tables 411a and 411b ###

### CAUTION: evalation of following commented-out 20 forms recreates
###          the contents of Tables 411a and 411b, but is time consuming;
###          use of the 20 uncommented `load' forms that follow results in
###          OLS_0p25_hf2_256 etc. being set quickly.
###
### OLS_0p25_hf2_256    <- do_OLS_pgram_mt(0.25,high_freq=1/3,N_ts=256)
### OLS_0p25_hf2_512    <- do_OLS_pgram_mt(0.25,high_freq=1/3,N_ts=512)
### OLS_0p25_hf2        <- do_OLS_pgram_mt(0.25,high_freq=1/3)
### OLS_0p25_hf2_512_K7 <- do_OLS_pgram_mt(0.25,high_freq=1/3,N_ts=512,K=7)
### OLS_0p25_hf2_K9     <- do_OLS_pgram_mt(0.25,high_freq=1/3,K=9)
###                                         .
### OLS_0p4_hf2_256     <- do_OLS_pgram_mt(0.4,high_freq=1/3,N_ts=256)
### OLS_0p4_hf2_512     <- do_OLS_pgram_mt(0.4,high_freq=1/3,N_ts=512)
### OLS_0p4_hf2         <- do_OLS_pgram_mt(0.4,high_freq=1/3)
### OLS_0p4_hf2_512_K7  <- do_OLS_pgram_mt(0.4,high_freq=1/3,N_ts=512,K=7)
### OLS_0p4_hf2_K9      <- do_OLS_pgram_mt(0.4,high_freq=1/3,K=9)
###                                         .
### OLS_0p45_hf2_256    <- do_OLS_pgram_mt(0.45,high_freq=1/3,N_ts=256)
### OLS_0p45_hf2_512    <- do_OLS_pgram_mt(0.45,high_freq=1/3,N_ts=512)
### OLS_0p45_hf2        <- do_OLS_pgram_mt(0.45,high_freq=1/3)
### OLS_0p45_hf2_512_K7 <- do_OLS_pgram_mt(0.45,high_freq=1/3,N_ts=512,K=7)
### OLS_0p45_hf2_K9     <- do_OLS_pgram_mt(0.45,high_freq=1/3,K=9)
###                                         .
### OLS_0p49_hf2_256    <- do_OLS_pgram_mt(0.49,high_freq=1/3,N_ts=256)
### OLS_0p49_hf2_512    <- do_OLS_pgram_mt(0.49,high_freq=1/3,N_ts=512)
### OLS_0p49_hf2        <- do_OLS_pgram_mt(0.49,high_freq=1/3)
### OLS_0p49_hf2_512_K7 <- do_OLS_pgram_mt(0.49,high_freq=1/3,N_ts=512,K=7)
### OLS_0p49_hf2_K9     <- do_OLS_pgram_mt(0.49,high_freq=1/3,K=9)

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p25_hf2_256.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p25_hf2_512.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p25_hf2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p25_hf2_512_K7.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p25_hf2_K9.Rdata"))

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p4_hf2_256.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p4_hf2_512.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p4_hf2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p4_hf2_512_K7.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p4_hf2_K9.Rdata"))

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p45_hf2_256.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p45_hf2_512.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p45_hf2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p45_hf2_512_K7.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p45_hf2_K9.Rdata"))

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p49_hf2_256.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p49_hf2_512.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p49_hf2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p49_hf2_512_K7.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/OLS_0p49_hf2_K9.Rdata"))

### Table 411a, row 1

round(c(var_biased(d_to_alpha(OLS_0p25_hf2_256$d_ests_pgram))/var_power_law_alpha_pgram(256),
        var_biased(d_to_alpha(OLS_0p4_hf2_256$d_ests_pgram))/var_power_law_alpha_pgram(256),
        var_biased(d_to_alpha(OLS_0p45_hf2_256$d_ests_pgram))/var_power_law_alpha_pgram(256),
        var_biased(d_to_alpha(OLS_0p49_hf2_256$d_ests_pgram))/var_power_law_alpha_pgram(256)),
      2)

### Table 411a, row 2

round(c(var_biased(d_to_alpha(OLS_0p25_hf2_512$d_ests_pgram))/var_power_law_alpha_pgram(512),
        var_biased(d_to_alpha(OLS_0p4_hf2_512$d_ests_pgram))/var_power_law_alpha_pgram(512),
        var_biased(d_to_alpha(OLS_0p45_hf2_512$d_ests_pgram))/var_power_law_alpha_pgram(512),
        var_biased(d_to_alpha(OLS_0p49_hf2_512$d_ests_pgram))/var_power_law_alpha_pgram(512)),
      2)

### Table 411a, row 3

round(c(var_biased(d_to_alpha(OLS_0p25_hf2$d_ests_pgram))/var_power_law_alpha_pgram(1024),
        var_biased(d_to_alpha(OLS_0p4_hf2$d_ests_pgram))/var_power_law_alpha_pgram(1024),
        var_biased(d_to_alpha(OLS_0p45_hf2$d_ests_pgram))/var_power_law_alpha_pgram(1024),
        var_biased(d_to_alpha(OLS_0p49_hf2$d_ests_pgram))/var_power_law_alpha_pgram(1024)),
      2)

### Table 411a, row 4

round(c(var_biased(d_to_alpha(OLS_0p25_hf2_256$d_ests_mt[,3]))/var_power_law_alpha_mt(256)$var_alpha,
        var_biased(d_to_alpha(OLS_0p4_hf2_256$d_ests_mt[,3]))/var_power_law_alpha_mt(256)$var_alpha,
        var_biased(d_to_alpha(OLS_0p45_hf2_256$d_ests_mt[,3]))/var_power_law_alpha_mt(256)$var_alpha,
        var_biased(d_to_alpha(OLS_0p49_hf2_256$d_ests_mt[,3]))/var_power_law_alpha_mt(256)$var_alpha),
      2)

### Table 411a, row 5

round(c(var_biased(d_to_alpha(OLS_0p25_hf2_512_K7$d_ests_mt[,4]))/var_power_law_alpha_mt(512)$var_alpha,
        var_biased(d_to_alpha(OLS_0p4_hf2_512_K7$d_ests_mt[,4]))/var_power_law_alpha_mt(512)$var_alpha,
        var_biased(d_to_alpha(OLS_0p45_hf2_512_K7$d_ests_mt[,4]))/var_power_law_alpha_mt(512)$var_alpha,
        var_biased(d_to_alpha(OLS_0p49_hf2_512_K7$d_ests_mt[,4]))/var_power_law_alpha_mt(512)$var_alpha),
      2)

### Table 411a, row 6

round(c(var_biased(d_to_alpha(OLS_0p25_hf2_K9$d_ests_mt[,5]))/var_power_law_alpha_mt(1024)$var_alpha,
        var_biased(d_to_alpha(OLS_0p4_hf2_K9$d_ests_mt[,5]))/var_power_law_alpha_mt(1024)$var_alpha,
        var_biased(d_to_alpha(OLS_0p45_hf2_K9$d_ests_mt[,5]))/var_power_law_alpha_mt(1024)$var_alpha,
        var_biased(d_to_alpha(OLS_0p49_hf2_K9$d_ests_mt[,5]))/var_power_law_alpha_mt(1024)$var_alpha),
      2)

### Table 411b, row 1

round(c(var(d_to_alpha(OLS_0p25_hf2_256$d_ests_mt[,3]))/var(d_to_alpha(OLS_0p25_hf2_256$d_ests_pgram)),
        var(d_to_alpha(OLS_0p4_hf2_256$d_ests_mt[,3]))/var(d_to_alpha(OLS_0p4_hf2_256$d_ests_pgram)),
        var(d_to_alpha(OLS_0p45_hf2_256$d_ests_mt[,3]))/var(d_to_alpha(OLS_0p45_hf2_256$d_ests_pgram)),
        var(d_to_alpha(OLS_0p49_hf2_256$d_ests_mt[,3]))/var(d_to_alpha(OLS_0p49_hf2_256$d_ests_pgram))),
      2)

### Table 411b, row 2

round(c(var(d_to_alpha(OLS_0p25_hf2_512_K7$d_ests_mt[,4]))/var(d_to_alpha(OLS_0p25_hf2_512$d_ests_pgram)),
        var(d_to_alpha(OLS_0p4_hf2_512_K7$d_ests_mt[,4]))/var(d_to_alpha(OLS_0p4_hf2_512$d_ests_pgram)),
        var(d_to_alpha(OLS_0p45_hf2_512_K7$d_ests_mt[,4]))/var(d_to_alpha(OLS_0p45_hf2_512$d_ests_pgram)),
        var(d_to_alpha(OLS_0p49_hf2_512_K7$d_ests_mt[,4]))/var(d_to_alpha(OLS_0p49_hf2_512$d_ests_pgram))),
      2)

### Table 411b, row 3

round(c(var(d_to_alpha(OLS_0p25_hf2_K9$d_ests_mt[,5]))/var(d_to_alpha(OLS_0p25_hf2$d_ests_pgram)),
        var(d_to_alpha(OLS_0p4_hf2_K9$d_ests_mt[,5]))/var(d_to_alpha(OLS_0p4_hf2$d_ests_pgram)),
        var(d_to_alpha(OLS_0p45_hf2_K9$d_ests_mt[,5]))/var(d_to_alpha(OLS_0p45_hf2$d_ests_pgram)),
        var(d_to_alpha(OLS_0p49_hf2_K9$d_ests_mt[,5]))/var(d_to_alpha(OLS_0p49_hf2$d_ests_pgram))),
      2)

### Figure 413 ###

fig_413 <- function(k=NULL,ts=ts_100)
{
  if(is.null(k))
  {
    xs <- 0:(length(ts_100)-1)
    ys <- ts
    y_lab <- "time series"
  }
  else
  {
    xs <- (k*17):((k*17)+31)
    ys <- ts[xs+1]
    y_lab <- paste("block",k)
  }
  plot(xs,ys,
       xlim=c(0,length(ts)),xlab=expression(italic(t)),
       ylim=c(-1,1),ylab=y_lab,
       typ="l",col="gray50",axes=FALSE,
       main="Figure 413")
  points(xs,ys,pch=16,cex=0.5)
  if(is.null(k)) for(n in 0:4) segments(c(n*17,n*17+31),c(-1,-1)+n*0.17,c(n*17,n*17+31),c(1,1)-n*0.17, col=c("black","gray20","gray40","gray60","gray80")[n+1],lwd=1.25)
  axis(1,at=seq(0,length(ts),20))
  axis(2,at=seq(-1,1,1),las=2)
  box(bty="l")
}

### Figure 413, plots from top to bottom

fig_413()
fig_413(0)
fig_413(1)
fig_413(2)
fig_413(3)
fig_413(4)

### Figure 415 ###

fig_415 <- function(N_S=c(64,128,256),N=1024,
                    tags=expression(italic(N[S] == 256),italic(N[S] == 128),italic(N[S] == 64)),
                    y_text=c(6.5,23.5,56),
                    y_max=64)
{
  compute_edofs <- function(N_S, N)
  {
    h_p <- c(hanning_taper(N_S),rep(0,N-N_S))  # taper padded with zeros ...
    possible_ns <- 1:N_S
    temp <- (N-N_S)/possible_ns + 1
    good_N_Bs <- temp%%1 == 0
    N_Bs <- temp[good_N_Bs]
    ns <- (N-N_S)/(N_Bs-1)
    N_ns <- length(ns)
    edofs <- rep(0,N_ns)
    overlap <- rep(0,N_ns)
    my_conv <- function(m) sum(h_p[1:N_S]*h_p[(1+m*n):(m*n+N_S)])^2
    for(k in 1:N_ns)
    {
      n <- ns[k]
      N_B <- N_Bs[k]
      ms <- 1:(N_B-1)
      edofs[k] <- 2*N_B/(1 + 2*sum((1-ms/N_B)*sapply(ms,my_conv)))
      overlap[k] <- (1-n/N_S)*100
    }
    return(list(x=rev(overlap),y=rev(edofs)))
  }
  plot(1:2,1:2,
       xlim=c(0,100),xlab="overlap percentage",
       ylim=c(0,y_max),ylab="EDOF",
       typ="n",col="gray50",axes=FALSE,
       main="Figure 415")
  for(n in 1:length(N_S))
  {
    xs_and_ys <- compute_edofs(N_S[n], N)
    lines(xs_and_ys$x,xs_and_ys$y,col="gray50")
    points(xs_and_ys$x,xs_and_ys$y,cex=0.5)
    i <- which.max(xs_and_ys$y) 
    points(xs_and_ys$x[i],xs_and_ys$y[i],pch=16,col="red",cex=0.5)
    text(103.5,y_text[n],tags[n],pos=2)
  }
  axis(1,at=c(0,50,100),labels=c("0%","50%","100%"))
  axis(1, at=seq(0,100,10), labels=FALSE, tcl=-0.25)
  axis(2,at=seq(0,y_max,10),las=2)
  box(bty="l")
}

### Figure 415

fig_415()

### Figure 417 ###

fig_417 <- function(N_s,tag,ts=ar2_1,coeffs=ar2_coeffs,innov_var=ar2_innov_var)
{
  true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
  plot(true_sdf$freqs,dB(true_sdf$sdf),
       xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-40,20),yaxs="i",ylab="AR(2) spectra  (dB)",
       typ="l",lwd=0.5,axes=FALSE,
       main="Figure 417")
  est_sdf <- wosa_sdf_est(ts,hanning_taper(N_s),center=FALSE,pad_factor=1024/N_s)
  lines(est_sdf$freqs, dB(est_sdf$sdfe))
  cc <- est_sdf$cc
  x_left <- 0.025
  x_mid <- x_left + cc$width/2
  y_cc <- -30.0
  lines(c(x_left,x_left+cc$width), c(y_cc,y_cc), lwd=0.5)
  lines(c(x_mid,x_mid), c(y_cc-cc$down,y_cc+cc$up), lwd=0.5)
  axis(1,at=seq(0,0.5,0.5))
  axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-40,20,20),las=2)
  axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
  text(0.49,15,tag,pos=2)
  box(bty="l")
}

### Figure 417, top row

fig_417(4,expression(italic(N[S] == 4)))
fig_417(16,expression(italic(N[S] == 16)))

### Figure 417, bottom row

fig_417(32,expression(italic(N[S] == 32)))
fig_417(64,expression(italic(N[S] == 64)))

### Figure 418 ###

fig_418 <- function(N_s,tag,ts=ar4_1,coeffs=ar4_coeffs,innov_var=ar4_innov_var)
{
  true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
  plot(true_sdf$freqs,dB(true_sdf$sdf),
       xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-60,20),yaxs="i",ylab="AR(4) spectra  (dB)",
       typ="l",lwd=0.5,axes=FALSE,
       main="Figure 418")
  est_sdf <- wosa_sdf_est(ts,hanning_taper(N_s),center=FALSE,pad_factor=1024/N_s)
  lines(est_sdf$freqs, dB(est_sdf$sdfe))
  cc <- est_sdf$cc
  x_left <- 0.025
  x_mid <- x_left + cc$width/2
  y_cc <- -30.0
  lines(c(x_left,x_left+cc$width), c(y_cc,y_cc), lwd=0.5)
  lines(c(x_mid,x_mid), c(y_cc-cc$down,y_cc+cc$up), lwd=0.5)
  axis(1,at=seq(0,0.5,0.5))
  axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-60,20,20),las=2)
  axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
  text(0.49,15,tag,pos=2)
  box(bty="l")
}

### Figure 418, top row

fig_418(64,expression(italic(N[S] == 64)))
fig_418(128,expression(italic(N[S] == 128)))

### Figure 418, bottom row

fig_418(256,expression(italic(N[S] == 256)))
fig_418(512,expression(italic(N[S] == 512)))

### Figure 419 ###

fig_419 <- function(the_tapers,tag_1,tag_2,ts=ar4_1,coeffs=ar4_coeffs,innov_var=ar4_innov_var)
{
  true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
  plot(true_sdf$freqs,dB(true_sdf$sdf),
       xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-60,20),yaxs="i",ylab="AR(4) spectra  (dB)",
       typ="l",lwd=0.5,axes=FALSE,
       main="Figure 419")
  mtse <- mt_sdf_est(ts,the_tapers,center=FALSE,pad_factor=1)
  lines(mtse$freqs, dB(mtse$sdfe))
  cc <- mtse$cc
  x_left <- 0.025
  x_mid <- x_left + cc$width/2
  y_cc <- -30.0
  lines(c(x_left,x_left+cc$width), c(y_cc,y_cc), lwd=0.5)
  lines(c(x_mid,x_mid), c(y_cc-cc$down,y_cc+cc$up), lwd=0.5)
  axis(1,at=seq(0,0.5,0.5))
  axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-60,20,20),las=2)
  axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
  text(0.49,15,tag_1,pos=2)
  text(0.01,-53.333,tag_2,pos=4)
  box(bty="l")
}

### Figure 419, top row

fig_419(t(as.matrix(taper("dpss",1024,n.taper=3,bandwidth=4))),expression(italic(K == 3)),expression(paste(italic(NW==4)," Slepian")))
fig_419(t(as.matrix(taper("sine",1024,n.taper=3))),expression(italic(K == 3)),"sinusoidal")

### Figure 419, bottom row

fig_419(t(as.matrix(taper("dpss",1024,n.taper=7,bandwidth=5.5))),expression(italic(K == 7)),expression(paste(italic(NW==5.5)," Slepian")))
fig_419(t(as.matrix(taper("sine",1024,n.taper=7))),expression(italic(K == 7)),"sinusoidal")

### Figure 420 ###

fig_420 <- function(N=1024,m_parzen=179,N_S=256,N_eigenvalues=16)
{
  the_taper <- hanning_taper(N)
  the_lw_2_sided <- sapply((-(N-1)):(N-1),function(tau) parzen_lag_window(tau,m_parzen))
  Q_lw <- the_taper %o% the_taper
  for(i in 1:N) Q_lw[i,] <- Q_lw[i,]*the_lw_2_sided[((N+1)-i):(2*N-i)]
  eigen_lw <- eigen(Q_lw,TRUE)$values[1:N_eigenvalues]
  N_dse <- 2*N/N_S - 1  # number of direct SDF estimates used to form WOSA (Hanning, 50% overlap)
  the_taper_wosa <- c(hanning_taper(N_S),rep(0,N-N_S))/sqrt(N_dse)
  Q_wosa <- the_taper_wosa %o% the_taper_wosa
  for(i in 1:(N_dse-1))
  {
    taper_shifted <- circular_shift(the_taper_wosa,i*N_S/2)
    Q_wosa <- Q_wosa + taper_shifted %o% taper_shifted
  }
  eigen_wosa <- eigen(Q_wosa,TRUE)$values[1:N_eigenvalues]
  xs <- 0:(N_eigenvalues-1)
  y_lim <- range(c(eigen_lw,eigen_wosa))
  plot(xs,eigen_lw,
       xlab=expression(italic(k)),
       ylim=y_lim,ylab="eigenvalue",
       typ="p",axes=FALSE,
       main="Figure 420")
  points(xs,eigen_wosa,pch="+")
  axis(1,at=seq(0,N_eigenvalues-1,5))
  axis(1,at=seq(0,N_eigenvalues-1,1),labels=FALSE,tcl=-0.25)
  axis(2,at=seq(0,round(y_lim[2],1),0.05),las=2)
  box(bty="l")
}

### Figure 420

fig_420()

### Figure 421 ###

fig_421 <- function(k,tag=" ",N=1024,m_parzen=179,N_S=256,y_lim=c(-0.066,0.075),y_at=0.05,y_text=0.04,wosa_p=FALSE)
{
  if(wosa_p)
  {
    N_dse <- 2*N/N_S - 1  # number of direct SDF estimates used to form WOSA (Hanning, 50% overlap)
    the_taper_wosa <- c(hanning_taper(N_S),rep(0,N-N_S))/sqrt(N_dse)
    Q_wosa <- the_taper_wosa %o% the_taper_wosa
    for(i in 1:(N_dse-1))
    {
      taper_shifted <- circular_shift(the_taper_wosa,i*N_S/2)
      Q_wosa <- Q_wosa + taper_shifted %o% taper_shifted
    }
    eigen_taper <- eigen(Q_wosa,TRUE)$vec[,k+1]
  }
  else
  {
    the_taper <- hanning_taper(N)
    the_lw_2_sided <- sapply((-(N-1)):(N-1),function(tau) parzen_lag_window(tau,m_parzen))
    Q_lw <- the_taper %o% the_taper
    for(i in 1:N) Q_lw[i,] <- Q_lw[i,]*the_lw_2_sided[((N+1)-i):(2*N-i)]
    eigen_taper <- eigen(Q_lw)$vec[,k+1]
  }
  plot(0:(N-1),eigen_taper*(if(eigen_taper[1]>0) 1 else -1),
       xlim=c(0,N),xaxs="i",xlab=expression(italic(t)),
       ylim=y_lim,ylab="eigenvector",
       typ="l",axes=FALSE,
       main="Figure 421")
  axis(1,at=seq(0,N,N/4))
  axis(2,at=y_at*c(-1,1),las=2)
  axis(2,at=0,las=2)
  text(N,y_text,tag,pos=2)
  box(bty="l")
}

### Figure 421, left-hand column

fig_421(0,tag=expression(italic(k==0)))
fig_421(1,tag=expression(italic(k==1)))
fig_421(2,tag=expression(italic(k==2)))
fig_421(3,tag=expression(italic(k==3)))
fig_421(4,tag=expression(italic(k==4)),y_text=-0.03)
fig_421(5,tag=expression(italic(k==5)))
fig_421(6,tag=expression(italic(k==6)),y_text=-0.03)

### Figure 421, right-hand column

fig_421(0,wosa_p=TRUE)
fig_421(1,wosa_p=TRUE)
fig_421(2,wosa_p=TRUE)
fig_421(3,wosa_p=TRUE)
fig_421(4,wosa_p=TRUE)
fig_421(5,wosa_p=TRUE)
fig_421(6,wosa_p=TRUE)

### Figure 422 ###

fig_422 <- function(L,tag,ts=ar4_1,m_parzen=179,N_S=256,wosa_p=FALSE)
{
  N <- length(ts)
  if(wosa_p)
  {
    N_dse <- 2*N/N_S - 1  # number of direct SDF estimates used to form WOSA (Hanning, 50% overlap)
    the_taper_wosa <- c(hanning_taper(N_S),rep(0,N-N_S))/sqrt(N_dse)
    Q_wosa <- the_taper_wosa %o% the_taper_wosa
    for(i in 1:(N_dse-1))
    {
      taper_shifted <- circular_shift(the_taper_wosa,i*N_S/2)
      Q_wosa <- Q_wosa + taper_shifted %o% taper_shifted
    }
    eigen_stuff <- eigen(Q_wosa)
    lw_or_wosa <- wosa_sdf_est(ts,hanning_taper(N_S),center=FALSE,pad_factor=N/N_S)
  }
  else
  {
    the_taper <- hanning_taper(N)
    the_lw_2_sided <- sapply((-(N-1)):(N-1),function(tau) parzen_lag_window(tau,m_parzen))
    Q_lw <- the_taper %o% the_taper
    for(i in 1:N) Q_lw[i,] <- Q_lw[i,]*the_lw_2_sided[((N+1)-i):(2*N-i)]
    eigen_stuff <- eigen(Q_lw)
    lw_or_wosa <- ts_to_lag_window_sdf_est(ts,m_parzen,lag_window=parzen_lag_window,taper=hanning_taper,center=FALSE) 
  }
  wmtse <- mt_sdf_est(ts,eigen_stuff$vec[,1:L],center=FALSE,weights=eigen_stuff$values[1:L]/sum(eigen_stuff$values[1:L]))
  plot(wmtse$freqs,dB(wmtse$sdf),
       xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
       ylim=c(-60,20),yaxs="i",ylab=if(wosa_p) "WOSA spectra  (dB)" else "lag window spectra  (dB)", 
       typ="l",lwd=0.5,axes=FALSE,
       main="Figure 422")
  lines(lw_or_wosa$freqs,dB(lw_or_wosa$sdf))
  axis(1,at=seq(0,0.5,0.5))
  axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-60,20,20),las=2)
  axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
  text(0.49,15,tag,pos=2)
  box(bty="l")
}

### Figure 422, left-hand column

fig_422(3,expression(italic(L==3)))
fig_422(3,expression(italic(L==3)),wosa_p=TRUE)

### Figure 422, right-hand column

fig_422(7,expression(italic(L==7)))
fig_422(6,expression(italic(L==6)),wosa_p=TRUE)

### Figure 426 ###

fig_426 <- function(ts,the_tapers,tag_1,tag_2,m=150,lw=parzen_lag_window,taper=slepian_taper(length(ts),2),delta_t=1/4)
{
  mtse <- mt_sdf_est(ts,the_tapers,center=TRUE,delta_t=delta_t)
  plot(mtse$freqs,dB(mtse$sdfe),
       xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
       ylim=c(-40,80),yaxs="i",ylab="dB",
       typ="l",axes=FALSE,
       main=paste("Figure 426",tag_2,sep=""))
  cc <- mtse$cc
  x_cc <- 0.25
  y_cc <- 25
  lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
  lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
  if(tag_2 == "(d)")
  {
    lwe <- ts_to_lag_window_sdf_est(ts,m,taper=taper,lag_window=lw,center=TRUE,delta_t=delta_t)
    lines(lwe$freqs, dB(lwe$sdfe), lwd=0.5, col="gray40")
    cc <- lwe$cc
    y_cc <- 15
    lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5, col="gray40")
    lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5, col="gray40")
    legend("bottomleft", legend=expression(paste("sinusoidal, ",italic(K == 12)),paste("Parzen, ",italic(m == 150))),col=c("black","gray40"),lwd=c(1.0,0.5),inset=0.02,cex=1.2,box.col="gray80")
  }
  axis(1,at=seq(0,2.0,0.5))
  axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-40,80,20),las=2)
  axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
  text(1.0,80,tag_1,pos=1)
  text(1.9,80,tag_2,pos=1)
  box(bty="l")
}

### Figure 426, upper plot

fig_426(ocean_wave,t(as.matrix(taper("dpss",length(ocean_wave),n.taper=6,bandwidth=4))),expression(paste("Slepian multitaper, ", italic(NW==4/Delta[t]),", ",italic(K==6))),"(a)")

### Figure 426, 2nd plot

fig_426(ocean_wave,t(as.matrix(taper("dpss",length(ocean_wave),n.taper=12,bandwidth=7))),expression(paste("Slepian multitaper, ", italic(NW==7/Delta[t]),", ",italic(K==12))),"(b)")

### Figure 426, 3rd plot

fig_426(ocean_wave,t(as.matrix(taper("sine",length(ocean_wave),n.taper=12))),expression(paste("sinusoidal multitaper, ",italic(K==12))),"(c)")

### Figure 426, 4th plot

fig_426(ocean_wave,t(as.matrix(taper("sine",length(ocean_wave),n.taper=12)))," ","(d)")

### Figure 428 ###

fig_428 <- function(ts,the_tapers,NW,tag_1,tag_2=" ",delta_t=1/4)
{
  amtse <- amt_sdf_est(ts,the_tapers,NW,center=TRUE,delta_t=delta_t)
  if(tag_1 == "(a)")
  {        
    plot(amtse$freqs,dB(amtse$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 428(a)")
  }
  else
    if(tag_1 == "(b)")
    {        
      plot(amtse$freqs,dB(amtse$ci_upper),
           xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
           ylim=c(-40,80),yaxs="i",ylab="dB",
           typ="l",axes=FALSE,
           main="Figure 428(b)")
      lines(amtse$freqs,dB(amtse$ci_lower))
    }
  if(tag_1 == "(c)")
  {        
    plot(amtse$freqs,amtse$edofs,
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(0,16),yaxs="i",ylab="EDOFs",
         typ="l",axes=FALSE,
         main="Figure 428(c)")
    axis(2,at=seq(0,16,4),las=2)
    axis(2,at=seq(0,16,2),label=FALSE,tcl=-0.25)
  }
  else
  {
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
  }
  axis(1,at=seq(0,2.0,0.5))
  axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
  text(1.9,80,tag_1,pos=1)
  text(1.0,80,tag_2,pos=1)
  box(bty="l")
}

### Figure 428, upper plot

fig_428(ocean_wave,t(as.matrix(taper("dpss",length(ocean_wave),n.taper=7,bandwidth=4))),4,"(a)",expression(paste("adaptive multitaper, ", italic(NW==4/Delta[t]),", ",italic(K==7))))

### Figure 428, middle plot

fig_428(ocean_wave,t(as.matrix(taper("dpss",length(ocean_wave),n.taper=7,bandwidth=4))),4,"(b)")

### Figure 428, lower plot

fig_428(ocean_wave,t(as.matrix(taper("dpss",length(ocean_wave),n.taper=7,bandwidth=4))),4,"(c)")

### Figure 429 ###

fig_429 <- function(delta_t=0.25)
{
  B_U_gaussian <- B_U(gaussian_lag_window(0:1023,23.67),slepian_taper(1024,2),delta_t=delta_t)
  BSs <- 32:128
  BHs <- sapply(BSs,function(a_BS) B_H(hanning_taper(a_BS),delta_t=delta_t))
  plot(BSs,BHs,
       xlab=expression(paste("block size (",italic(N[S]),")")),
       ylab="bandwidth",
       typ="p",axes=FALSE,
       main="Figure 429")
  abline(h=B_U_gaussian,lty="dashed")
  temp <- abs(BHs - B_U_gaussian)
  i_closest <- which(min(temp) == temp)
  points(BSs[i_closest], BHs[i_closest], pch=16)
  axis(1,at=seq(40,120,20))
  axis(1,at=seq(30,130,10),label=FALSE,tcl=-0.25)
  axis(2,at=seq(0.1,0.25,0.05),las=2)
  box(bty="l")
}

### Figure 429

fig_429()

### Figure 430 ###

fig_430 <- function(ts,m=23.666,lw=gaussian_lag_window,lw_taper=slepian_taper(length(ts),2),wosa_taper=hanning_taper(61),delta_t=1/4)
{
  N <- length(ts)
  N_S <- length(wosa_taper)
  stops <- seq(N_S,N,floor(N_S/2))
  stops[length(stops)] <- N
  starts <- stops - N_S + 1
  wosae <- wosa_sdf_est(ts,wosa_taper,starts,center=TRUE,pad_factor=1024/N_S,delta_t=delta_t)
  plot(wosae$freqs,dB(wosae$sdfe),
       xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
       ylim=c(-40,80),yaxs="i",ylab="dB",
       typ="l",axes=FALSE,
       main="Figure 430")
  x_cc <- 0.25
  y_cc <- 25
  cc <- wosae$cc
  lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
  lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
  ##
  lwe <- ts_to_lag_window_sdf_est(ts,m,taper=lw_taper,lag_window=lw,center=TRUE,delta_t=delta_t)
  lines(lwe$freqs,dB(lwe$sdfe),lwd=0.5,col="gray40")
  y_cc <- 15
  cc <- lwe$cc
  lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5, col="gray40")
  lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5, col="gray40")
  ##
  axis(1,at=seq(0,2.0,0.5))
  axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
  axis(2,at=seq(-40,80,20),las=2)
  axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
  box(bty="l")
  legend("topright", legend=c("WOSA","Gaussian lag window"),col=c("black","gray40"),lwd=c(1.0,0.5),inset=0.02,box.col="gray80")
}

### Figure 430 ###

fig_430(ocean_wave)

### Figure 431 ###

fig_431 <- function(ff,tag="(a)",K=7,M=140,x_upper_lim_log=log10(0.06))
{
  mtse <- mt_sdf_est(ff,t(as.matrix(taper("sine",length(ff),n.taper=K))),center=TRUE)
  if(tag == "(a)")
  {
    plot(mtse$freqs,dB(mtse$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/minute)")),
         ylim=c(-30,-5),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 431(a)")
    x_cc <- 0.3
    y_cc <- -25.0
    lines(c(x_cc,x_cc),y_cc+c(mtse$cc$up,-mtse$cc$down),lwd=0.5)
    lines(x_cc+c(-mtse$cc$width/2,mtse$cc$width/2),c(y_cc,y_cc),lwd=0.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    x_text <- 0.001
  }
  else
  {
    log10_freqs <- log10(mtse$freqs[-1])         
    dB_mtse <- dB(mtse$sdfe[-1])        
    M_0 <- ceiling(K/2)
    plot(log10_freqs[1:(M_0-1)],dB_mtse[1:(M_0-1)],
         xlim=c(log10_freqs[1],x_upper_lim_log),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/minute)")),
         ylim=c(-30,-5),yaxs="i",ylab="dB",
         typ="l",lwd=0.5,col="gray40",axes=FALSE,
         main="Figure 431(b)")
    lines(log10_freqs[M_0:M], dB_mtse[M_0:M],lwd=1)
    lines(log10_freqs[-(1:M)], dB_mtse[-(1:M)], lwd=0.5, col="gray40")
    reg_coeffs <- coef(lm(dB_mtse[M_0:M] ~ log10_freqs[M_0:M]))
    reg_xs <- log10(c(M_0,M)/length(ff))
    reg_ys <- reg_coeffs[1] + reg_coeffs[2]*reg_xs
    lines(reg_xs,reg_ys,lwd=1.75)
    ttn <- 2:9
    axis(1, at=log10(c(1/100000,1/10000,1/1000,1/100,1/10,1)), labels=expression(10^-5,10^-4,10^-3,10^-2,10^-1,10^0))
    axis(1, at=log10(c(ttn/100000,ttn/10000,ttn/1000,ttn/100,ttn/10)), labels=FALSE, tcl=-0.25)
    temp <- range(log10_freqs)
    x_text <- temp[1] + diff(temp) * 0.002
  }
  axis(2,at=seq(-30,-10,10),las=2)
  axis(2,at=seq(-30,-5,5),label=FALSE,tcl=-0.25)
  text(x_text,-6.5625,tag,pos=4)
  box(bty="l")
}

### Figure 431, left-hand plot

fig_431(ac_fractional_frequencies)

### Figure 431, right-hand plot

fig_431(ac_fractional_frequencies,"(b)")

### Table 432, top row

fd_ff <- diff(ac_fractional_frequencies)
round(do_it_innov_var(fd_ff-mean(fd_ff),Ls=2:25,Ks=2:25)$mt[c(4,9,14,19,24)],5)

### Table 432, bottom row

round(do_it_innov_var(ac_fractional_frequencies-mean(ac_fractional_frequencies),Ls=2:25,Ks=2:25)$mt[c(4,9,14,19,24)],5)

### Figure 443 ###

fig_443 <- function(ws)
{
  xs <- ws$V1
  ys <- ws$V2
  plot(xs,ys,
       xlim=c(0,100),xaxs="i",xlab=expression(paste(italic(lambda[k]),"  (km)")),
       ylim=c(0,0.2),yaxs="i",ylab="SDF",
       typ="p",pch=16,axes=FALSE,
       main="Figure 443")
  lines(xs,ys)
  axis(1,at=seq(0,100,20))
  axis(2,at=seq(0,0.2,0.05),las=2)
  box(bty="l")
}

### Figure 443

fig_443(wavelength_spec)


