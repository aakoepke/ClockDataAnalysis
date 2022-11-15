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

