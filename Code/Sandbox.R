##### Sandbox #####

#this is a script to try things out before they go in real scripts


#test out the functions in AVARSpectralEstimate file

X.t <- X.t_missing <- rnorm(2048)
X.t_missing[c(100:200, 700:950, 1500:1600)] <- NA

plot(X.t)

#get spectral estimate
spec.est <- spec.pgram(X.t, demean = TRUE)
test_AVAR <- AVAR_bpvar(spectral_est = spec.est$spec, taus = 2^(0:8))
test_AVAR

spec.est <- multitaper_est(X.t_missing,W=0.005, K = 5)


length(spec.est$spectrum)

t.vec <- 1:2048
test_var <- var_AVAR_trfunc(tapers = spec.est$tapers, t.n =  t.vec[which(is.na(X.t_missing))], X.t = na.omit(X.t_missing),taus = 2^(0:8))

rm(A.prime, A.prime.matlab, clock_MTSE, clock_MTSE_omitted, dist.mat, eig_vecs, eigdec, OSS, oss.mat, u.mat, uk.mat, uvec.mat, vk.mat, combined.ts, random.walk)