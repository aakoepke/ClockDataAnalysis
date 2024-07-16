###### Interpolation Error - Tara Fortier ####
library(openxlsx)
library(readr)

X20240426_Yb_Freq_Shifts <- read_delim("Code/Cait/20240426_Yb_Freq_Shifts.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE, skip = 7)


convertToDateTime(head(X20240426_Yb_Freq_Shifts$t_mjd), origin = "1858-11-17",tz="UTC")

X20240426_Yb_Freq_Shifts$t_mjd


N.long <- 300 + 50
X.t <- rnorm(N.long)
X.t[c(1:20, 90:115)] <- NA
plot(X.t, pch = 19)
