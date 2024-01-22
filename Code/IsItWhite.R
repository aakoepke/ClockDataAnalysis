############ Is it White? ###########

library(readr)
#read in clock data
primerData1 <- read_delim("Data/primerData1.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE, skip = 9)
plot(primerData1$C)


dat180403 <- read_delim("//cfs2w.nist.gov/unix$/776unix/cmb15/ClockDataAnalysis/Data/180403 optical analysis TiS_fixed.dat",
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE, skip = 2)
colnames(dat180403)=c("MJD","AlYb","SrYb","AlSr")
dat
#look at histogram
## distribution normal looking?
hist(na.omit(dat180403$AlYb))
hist(na.omit(dat180309$SrYb))
hist(na.omit(dat180309$AlSr))

## Uncorrelated? 
acf(na.omit(dat180403$AlYb))
acf(na.omit(diff(dat180403$AlYb)))
acf(na.omit(dat180403$AlSr)) 
#there is clearly autocorrelation in these data which is the first sign to me that this is not white noise

## ACF of first differences just for fun, still not white looking
acf(diff(na.omit(dat180403$AlYb)))
acf(diff(na.omit(dat180403$SrYb)))
acf(diff(na.omit(dat180403$AlSr)))


#look at acf

wn1 <- rnorm(1000, mean = 0, sd = 1)
wn2 <- rnorm(1000, mean = 0, sd = 2)
plot(wn1)
plot(wn2)

plot(wn1/wn2)
hist(wn1/wn2, breaks = 200)
