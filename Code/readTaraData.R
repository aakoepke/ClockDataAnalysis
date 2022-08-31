library(readr)
library(ggplot2)
library(dplyr)
###############################################################################################
### 04/03/2018 analyzed on 07/10/18 frac.diff. = (measured-ideal)/ideal * 1e15
### Columns: MJD  Al vs Yb frac.diff. from(2.1628871275166585754)  Sr vs Yb frac.diff. from(1.2075070393433381221)  Al vs Sr frac.diff. from(2.6117014317814574243)
###############################################################################################
dat180403 <- read_delim("Data/180403 optical analysis TiS_fixed.dat",
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE, skip = 2)
colnames(dat180403)=c("MJD","AlYb","SrYb","AlSr")


###############################################################################################
### 03/09/2018 analyzed on 07/10/18 frac.diff. = (measured-ideal)/ideal * 1e15
### Columns: MJD  Al vs Yb frac.diff. from(2.1628871275166585754)  Sr vs Yb frac.diff. from(1.2075070393433381221)  Al vs Sr frac.diff. from(2.6117014317814574243)
###############################################################################################
dat180309 <- read_delim("Data/180309 optical analysis TiS_fixed.dat",
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE, skip = 2)
colnames(dat180309)=c("MJD","AlYb","SrYb","AlSr")

###############################################################################################
### 03/05/2018 analyzed on 07/10/18 frac.diff. = (measured-ideal)/ideal * 1e15
### Columns: MJD  Al vs Yb frac.diff. from(2.1628871275166585754)  Sr vs Yb frac.diff. from(1.2075070393433381221)  Al vs Sr frac.diff. from(2.6117014317814574243)
###############################################################################################
dat180305 <- read_delim("Data/180305 optical analysis TiS_fixed.dat",
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE, skip = 2)
colnames(dat180305)=c("MJD","AlYb","SrYb","AlSr")

###############################################################################################
### 03/06/2018 analyzed on 07/10/18 frac.diff. = (measured-ideal)/ideal * 1e15
### Columns: MJD  Al vs Yb frac.diff. from(2.1628871275166585754)  Sr vs Yb frac.diff. from(1.2075070393433381221)  Al vs Sr frac.diff. from(2.6117014317814574243)
###############################################################################################
dat180306 <- read_delim("Data/180306 optical analysis TiS_fixed.dat",
                  delim = "\t", escape_double = FALSE,
                  col_names = FALSE, trim_ws = TRUE, skip = 2)
colnames(dat180306)=c("MJD","AlYb","SrYb","AlSr")

###############################################################################################
### combine all, MJD will keep separate
### I think we should start analyzing data just one day at a time, but might be interesting to look at multiday data at some point
###############################################################################################

alldat=bind_rows(dat180305,dat180306,dat180309,dat180403)

###############################################################################################
###############################################################################################
### Look at data
###############################################################################################
###############################################################################################

ggplot(alldat,aes(MJD,AlYb))+
  geom_point()

### plot all 3 together, not sure how useful this is
# library(tidyr)
# alldat_long <- gather(alldat, ratio, measurement, AlYb:AlSr, factor_key=TRUE)
# ggplot(alldat_long,aes(MJD,measurement,col=ratio))+
#   geom_point()
