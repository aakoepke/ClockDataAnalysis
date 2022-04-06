# some code here: http://faculty.washington.edu/dbp/sauts/R-code/chapter-08.R

library(readr)
library(ggplot2)

primerData1 <- read_csv("Data/primerData1.txt", col_names = FALSE, skip = 10)
colnames(primerData1)="X1tilde"
primerData1$t = 1:dim(primerData1)[1]

ggplot(primerData1,aes(t,X1tilde))+
  geom_point()

test=spectrum(primerData1$X1tilde)

res=data.frame(f=test$freq,spec=test$spec)

ggplot(res,aes(f,spec))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()


arSpec=spectrum(primerData1$X1tilde,method = "ar")
arSpec

arRes=data.frame(f=arSpec$freq,spec=arSpec$spec)

ggplot(arRes,aes(f,spec))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()



################ second dataset

primerData2 <- read_csv("Data/primerData2.txt", col_names = FALSE, skip = 10)
colnames(primerData2)="X2tilde"
primerData2$t = 1:dim(primerData2)[1]

ggplot(primerData2,aes(t,X2tilde))+
  geom_point()

test=spectrum(primerData2$X2tilde)

res=data.frame(f=test$freq,spec=test$spec)

ggplot(res,aes(f,spec))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()


arSpec=spectrum(primerData2$X2tilde,method = "ar")
arSpec

arRes=data.frame(f=arSpec$freq,spec=arSpec$spec)

ggplot(arRes,aes(f,spec))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()
