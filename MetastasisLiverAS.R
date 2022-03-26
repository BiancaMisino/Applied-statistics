

library(readxl)
library(carData)


dim(Database_TRG_anonymized_06_03_22_Cliniche)
n <- dim(Database_TRG_anonymized_06_03_22_Cliniche)[1]
p <- dim(Database_TRG_anonymized_06_03_22_Cliniche)[2]


DB_cliniche2<-na.omit(Database_TRG_anonymized_06_03_22_Cliniche)

DB_cliniche2_label<-DB_cliniche2[,1:4]
DB_cliniche2<-DB_cliniche2[,-(1:4)]
DB_cliniche2_num<-DB_cliniche2[,1:24]
DB_cliniche2_num<-DB_cliniche2[,-c(2,4,5,6,8,11:18,19,21,22,23,24)]


  
x11()
pairs(DB_cliniche2_num)
  
#exploration

x11()
boxplot(DB_cliniche2_num,col='gold')

x11()

boxplot(scale(x=DB_cliniche2_num,center = F, scale=F), las=2, col='gold',outline=FALSE)

##c'è outlier in CEA, leviamolo



NDB <- DB_cliniche2_num[-43,]


x11()

boxplot(scale(x=NDB,center = F, scale=F), las=2, col='gold')




#Q <- quantile(DB_cliniche2_num$CEA, probs=c(.25, .75), na.rm = FALSE)
#iqr <- IQR(DB_cliniche2_num$CEA)
#up <-  Q[2]+1.5*iqr # Upper Range  
#low<- Q[1]-1.5*iqr # Lower Range

#eliminated<- subset(DB_cliniche2_num, DB_cliniche2_num$CEA > (Q[1] - 1.5*iqr) & DB_cliniche2_num$CEA < (Q[2]+1.5*iqr)
#DB_cliniche2_num(eliminated, Interval CT-Surgery,days breaks, outlier.tagging = TRUE) 


# We perform the PCA on std data

picia <- princomp(NDB,cor=TRUE,scores = TRUE)
summary(picia)


# loadings (recall: coefficients of the linear combination of the original 
#           variables that defines each principal component)
load.tour <- picia$loadings
load.tour

load.tour[,1:6]


# graphical representation of the loadings of the first six principal components
x11()
par(mfcol = c(2,3))
for(i in 1:6) barplot(load.tour[,i], ylim = c(-1, 1), main=paste("PC",i))## bars proportional to loadings 



##Explained variance

x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(picia, las=2, main='Principal Components', ylim=c(0,7))
abline(h=1, col='blue')
barplot(sapply(NDB.sd,sd)^2, las=2, main='Original Variables', ylim=c(0,7), ylab='Variances')
plot(cumsum(picia$sde^2)/sum(picia$sde^2), type='b', axes=F, xlab='Number of components', ylab='Contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(NDB.sd),labels=1:ncol(NDB.sd),las=2)


scores.NDB <- picia$scores
scores.NDB

x11()
plot(scores.NDB[,1:4])
abline(h=0, v=0, lty=2, col='grey')

x11()
layout(matrix(c(1,2),2))
boxplot(NDB.sd, las=2, col='gold', main='Standardized variables')
scores.NDB <- data.frame(scores.NDB)
boxplot(scores.NDB, las=2, col='gold', main='Principal components')

x11()
biplot(picia)

graphics.off()






