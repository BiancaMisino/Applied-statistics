
#import the library
library( car )
library( ellipse )
library( faraway )
library( leaps )
library(MASS)
library( GGally)
library(rgl)
library(carData)
library(mvtnorm)
library(ISLR)

#levo le cliniche ridondanti

DB<-Database_finale[,-c(1,2,3,4,5,11,24,25)]

#rendo numeriche e categoriche le variabili:

#SEX
DB[DB[,3]=='M' | DB[,3]=='m',3]<-as.character(1)
DB[DB[,3]=='F' | DB[,3]=='f',3]<-as.character(0)
DB[,3]<-as.numeric(unlist(DB[,3]))
DB[,3]<-factor(DB$Sex)

#SYNC
DB[DB[,5]=='sync',5]<-as.character(1)
DB[DB[,5]=='met' | DB[,5]=='Met',5]<-as.character(0)
DB[,5]<-as.numeric(unlist(DB[,5]))
DB[,5]<-factor(DB$sync)


#diametro in mm
DB[,7]<-as.numeric(unlist(DB[,7]))

#bilater
DB[,9]<-as.numeric(unlist(DB[,9]))
DB[,9]<-factor(DB$bilater)

#primary tumor site
var1=which(DB$Primary_tumor_site=='right') 
var2=which(DB$Primary_tumor_site=='left')
var3=which(DB$Primary_tumor_site=='rectum')

DB[var1,4]<-as.character(2)
DB[var2,4]<-as.character(0)
DB[var3,4]<-as.character(1)
DB[,4]<-as.numeric(unlist(DB[,5]))
DB[,4]<-factor(DB$Primary_tumor_site)

#numero metastasi classi
var1=which(DB$Numero_metastasi_classi=='Una') 
var2=which(DB$Numero_metastasi_classi=='due-tre')
var3=which(DB$Numero_metastasi_classi=='quattro-nove')
var4=which(DB$Numero_metastasi_classi=='10+') 

DB[var1,6]<-as.character(0)
DB[var2,6]<-as.character(1)
DB[var3,6]<-as.character(2)
DB[var4,6]<-as.character(3)
DB[,6]<-as.numeric(unlist(Database_finale[,6]))
DB[,6]<-factor(DB$Numero_metastasi_classi)

#risposta radiologica
DB[DB[,19]=='PR',19]<-as.character(1)
DB[DB[,19]=='SD',19]<-as.character(0)
DB[,19]<-as.numeric(unlist(DB[,19]))

#M_HU_CONVENTIONAL
DB[,70]<-as.numeric(unlist(DB[,70]))

#interval grater 30 days
DB[,1]<-factor(DB$`Interval _greater_30 days`)

#variabili da 10 a 20
DB[,10]<-factor(DB$KRAS)
DB[,11]<-factor(DB$NRAS)
DB[,12]<-factor(DB$BRAF)
DB[,13]<-factor(DB$OXALIPLATINO)
DB[,14]<-factor(DB$IRINOTECAN)
DB[,15]<-factor(DB$ANTI_VEGF)
DB[,16]<-factor(DB$ANTI_EGFR)
DB[,17]<-factor(DB$Linee_di_chemioterapia)
DB[,18]<-factor(DB$greater_6_cicli)
DB[,19]<-factor(DB$Risposta_radiologica)
DB[,20]<-factor(DB$TRG)
                  



























help(DB_cli_cor)
names(DB_cli_cor)
dim(DB_cli_cor)

unlist(DB_cli_cor)
as.double(DB_cli_cor)



###  Subset Selection Methods
help(regsubsets)

# Forward and Backward Stepwise Selection
regfit.fwd <- regsubsets(TRG~.,data=DB_cli_cor,method="forward")
summary(regfit.fwd)




x11(height=7,width=14)
par(mfrow=c(1,3))
plot(summary(regfit.fwd)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(summary(regfit.fwd)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(summary(regfit.fwd)$rss,xlab="Number of Variables",ylab="RSS",type="b")

x11()
plot(regfit.fwd,scale="r2",main="Forward Stepwise Selection")

x11()
plot(regfit.fwd,scale="adjr2",main="Forward Stepwise Selection")



regfit.bwd <- regsubsets(Salary~.,data=Hitters,nvmax=19,method="backward")
summary(regfit.bwd)

x11(height=7,width=14)
par(mfrow=c(1,3))
plot(summary(regfit.bwd)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(summary(regfit.bwd)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(summary(regfit.bwd)$rss,xlab="Number of Variables",ylab="RSS",type="b")

x11()
plot(regfit.bwd,scale="r2",main="Backward Stepwise Selection")

x11()
plot(regfit.bwd,scale="adjr2",main="Backward Stepwise Selection")


coef(regfit.full,7) # Exhaustive search
coef(regfit.fwd,7)  # Forward Stepwise Selection
coef(regfit.bwd,7)  # Backward Stepwise Selection

graphics.off()