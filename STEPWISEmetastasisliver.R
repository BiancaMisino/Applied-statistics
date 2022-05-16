
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
library(dummy)
#levo le cliniche ridondanti

DB<-Database_finale[,-c(1,2,3,4,5,11,24,25)]

#rendo numeriche e categoriche le variabili:

#SEX
DB[DB[,3]=='M' | DB[,3]=='m',3]<-as.character(1)
DB[DB[,3]=='F' | DB[,3]=='f',3]<-as.character(0)
DB[,3]<-as.numeric(unlist(DB[,3]))
DB[,3]<-as.factor(DB$Sex)

#SYNC
DB[DB[,5]=='sync',5]<-as.character(1)
DB[DB[,5]=='met' | DB[,5]=='Met',5]<-as.character(0)
DB[,5]<-as.numeric(unlist(DB[,5]))
DB[,5]<-as.factor(DB$sync)


#diametro in mm
DB[,7]<-as.numeric(unlist(DB[,7]))

#bilater
DB[,9]<-as.numeric(unlist(DB[,9]))
DB[,9]<-as.factor(DB$bilater)

#primary tumor site
DB$Primary_tumor_site=as.factor(DB$Primary_tumor_site)
right_primary_tumor_site<-ifelse(DB$Primary_tumor_site=='right',1,0)
left_primary_tumor_site<-ifelse(DB$Primary_tumor_site=='left',1,0)

DB<-data.frame(DB,Left_primary_tumor_site=left_primary_tumor_site,Right_primary_tumor_site=right_primary_tumor_site)

#numero metastasi classi->dummy variables
DB$Numero_metastasi_classi=as.factor(DB$Numero_metastasi_classi)
classe1_Numero_metastasi<-ifelse(DB$Numero_metastasi_classi=='Una',1,0)
classe2_Numero_metastasi<-ifelse(DB$Numero_metastasi_classi=='due-tre',1,0)
classe3_Numero_metastasi<-ifelse(DB$Numero_metastasi_classi=='quattro-nove',1,0)

DB<-data.frame(DB,classe1_Numero_metastasi=classe1_Numero_metastasi,classe2_Numero_metastasi=classe2_Numero_metastasi,classe3_Numero_metastasi=classe3_Numero_metastasi)

#risposta radiologica
DB[DB[,19]=='PR',19]<-as.character(1)
DB[DB[,19]=='SD',19]<-as.character(0)
DB[,19]<-as.numeric(unlist(DB[,19]))

#M_HU_CONVENTIONAL
DB[,70]<-as.numeric(unlist(DB[,70]))

#interval greater 30 days
DB[,1]<-factor(DB$Interval._greater_30.days)

#linee di chemioterapia
DB[DB[,17]==1 ,17]<-0
DB[DB[,17]==2,17]<-1

DB[,17]<-as.factor(DB$Linee_di_chemioterapia)

#TGR:0->(1,3), 1->5
DB[DB[,20]==1 | DB[,20]==3 ,20]<-0
DB[DB[,20]==5,20]<-1



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

                  


############# STEPWISE CLINICHE+CORE CON RISPOSTA RADIOLOGICA##########################
library(leaps)
library(ISLR)



DB_core_cliniche<-DB[,c(1:68,117:121)]

#elimino variabili del core correlate :
DB_core_cliniche<-DB_core_cliniche[,-c(25,26,27,29,31,38,41,43,47,49,50,51,53,54,55,58,61,64,65,66,67)]

#elimino variabili NA e quelle delle dummy
DB_core_cliniche<-DB_core_cliniche[,-c(4,6,10,11,12)]

#standardizzo:
DB_core_cliniche[,c(2,5,6,16:42)]<-scale(DB_core_cliniche[,c(2,5,6,16:42)],center=FALSE)

regfit.step <- regsubsets(TRG~.,data=DB_core_cliniche,nvmax=47,method="forward")
summary(regfit.step)

which.max(summary(regfit.step)$adjr2)
coef(regfit.step,19)

max(summary(regfit.step)$adjr2)


x11(height=7,width=14)
par(mfrow=c(1,3))
plot(summary(regfit.step)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(summary(regfit.step)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(summary(regfit.step)$rss,xlab="Number of Variables",ylab="RSS",type="b")

x11()
plot(regfit.step,scale="r2",main="Stepwise Selection")

x11()
plot(regfit.step,scale="adjr2",main="Stepwise Selection")
reg.summary <- summary(regfit.step)
which.max(regfit.step.summary$adjr2)
coef(regfit.step,19)


############# STEPWISE CLINICHE+CORE SENZA RISPOSTA RADIOLOGICA##########################


library(leaps)
library(ISLR)



DB_core_cliniche2<-DB[,c(1:68,117:121)]

#elimino variabili del core correlate :
DB_core_cliniche2<-DB_core_cliniche2[,-c(25,26,27,29,31,38,41,43,47,49,50,51,53,54,55,58,61,64,65,66,67)]

#elimino variabili NA e quelle delle dummy
DB_core_cliniche2<-DB_core_cliniche2[,-c(4,6,10,11,12,19)]

#standardizzo:
DB_core_cliniche2[,c(2,5,6,15:41)]<-scale(DB_core_cliniche2[,c(2,5,6,15:41)],center=FALSE)

regfit.step2 <- regsubsets(TRG~.,data=DB_core_cliniche2,nvmax=47,method="forward")
summary(regfit.step2)

which.max(summary(regfit.step2)$adjr2)
coef(regfit.step2,26)

max(summary(regfit.step2)$adjr2)


x11(height=7,width=14)
par(mfrow=c(1,3))
plot(summary(regfit.step2)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(summary(regfit.step2)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(summary(regfit.step2)$rss,xlab="Number of Variables",ylab="RSS",type="b")

x11()
plot(regfit.step,scale="r2",main="Stepwise Selection")

x11()
plot(regfit.step,scale="adjr2",main="Stepwise Selection")
reg.summary <- summary(regfit.step)
which.max(reg.summary$adjr2)
coef(regfit.step,19)







