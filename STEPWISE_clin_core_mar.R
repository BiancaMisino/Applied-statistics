setwd("C:/Users/sarar/Desktop/POLI/APPLIED STATISTIC/Project/Applied-statistics-main3")

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
library(glmnet)
library(tidyverse)
library(caret)
library(clusterSim)
library(MLmetrics)
library(ROCR)
library(readxl)
library(leaps)

#importo dataset
Database_finale <- read_excel("Database_finale.xlsx")

#table(DB$TRG)
#1   3   5 
#24  28 117 

#collinearità radiomiche margin
Database_finale<-Database_finale[,-c(81,82,83,85,87,94,95,97,98,101,103,104,107,108,109,110,116,118,120,121,123)]
#collinearità radiomiche core
Database_finale<-Database_finale[,-c(33,34,35,37,39,40,46,49,51,55,58,59,60,61,62,66,69,72,73,74,75)]

##### SCELTA DEL DATASET
#DB_cliniche <- Database_finale[,-c(29:82)]
#DB_cl_core <- Database_finale[,-c(56:82)]
DB_cl_core_mar <- Database_finale

#levo le cliniche ridondanti 
#DB<-DB_cliniche[,-c(1,2,3,4,5,11,24,25)]
#DB<-DB_cl_core[,-c(1,2,3,4,5,11,24,25)]
DB <- DB_cl_core_mar[,-c(1,2,3,4,5,11,24,25)]


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

#HU_CONVENTIONAL
DB[,48]<-as.numeric(unlist(DB[,48]))

#linee di chemioterapia
#linee di chemioterapia
DB[DB[,17]==1 ,17]<-0
DB[DB[,17]==2,17]<-1

#elimino le variabili con NA
DB<-DB[,-c(10,11,12)]

# sistemo dummy
DBP <- model.matrix(TRG ~., DB)[,-1]
DBP <- cbind(DBP,DB[,17])

#standardizzo:
DBP[,-c(1,3:6,9:14, 15,16,71)]<-scale(DBP[,-c(1,3:6,9:14, 15,16,71)],center=FALSE)


########################################################################################
######################### CASO A TRG:(1,3)->0 & (5)->1 #############################################                  

#TGR:0->(1,3), 1->5
DBP[DBP[,71]==1 | DBP[,71]==3 ,71]<-0
DBP[DBP[,71]==5,71]<-1
DBP$TRG<-factor(DBP$TRG,c(0,1),c(0,1))

# calcolo la percentuale di TRG=1 sul totale dei dati
n<- table(DBP$TRG)
n1<-n[names(n)==1]
ratio_of_1<-round(n1/169,digits = 2)
ratio_of_1


############# STEPWISE CON RISPOSTA RADIOLOGICA##########################

#variables selection
reg0<-glm(TRG~1,data=DBP,family = binomial)
reg1<-glm(TRG~.,data=DBP,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$formula


#model selection in 10-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DBP[1:127,]),replace=TRUE)
folds_to <- sample(1:k,nrow(DBP[128:169,]),replace=TRUE)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to
table(folds)

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG ~ Risposta_radiologica + ANTI_VEGF + CEA +
               `\`M_SHAPE_Sphericity (only for 3D ROI (nz>1)\`` + 
               M_HISTO_Entropy_log2 + GLRLM_SRLGE + M_GLRLM_LGRE +
               M_CONVENTIONAL_HUmean + OXALIPLATINO + M_GLZLM_SZE + 
               M_GLZLM_GLNU + Diametro_in_mm + `\`SHAPE_Compacity only for 3D ROI (nz>1)\``,
             data=DBP[folds!=j,], family=binomial)
  
  probabilities <- predict(model,newdata=DBP[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > ratio_of_1, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DBP[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DBP[folds==j,]$TRG)
}

# modello:
best.fit$formula

accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)


############# STEPWISE SENZA RISPOSTA RADIOLOGICA##########################

# tolgo la risposta radiologica
DBP<-DBP[,-c(16)]

#variables selection
reg0<-glm(TRG~1,data=DBP,family = binomial)
reg1<-glm(TRG~.,data=DBP,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$formula

#model selection in 10-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(123)
folds_mi <- sample(1:k,nrow(DBP[1:127,]),replace=TRUE)
folds_to <- sample(1:k,nrow(DBP[128:169,]),replace=TRUE)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to
table(folds)

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG ~ CEA + ANTI_VEGF + `\`M_SHAPE_Sphericity (only for 3D ROI (nz>1)\`` + 
               Primary_tumor_site2 + Diametro_in_mm + GLRLM_SRLGE + M_GLZLM_SZE + 
               `\`SHAPE_Volume (mL)\`` + M_HISTO_Entropy_log2 + M_GLRLM_LGRE + 
               M_GLRLM_SRHGE + OXALIPLATINO + NGLDM_Coarseness,
             data=DBP[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DBP[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > ratio_of_1, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DBP[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DBP[folds==j,]$TRG)
}

#modello:
best.fit$formula

accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)



########################################################################################
######################### CASO A TRG:(1)->0 & (3,5)->1 #############################################                  

#TGR:0->(1), 1->(3,5)
DBP[DBP[,71]==1,71]<-0
DBP[DBP[,71]==5| DBP[,71]==3,71]<-1
DBP$TRG<-factor(DBP$TRG,c(0,1),c(0,1))

# calcolo la percentuale di TRG=1 sul totale dei dati
n<- table(DBP$TRG)
n1<-n[names(n)==1]
ratio_of_1<-round(n1/169,digits = 2)
ratio_of_1


############# STEPWISE CON RISPOSTA RADIOLOGICA##########################


#variables selection
reg0<-glm(TRG~1,data=DBP,family = binomial)
reg1<-glm(TRG~.,data=DBP,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$formula

#model selection in 10-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(123)
folds_mi <- sample(1:k,nrow(DBP[1:127,]),replace=TRUE)
folds_to <- sample(1:k,nrow(DBP[128:169,]),replace=TRUE)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to
table(folds)

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG ~ Risposta_radiologica + GLZLM_LGZE + M_CONVENTIONAL_HUstd + 
               Linee_di_chemioterapia + GLZLM_ZP,
             data=DBP[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DBP[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > ratio_of_1, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DBP[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DBP[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)


############# STEPWISE SENZA RISPOSTA RADIOLOGICA##########################

# tolgo la risposta radiologica
DBP<-DBP[,-c(16)]

#variables selection
reg0<-glm(TRG~1,data=DBP,family = binomial)
reg1<-glm(TRG~.,data=DBP,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$formula

#model selection in 10-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(13)
folds_mi <- sample(1:k,nrow(DBP[1:127,]),replace=TRUE)
folds_to <- sample(1:k,nrow(DBP[128:169,]),replace=TRUE)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to
table(folds)

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG ~ GLRLM_SRLGE + Linee_di_chemioterapia + 
               `\`GLCM_Homogeneity (=Inverse difference)\`` + 
               OXALIPLATINO + greater_6_cicli + CEA,
             data=DBP[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DBP[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > ratio_of_1, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DBP[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DBP[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)

