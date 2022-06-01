setwd("C:/Users/sarar/Desktop/POLI/APPLIED STATISTIC/Project/Applied-statistics-main5")


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
library(ROCR)
library(glmnet)
library(tidyverse)
library(MLmetrics)
library(caret)
library(readxl)
library(glmnet)
library(leaps)
library(ISLR)
library(PRROC)


#importo dataset
Database_finale <- read_excel("Database_finale.xlsx")

#levo le cliniche ridondanti

DB<-Database_finale[,-c(1,2,3,4,5,11,24,25)]

#rendo numeriche e categoriche le variabili:

#SEX
DB[DB[,3]=='M' | DB[,3]=='m',3]<-as.character(1)
DB[DB[,3]=='F' | DB[,3]=='f',3]<-as.character(0)
DB[,3]<-as.numeric(unlist(DB[,3]))


#SYNC
DB[DB[,5]=='sync',5]<-as.character(1)
DB[DB[,5]=='met' | DB[,5]=='Met',5]<-as.character(0)
DB[,5]<-as.numeric(unlist(DB[,5]))



#diametro in mm
DB[,7]<-as.numeric(unlist(DB[,7]))

#bilater
DB[,9]<-as.numeric(unlist(DB[,9]))


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
DB[DB[,19]=='PR',19]<-as.character(0)
DB[DB[,19]=='SD',19]<-as.character(1)
DB[,19]<-as.numeric(unlist(DB[,19]))

#M_HU_CONVENTIONALmean
DB[,70]<-as.numeric(unlist(DB[,70]))

#M_HU_CONVENTIONALmin
DB[,69]<-as.numeric(unlist(DB[,69]))

#linee di chemioterapia
DB[DB[,17]==1 ,17]<-0
DB[DB[,17]==2,17]<-1

#elimino variabili del core correlate :
DB<-DB[,-c(25,26,27,29,31,38,41,43,47,49,50,51,53,54,55,58,61,64,65,66,67)]

#elimino le variabili del margin correlate:
DB<-DB[,-c(52,53,54,56,58,65,66,68,69,72,74,75,78,79,80,81,87,89,91,93,94)]

#elimino variabili NA e quelle delle dummy
DB<-DB[,-c(4,6,10,11,12)]

#standardizzo:
DB[,c(2,5,6,16:69)]<-scale(DB[,c(2,5,6,16:69)],center=FALSE)


##########################################################################################
########### CLINICHE #############################################################

DB_cliniche<-DB[,c(1:15,70:74)]


######################### CASO A TRG:(1,3)->0 & (5)->1 #############################################                  
#TGR:0->(1,3), 1->5
DB_cliniche[DB_cliniche[,15]==1 | DB_cliniche[,15]==3 ,15]<-0
DB_cliniche[DB_cliniche[,15]==5,15]<-1
DB_cliniche$TRG<-factor(DB_cliniche$TRG,c(0,1),c(0,1))

# calcolo la percentuale di TRG=1 sul totale dei dati
n<- table(DB_cliniche$TRG)
n1<-n[names(n)==1]
ratio_of_1<-round(n1/169,digits = 2)
ratio_of_1

#############  CON RISPOSTA RADIOLOGICA##########################

#variables selection
x <- model.matrix(TRG ~., DB_cliniche)[,-1]
y <- factor(DB_cliniche$TRG) 


cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
cv.lasso$lambda.min
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables


model<-glm(TRG ~ Interval._greater_30.days+Sex+sync+Diametro_in_mm+CEA+bilater+OXALIPLATINO+
                 IRINOTECAN+ANTI_VEGF+ANTI_EGFR+Linee_di_chemioterapia+greater_6_cicli+Risposta_radiologica+
                 Left_primary_tumor_site+classe1_Numero_metastasi,
                 data=DB_cliniche,family = binomial)


#intervalli di confidenza dei beta
CI<-confint(model, level = 0.95)

#Odds ratio
odds.ratio<-matrix(NA,1,model$rank-1)
odds.ratio<-exp(model$coefficients)

model$coefficients
CI
odds.ratio

#curva roc
fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)


k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG~ Interval._greater_30.days+Sex+sync+Diametro_in_mm+CEA+bilater+OXALIPLATINO+
               IRINOTECAN+ANTI_VEGF+ANTI_EGFR+Linee_di_chemioterapia+greater_6_cicli+Risposta_radiologica+
               Left_primary_tumor_site+classe1_Numero_metastasi,
             data=DB_cliniche[folds!=j,], family=binomial)
  
  probabilities <- predict(model,newdata=DB_cliniche[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > ratio_of_1, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DB_cliniche[folds==j,]$TRG)
}


accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)



#############  SENZA RISPOSTA RADIOLOGICA##########################

# tolgo la risposta radiologica, colonna 14
DB_cliniche<-DB_cliniche[,-c(14)]

#variables selection
x <- model.matrix(TRG ~., DB_cliniche)[,-1]
y <- factor(DB_cliniche$TRG) 


cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
cv.lasso$lambda.min
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model<-glm(TRG ~ Interval._greater_30.days+ Age+sync+ Diametro_in_mm + CEA + bilater+                              
             OXALIPLATINO+IRINOTECAN+ANTI_VEGF+ANTI_EGFR+Linee_di_chemioterapia+greater_6_cicli+
             Left_primary_tumor_site+Right_primary_tumor_site+classe1_Numero_metastasi+
             classe2_Numero_metastasi+classe3_Numero_metastasi,
                 data=DB_cliniche,family = binomial)
#curva roc
fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

#intervalli di confidenza dei beta
CI<-confint(model, level = 0.95)

#Odds ratio
odds.ratio<-matrix(NA,1,model$rank-1)
odds.ratio<-exp(model$coefficients)

model$coefficients
CI
odds.ratio


k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG~ Interval._greater_30.days+ Age+sync+ Diametro_in_mm + CEA + bilater+                              
               OXALIPLATINO+IRINOTECAN+ANTI_VEGF+ANTI_EGFR+Linee_di_chemioterapia+greater_6_cicli+
               Left_primary_tumor_site+Right_primary_tumor_site+classe1_Numero_metastasi+
               classe2_Numero_metastasi+classe3_Numero_metastasi,
             data=DB_cliniche[folds!=j,], family=binomial)
  
  probabilities <- predict(model,newdata=DB_cliniche[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > ratio_of_1, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DB_cliniche[folds==j,]$TRG)
}


accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)



######################### CASO A TRG:(1)->0 & (3,5)->1 #############################################                  
#TGR:0->(1), 1->(3,5)
DB_cliniche[DB_cliniche[,15]==1,15]<-0
DB_cliniche[DB_cliniche[,15]==5| DB_cliniche[,15]==3 ,15]<-1
DB_cliniche$TRG<-factor(DB_cliniche$TRG,c(0,1),c(0,1))

# calcolo la percentuale di TRG=1 sul totale dei dati
n<- table(DB_cliniche$TRG)
n1<-n[names(n)==1]
ratio_of_1<-round(n1/169,digits = 2)
ratio_of_1

#############  CON RISPOSTA RADIOLOGICA##########################

#variables selection
x <- model.matrix(TRG ~., DB_cliniche)[,-1]
y <- factor(DB_cliniche$TRG) 


cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
cv.lasso$lambda.min
x11()
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min

c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables


model<-glm(TRG ~ Interval._greater_30.days+Age+Sex+sync+Diametro_in_mm+bilater+OXALIPLATINO+
             IRINOTECAN+ANTI_VEGF+ANTI_EGFR+Linee_di_chemioterapia+greater_6_cicli+
             Risposta_radiologica+Left_primary_tumor_site+Right_primary_tumor_site+
             classe2_Numero_metastasi+classe3_Numero_metastasi,
           data=DB_cliniche,family = binomial)


#intervalli di confidenza dei beta
CI<-confint(model, level = 0.95)

#Odds ratio
odds.ratio<-matrix(NA,1,model$rank-1)
odds.ratio<-exp(model$coefficients)

model$coefficients
CI
odds.ratio


k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG~ Interval._greater_30.days+Age+Sex+sync+Diametro_in_mm+bilater+OXALIPLATINO+
               IRINOTECAN+ANTI_VEGF+ANTI_EGFR+Linee_di_chemioterapia+greater_6_cicli+
               Risposta_radiologica+Left_primary_tumor_site+Right_primary_tumor_site+
               classe2_Numero_metastasi+classe3_Numero_metastasi,
             data=DB_cliniche[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > ratio_of_1, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DB_cliniche[folds==j,]$TRG)
}


accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)



#############  SENZA RISPOSTA RADIOLOGICA##########################

# tolgo la risposta radiologica, colonna 14
DB_cliniche<-DB_cliniche[,-c(14)]

#variables selection
x <- model.matrix(TRG ~., DB_cliniche)[,-1]
y <- factor(DB_cliniche$TRG) 


cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
cv.lasso$lambda.min
x11()
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)

c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model<-glm(TRG ~ Interval._greater_30.days+Age+Sex+sync+Diametro_in_mm+
             CEA+bilater+OXALIPLATINO+IRINOTECAN+ANTI_VEGF+ANTI_EGFR+
             Linee_di_chemioterapia+greater_6_cicli+Left_primary_tumor_site+
             Right_primary_tumor_site+classe2_Numero_metastasi+classe3_Numero_metastasi,
           data=DB_cliniche,family = binomial)


#intervalli di confidenza dei beta
CI<-confint(model, level = 0.95)

#Odds ratio
odds.ratio<-matrix(NA,1,model$rank-1)
odds.ratio<-exp(model$coefficients)

model$coefficients
CI
odds.ratio


k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG~ Interval._greater_30.days+Age+Sex+sync+Diametro_in_mm+
               CEA+bilater+OXALIPLATINO+IRINOTECAN+ANTI_VEGF+ANTI_EGFR+
               Linee_di_chemioterapia+greater_6_cicli+Left_primary_tumor_site+
               Right_primary_tumor_site+classe2_Numero_metastasi+classe3_Numero_metastasi ,
             data=DB_cliniche[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > ratio_of_1, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DB_cliniche[folds==j,]$TRG)
}


accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)


