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
library(PRROC)

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

#M_HU_CONVENTIONAL
DB[,70]<-as.numeric(unlist(DB[,70]))

#linee di chemioterapia
DB[DB[,17]==1 ,17]<-0
DB[DB[,17]==2,17]<-1

########### CLINICHE+CORE  ###############

DB_core_cliniche<-DB[,c(1:68,117:121)]

#elimino variabili del core correlate :
DB_core_cliniche<-DB_core_cliniche[,-c(25,26,27,29,31,38,41,43,47,49,50,51,53,54,55,58,61,64,65,66,67)]

#elimino variabili NA e quelle delle dummy
DB_core_cliniche<-DB_core_cliniche[,-c(4,6,10,11,12)]

#standardizzo:
DB_core_cliniche[,c(2,5,6,16:42)]<-scale(DB_core_cliniche[,c(2,5,6,16:42)],center=FALSE)


######################### CASO A TRG:(1,3)->0 & (5)->1 #############################################                  
#TGR:0->(1,3), 1->5
DB_core_cliniche[DB_core_cliniche[,15]==1 | DB_core_cliniche[,15]==3 ,15]<-0
DB_core_cliniche[DB_core_cliniche[,15]==5,15]<-1
DB_core_cliniche$TRG<-factor(DB_core_cliniche$TRG,c(0,1),c(0,1))

############# LASSO CON RISPOSTA RADIOLOGICA##########################
library(leaps)
library(ISLR)


#variables selection
x <- model.matrix(TRG ~., DB_core_cliniche)[,-1]
y <- factor(DB_core_cliniche$TRG) 

library(glmnet)
set.seed(12)
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model<-glm(TRG ~ Interval._greater_30.days + sync + Diametro_in_mm 
              + OXALIPLATINO + ANTI_VEGF + ANTI_EGFR + Risposta_radiologica + HISTO_ExcessKurtosis
              + GLRLM_SRLGE + Left_primary_tumor_site + classe1_Numero_metastasi
             ,data=DB_core_cliniche,family = binomial)
CI<-confint(model, level = 0.95)
odds.ratio<-matrix(NA,1,model$rank - 1)
for (i in 1 :model$rank - 1){
  odds.ratio[i]=exp( coef( model )[ i ] )
}
model$coefficients
odds.ratio
CI
fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_core_cliniche$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_core_cliniche[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_core_cliniche[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG~ Interval._greater_30.days + sync + Diametro_in_mm 
             + OXALIPLATINO + ANTI_VEGF + ANTI_EGFR + Risposta_radiologica + HISTO_ExcessKurtosis
             + GLRLM_SRLGE + Left_primary_tumor_site + classe1_Numero_metastasi,
             data=DB_core_cliniche[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_core_cliniche[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_core_cliniche[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DB_core_cliniche[folds==j,]$TRG)
}


accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)

#ne tengo 10


############# LASSO SENZA RISPOSTA RADIOLOGICA##########################

library(leaps)
library(ISLR)

#levo la risposta radiologica
DB_core_cliniche2<-DB_core_cliniche[,-14]


#variables selection
x2 <- model.matrix(TRG ~., DB_core_cliniche2)[,-1]
y2 <- factor(DB_core_cliniche2$TRG) 


cv.lasso2 <- cv.glmnet(x2, y2, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso2,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model2<-glm(TRG ~ Interval._greater_30.days + sync + Diametro_in_mm 
                + CEA + OXALIPLATINO + ANTI_VEGF + ANTI_EGFR + Linee_di_chemioterapia
                + HISTO_ExcessKurtosis + GLRLM_SRLGE + Left_primary_tumor_site,
                 data=DB_core_cliniche2, family = binomial)
CI<-confint(model2, level = 0.95)
odds.ratio<-matrix(NA,1,model2$rank - 1)
for (i in 1 :model2$rank - 1){
  odds.ratio[i]=exp( coef( model2 )[ i ] )
}

model2$coefficients
odds.ratio
CI
fit2.2<-model2$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2.2, weights.class0=as.numeric(paste(DB_core_cliniche$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_core_cliniche[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_core_cliniche[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv2<-matrix(NA,1,k)
auc.cv2<-matrix(NA,1,k)

for (j in 1:k) {
  
  model2<-glm(TRG~ Interval._greater_30.days + sync + Diametro_in_mm 
              + CEA + OXALIPLATINO + ANTI_VEGF + ANTI_EGFR + Linee_di_chemioterapia
              + HISTO_ExcessKurtosis + GLRLM_SRLGE + Left_primary_tumor_site,
               data=DB_core_cliniche2[folds!=j,], family=binomial)
  probabilities <- predict(model2,newdata=DB_core_cliniche2[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_core_cliniche2[folds==j,]$TRG))
  accuracy.cv2[j]<-cm$overall[1]
  auc.cv2[j]<-AUC(predicted,DB_core_cliniche2[folds==j,]$TRG)
}


accuracy.cv2
auc.cv2

mean(accuracy.cv2)
mean(auc.cv2)



######################### CASO A TRG:(1)->0 & (3,5)->1 #############################################                  
#TGR:0->(1,3), 1->5
DB_core_cliniche[DB[,20]==1 ,15]<-0
DB_core_cliniche[DB[,20]==5 | DB[,20]==3,15]<-1
DB_core_cliniche$TRG<-factor(DB_core_cliniche$TRG,c(0,1),c(0,1))


############# LASSO CON RISPOSTA RADIOLOGICA##########################
library(leaps)
library(ISLR)


#variables selection
x <- model.matrix(TRG ~., DB_core_cliniche)[,-1]
y <- factor(DB_core_cliniche$TRG) 

library(glmnet)

cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model<-glm(TRG ~ Risposta_radiologica + GLRLM_SRLGE + 
                 GLZLM_LGZE + left_primary_tumor_site,
             data=DB_core_cliniche,family = binomial)
CI<-confint(model, level = 0.95)
odds.ratio<-matrix(NA,1,model$rank - 1)
for (i in 1 :model$rank - 1){
  odds.ratio[i]=exp( coef( model )[ i ] )
}
model$coefficients
odds.ratio
CI

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_core_cliniche[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_core_cliniche[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG~ Risposta_radiologica + GLRLM_SRLGE + 
               GLZLM_LGZE + Left_primary_tumor_site,
               data=DB_core_cliniche[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_core_cliniche[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_core_cliniche[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(predicted,DB_core_cliniche[folds==j,]$TRG)
}


accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)

#ne tengo 10


############# LASSO SENZA RISPOSTA RADIOLOGICA##########################

library(leaps)
library(ISLR)

#levo la risposta radiologica
DB_core_cliniche2<-DB_core_cliniche[,-14]


#variables selection
x2 <- model.matrix(TRG ~., DB_core_cliniche2)[,-1]
y2 <- factor(DB_core_cliniche2$TRG) 


cv.lasso2 <- cv.glmnet(x2, y2, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso2,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model2<-glm(TRG ~ GLRLM_SRLGE + GLZLM_LGZE + Left_primary_tumor_site 
             ,data=DB_core_cliniche2,family = binomial)
CI<-confint(model2, level = 0.95)
odds.ratio<-matrix(NA,1,model2$rank - 1)
for (i in 1 :model2$rank - 1){
  odds.ratio[i]=exp( coef( model2 )[ i ] )
}

model2$coefficients
odds.ratio
CI

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_core_cliniche[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_core_cliniche[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv2<-matrix(NA,1,k)
auc.cv2<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG~ GLRLM_SRLGE + GLZLM_LGZE + Left_primary_tumor_site ,
               data=DB_core_cliniche2[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_core_cliniche2[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_core_cliniche2[folds==j,]$TRG))
  accuracy.cv2[j]<-cm$overall[1]
  auc.cv2[j]<-AUC(predicted,DB_core_cliniche2[folds==j,]$TRG)
}


accuracy.cv2
auc.cv2

mean(accuracy.cv2)
mean(auc.cv2)


