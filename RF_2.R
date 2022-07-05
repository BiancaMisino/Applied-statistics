#### RANDOM FOREST

library(randomForest)
library( car )
library( ellipse )
library( faraway )
library( leaps )
library(MASS)
library( GGally)
library(carData)
library(ISLR)
library(tidyverse)
library(caret)
library(clusterSim)
library(MLmetrics)
library(ROCR)
library(bootstrap)
library(stats)



# STEP 1.3
# Per-Patient analysis con e senza 'risposta radiologica'
# unica popolazione mista (milano + torino)
# cross-validation 5-fold 80/20

# da importare i dataset:
# DB_cliniche_core, DB_cliniche_core_margin

## CON RISPOSTA RADIOLOGICA

#NON funziona -> ERRORE SU: DB_cliniche_core[folds!=j,] 
# 5-fold CV
k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_core[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_core[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)


for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  modello_rf<-randomForest(TRG ~., data=DB_cliniche_core[folds!=j,], ntree = 200, proximity =TRUE)
  predicted <- predict(modello_rf, DB_cliniche_core[folds==j,])    
  
  cm <- confusionMatrix(as.factor(predicted), as.factor(DB_cliniche_core[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
}
accuracy.cv

mean(accuracy.cv)

#Funziona dividendo in train e test 
# cliniche+ core
set.seed(1)
training.samples <- DB_cliniche_core$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB_cliniche_core[training.samples, ]
test.data <- DB_cliniche_core[-training.samples, ]

x1 <- model.matrix(TRG ~., train.data)[,-1]
y1 <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)


modello_rf<-randomForest(x= x1,y =y1, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf


predicted <- predict(modello_rf)

confusionMatrix(as.factor(predicted), as.factor(y))

### cliniche + core + margin

x <- model.matrix(TRG ~., DB_cliniche_core_margin)[,-1]
y <- factor(DB_cliniche_core_margin$TRG) 


modello_rf<-randomForest(x= x,y =y, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf

predicted2 <- predict(modello_rf, x.test)

confusionMatrix(as.factor(predicted2), as.factor(y.test))

#SENZA RISPOSTA RADIOLOGICA
DB_cc_senzaRR <- DB_cliniche_core[,-c(14)]
DB_ccm_senzaRR <- DB_cliniche_core_margin[,-c(14)]
# cliniche+ core
set.seed(1)
training.samples <- DB_cc_senzaRR$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB_cc_senzaRR[training.samples, ]
test.data <- DB_cc_senzaRR[-training.samples, ]

x1 <- model.matrix(TRG ~., train.data)[,-1]
y1 <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)


modello_rf<-randomForest(x= x1,y =y1, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf

predicted <- predict(modello_rf, x.test)

confusionMatrix(as.factor(predicted), as.factor(y.test))

### cliniche + core + margin

set.seed(1)
training.samples <- DB_ccm_senzaRR$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB_ccm_senzaRR[training.samples, ]
test.data <- DB_ccm_senzaRR[-training.samples, ]

x1 <- model.matrix(TRG ~., train.data)[,-1]
y1 <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)


modello_rf<-randomForest(x= x1,y =y1, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf

predicted2 <- predict(modello_rf, x.test)

confusionMatrix(as.factor(predicted2), as.factor(y.test))

#OUTPUT IDENTICI PER TUTTI I CASI SOPRA
# Confusion Matrix and Statistics
# 
# Reference
# Prediction  1  3  5
# 1  0  0  0
# 3  0  0  0
# 5  5  5 23
# 
# Overall Statistics
# 
# Accuracy : 0.697           
# 95% CI : (0.5129, 0.8441)
# No Information Rate : 0.697           
# P-Value [Acc > NIR] : 0.5844          
#            
# Statistics by Class:
# 
#                      Class: 1 Class: 3 Class: 5
# Sensitivity            0.0000   0.0000    1.000
# Specificity            1.0000   1.0000    0.000
# Pos Pred Value            NaN      NaN    0.697
# Neg Pred Value         0.8485   0.8485      NaN
# Prevalence             0.1515   0.1515    0.697
# Detection Rate         0.0000   0.0000    0.697
# Detection Prevalence   0.0000   0.0000    1.000
# Balanced Accuracy      0.5000   0.5000    0.500

#######PARTE SOFI E BIANCA  #############

############# CASO TRG 1 VS 3 VS 5 ##########################

#CON RISPOSTA RADIOLOGICA
x <- model.matrix(TRG ~., DB_cliniche_core_margin)[,-1]
y <- factor(DB_cliniche_core_margin$TRG) 


modello_rf<-randomForest(x= x,y =y, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf

predicted <- predict(modello_rf)

confusionMatrix(as.factor(predicted), as.factor(y))


#SENZA RISPOSTA RADIOLOGICA
#LEVO RISPOSTA RADIOLOGICA
DB_senzaRR <- DB_cliniche_core_margin[,-c(14)]
x_senzaRR<- model.matrix(TRG ~., DB_senzaRR)[,-c(1)]
y <- factor(DB_senzaRR$TRG) 


modello_rf_senzaRR<-randomForest(x= x_senzaRR,y =y, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf_senzaRR

predicted2 <- predict(modello_rf_senzaRR)

confusionMatrix(as.factor(predicted2), as.factor(y))



#############  CASO  A 1-3 VS 5 #################################
#TGR:0->(1,3), 1->5
DB_A<-DB_cliniche_core_margin
DB_A[DB_A[,15]==1 | DB_A[,15]==3 ,15]<-0
DB_A[DB_A[,15]==5,15]<-1

#CON RISPOSTA RADIOLOGICA

x_A <- model.matrix(TRG ~., DB_A)[,-1]
y_A <- factor(DB_A$TRG) 


modello_rf_A<-randomForest(x= x_A,y =y_A, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf_A

predicted_A <- predict(modello_rf_A)

confusionMatrix(as.factor(predicted_A), as.factor(y_A))


#SENZA RISPOSTA RADIOLOGICA
#TGR:0->(1,3), 1->5
DB_AsenzaRR<-DB_senzaRR
DB_AsenzaRR[DB_AsenzaRR[,14]==1 | DB_AsenzaRR[,14]==3 ,14]<-0
DB_AsenzaRR[DB_AsenzaRR[,14]==5,14]<-1

x_AsenzaRR <- model.matrix(TRG ~., DB_AsenzaRR)[,-1]
y_AsenzaRR <- factor(DB_AsenzaRR$TRG) 


modello_rf_AsenzaRR<-randomForest(x= x_AsenzaRR,y =y_AsenzaRR, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf_AsenzaRR

predicted_ARR <- predict(modello_rf_AsenzaRR)

confusionMatrix(as.factor(predicted_ARR), as.factor(y_AsenzaRR))




# STEP 2.2
# Per-Lesion analysis 
# SOLO variabili RADIOMICHE
# analisi di tutte le lesioni: 378 (milano:319, torino:59)

#Ripulire tutto enviroment e caricare dataset con tutte le lesioni

Database_finale<-DATABASE_PER_LESION

Database_finale<-subset(Database_finale[,-c(1,2,3,5,11,24,25)])


attach(Database_finale)
##trasformo il sesso 1->M , 0->F
var1=which(Sex=='M') 
var2=which(Sex=='F')
var3=which(Sex=='f')
var4=which(Sex=='m')

Database_finale[var1,4]<-as.character(1)
Database_finale[var4,4]<-as.character(1)
Database_finale[var2,4]<-as.character(0)
Database_finale[var3,4]<-as.character(0)


Database_finale[,4]<-as.numeric(unlist(Database_finale[,4]))

###divido primary tumor site in left->0 rectum->1 right->2
#var1=which(Primary_tumor_site=='right') 
#var2=which(Primary_tumor_site=='left')
#var3=which(Primary_tumor_site=='rectum')

#Database_finale[var1,5]<-as.character(2)
#Database_finale[var2,5]<-as.character(0)
#Database_finale[var3,5]<-as.character(1)
#Database_finale[,5]<-as.numeric(unlist(Database_finale[,5]))

###creo variabile dummy per primary tumor site
Primary_tumor_site<-as.factor(Primary_tumor_site)
primary_tumor_site_left = ifelse(Primary_tumor_site == 'left', 1, 0)
primary_tumor_site_right = ifelse(Primary_tumor_site == 'right', 1, 0)
primary_tumor_site_rectum= ifelse(Primary_tumor_site=='rectum',1,0)

Database_finale<-data.frame(Database_finale,primary_tumor_site_left=primary_tumor_site_left,primary_tumor_site_right=primary_tumor_site_right,primary_tumor_site_rectum=primary_tumor_site_rectum)

##ora divido il sync in Met->0 e sync->1
var1=which(sync=='Met') 
var2=which(sync=='sync')
var3=which(sync=='met')

Database_finale[var1,6]<-as.character(0)
Database_finale[var2,6]<-as.character(1)
Database_finale[var3,6]<-as.character(0)
Database_finale[,6]<-as.numeric(unlist(Database_finale[,6]))

###ora divido le metastasi in classi
## 'Una'-> 0
## 'due-tre'->1
##'quattro-nove'->2
## '10+'->3
#var1=which(Numero_metastasi_classi=='Una') 
#var2=which(Numero_metastasi_classi=='due-tre')
#var3=which(Numero_metastasi_classi=='quattro-nove')
#var4=which(Numero_metastasi_classi=='10+') 

#Database_finale[var1,7]<-as.character(0)
#Database_finale[var2,7]<-as.character(1)
#Database_finale[var3,7]<-as.character(2)
#Database_finale[var4,7]<-as.character(3)
#Database_finale[,7]<-as.numeric(unlist(Database_finale[,7]))

#faccio dummy di numero metastasi classi

Numero_metastasi_classi<-as.factor(Numero_metastasi_classi)
numerometastasi_1 = ifelse(Numero_metastasi_classi == 'Una', 1, 0)
numerometastasi_2to3 = ifelse(Numero_metastasi_classi =='due-tre' , 1, 0)
numerometastasi_4to9 = ifelse(Numero_metastasi_classi=='quattro-nove',1,0)
numerometastasi_up10= ifelse(Numero_metastasi_classi=='10+',1,0)

Database_finale<-data.frame(Database_finale,numerometastasi_1=numerometastasi_1,numerometastasi_2to3=numerometastasi_2to3,numerometastasi_4to9=numerometastasi_4to9,numerometastasi_up10=numerometastasi_up10)



##diam in mm la rendo numerica
Database_finale[,8]<-as.numeric(unlist(Database_finale[,8]))

##binary la rendo numerica
Database_finale[,10]<-as.numeric(unlist(Database_finale[,10]))

##rendo numerica risp rad PR->0 SD->1

var1<-which(Risposta_radiologica=='PR')
var2<-which(Risposta_radiologica=='SD')

Database_finale[var1,20]<-as.character(0)
Database_finale[var2,20]<-as.character(1)
Database_finale[,20]<-as.numeric(unlist(Database_finale[,20]))

##rendo numerica M_CONVENTIONAL_HUmin
Database_finale[,70]<-as.numeric(unlist(Database_finale[,70]))
###traslo di 1 le linee di chemio-terapia
Database_finale[Database_finale[,18]==1,18] <-0
Database_finale[Database_finale[,18]==2,18] <-1

rm(var1)
rm(var2)
rm(var3)
rm(var4)
rm(numerometastasi_1)
rm(numerometastasi_2to3)
rm(numerometastasi_4to9)
rm(numerometastasi_up10)
rm(primary_tumor_site_left)
rm(primary_tumor_site_right)
rm(primary_tumor_site_rectum)

summary(Database_finale)

##non ci sono più character--> posso procedere con le analisi dividendo tra categoriche e non

###tolgo VOI che non ha senso
Database_finale<-Database_finale[,-c(1)]
##elenco var categoriche 1,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20

#Database_finale[,1]<-factor(`Interval _greater_30 days`)
#Database_finale[,3]<-factor(Sex)

#Database_finale[,5]<-factor(sync)
#Database_finale[,6]<-factor(Numero_metastasi_classi)
#Database_finale[,9]<-factor(bilater)
#Database_finale[,10]<-factor(KRAS,exclude=NA)
#Database_finale[,11]<-factor(NRAS,exclude=NA)
#Database_finale[,12]<-factor(BRAF,exclude=NA)
#Database_finale[,13]<-factor(OXALIPLATINO)
#Database_finale[,14]<-factor(IRINOTECAN)
#Database_finale[,15]<-factor(ANTI_VEGF)
#Database_finale[,16]<-factor(ANTI_EGFR)
#Database_finale[,17]<-factor(Linee_di_chemioterapia)
#Database_finale[,18]<-factor(greater_6_cicli)
#Database_finale[,19]<-factor(Risposta_radiologica)
#Database_finale[,20]<-factor(TRG)

#PREPARO I 3 DATASET SEPARATI

###da ora in poi il mio dataset è ok per procedere con le analisi

##incomincio solo cliniche
##tolgo le var che ho reso dummy

#####rimuovo le non  dummy
NRAS<-Database_finale[,3]
DB_cliniche<-Database_finale[,c(1:3,5,7:20,117:123)]


##rimuovo gli NA che non sappiamo trattare
DB_cliniche<-DB_cliniche[,-c(8,9,10)]

###tolgo quelle lin dipendenti-> nelle dummy le avevo aggiunte ma non aveva senso -> le ritolgo
DB_cliniche<-DB_cliniche[,-c(18,22)]

DB_cliniche[,2]<-scale(DB_cliniche[,2],center=FALSE,scale=TRUE)
DB_cliniche[,5]<-scale(DB_cliniche[,5],center=FALSE,scale=TRUE)
DB_cliniche[,6]<-scale(DB_cliniche[,6],center=FALSE,scale=TRUE)
#LE CLINICHE SONO PRONTE
summary(DB_cliniche)


DB_core<-Database_finale[,c(21:68)]
###tolgo le correlate del core
DB_core<-DB_core[,-c(5,6,7,9,11,18,21,23,27,29,30,31,33,34,35,38,41,44,45,46,47)]

DB_core<-data.frame(scale(DB_core,center=FALSE,scale=TRUE))
#IL CORE CHECK
summary(DB_core)

#IL MARGIN
DB_margin<-Database_finale[,c(69:116)]

###tolgo le correlate del core
DB_margin<-DB_margin[,-c(5,6,7,9,11,18,19,21,22,25,27,28,31,32,33,34,40,42,44,45,47)]

DB_margin<-data.frame(scale(DB_margin,center=FALSE,scale=TRUE))
#IL MARGIN CHECK
summary(DB_margin)

#DATASET FINALE:
DB_cliniche_core_margin<-data.frame(DB_cliniche,DB_core,DB_margin)

#mando Random

############# CASO TRG 1 VS 3 VS 5 ##########################

#CON RISPOSTA RADIOLOGICA
x <- model.matrix(TRG ~., DB_cliniche_core_margin)[,-1]
y <- factor(DB_cliniche_core_margin$TRG) 


modello_rf<-randomForest(x= x,y =y, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf

predicted <- predict(modello_rf)

confusionMatrix(as.factor(predicted), as.factor(y))


#SENZA RISPOSTA RADIOLOGICA
#LEVO RISPOSTA RADIOLOGICA
DB_senzaRR <- DB_cliniche_core_margin[,-c(14)]
x_senzaRR<- model.matrix(TRG ~., DB_senzaRR)[,-c(1)]
y <- factor(DB_senzaRR$TRG) 


modello_rf_senzaRR<-randomForest(x= x_senzaRR,y =y, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf_senzaRR

predicted2 <- predict(modello_rf_senzaRR)

confusionMatrix(as.factor(predicted2), as.factor(y))



#############  CASO  A 1-3 VS 5 #################################
#TGR:0->(1,3), 1->5
DB_A<-DB_cliniche_core_margin
DB_A[DB_A[,15]==1 | DB_A[,15]==3 ,15]<-0
DB_A[DB_A[,15]==5,15]<-1

#CON RISPOSTA RADIOLOGICA

x_A <- model.matrix(TRG ~., DB_A)[,-1]
y_A <- factor(DB_A$TRG) 


modello_rf_A<-randomForest(x= x_A,y =y_A, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf_A

predicted_A <- predict(modello_rf_A)

confusionMatrix(as.factor(predicted_A), as.factor(y_A))


#SENZA RISPOSTA RADIOLOGICA
#TGR:0->(1,3), 1->5
DB_AsenzaRR<-DB_senzaRR
DB_AsenzaRR[DB_AsenzaRR[,14]==1 | DB_AsenzaRR[,14]==3 ,14]<-0
DB_AsenzaRR[DB_AsenzaRR[,14]==5,14]<-1

x_AsenzaRR <- model.matrix(TRG ~., DB_AsenzaRR)[,-1]
y_AsenzaRR <- factor(DB_AsenzaRR$TRG) 


modello_rf_AsenzaRR<-randomForest(x= x_AsenzaRR,y =y_AsenzaRR, ntree = 200,importance=TRUE, proximity =TRUE)
modello_rf_AsenzaRR

predicted_ARR <- predict(modello_rf_AsenzaRR)

confusionMatrix(as.factor(predicted_ARR), as.factor(y_AsenzaRR))





