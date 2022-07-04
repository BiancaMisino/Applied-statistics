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
  x_train <- model.matrix(TRG ~., data=DB_cliniche_core[folds!=j,])[,-1]
  y_train <- factor(DB_cliniche_core[folds!=j,]$TRG) 
  x_test <- model.matrix(TRG ~., DB_cliniche_core[folds==j,])[,-1]
  y_test <- factor(DB_cliniche_core[folds==j,]$TRG) 
  
  modello_rf<-randomForest(TRG ~., data=DB_cliniche_core[folds!=j,], ntree = 200, proximity =TRUE)
  predicted <- predict(modello_rf, x_test)    
  
  cm <- confusionMatrix(as.factor(predicted), as.factor(y_test))
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

predicted <- predict(modello_rf, x.test)

confusionMatrix(as.factor(predicted), as.factor(y.test))

### cliniche + core + margin

set.seed(1)
training.samples <- DB_cliniche_core_margin$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB_cliniche_core_margin[training.samples, ]
test.data <- DB_cliniche_core_margin[-training.samples, ]

x1 <- model.matrix(TRG ~., train.data)[,-1]
y1 <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)


modello_rf<-randomForest(x= x1,y =y1, ntree = 200,importance=TRUE, proximity =TRUE)
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


# STEP 2.2
# Per-Lesion analysis 
# SOLO variabili RADIOMICHE
# analisi di tutte le lesioni: 378 (milano:319, torino:59)




