## LASSO
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

#table(DB$TRG)
#1   3   5 
#24  28 117 

#ELIMINO VARIABILI COLLINEARI
#collinearità radiomiche margin
Database_finale<-Database_finale[,-c(81,82,83,85,87,94,95,97,98,101,103,104,107,108,109,110,116,118,120,121,123)]
#collinearità radiomiche core
Database_finale<-Database_finale[,-c(33,34,35,37,39,40,46,49,51,55,58,59,60,61,62,66,69,72,73,74,75)]

#levo le cliniche ridondanti
DB<-Database_finale[,-c(1,2,3,4,5,11,24,25)]

#PULIZIA CLINICHE
#rendo numeriche le variabili categoriche:
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
# DB[,20]<-factor(DB$TRG) da fare dopo la categorizzazione

#M_CONVENTIONAL_HU_MIN (radiomica_margin)
DB[,48]<-as.numeric(unlist(DB[,48]))

#DB <- DB[-45,] (outlier CEA)

#elimino le variabili con NA
DB<-DB[,-c(10,11,12)]

#write.table(DB, "DB_pulito.csv")

# SPLIT PER VALUTARE I MODELLI SUL SOLITO TEST
# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
DB_train  <- DB[training.samples, ] # 1->21, 3->21, 5->94 
DB_test <- DB[-training.samples, ]


####################################################################################
################CASO Dicotomico A (CLINICHE)
DB_cA <- DB_train[,-c(18:71)]
DBTEST_cA <- DB_test[,-c(18:71)]
################TRG ={0,1} (1+3,5) 

DB_cA[DB_cA[,17]<=3,17] <-0
DB_cA[DB_cA[,17]>=5,17] <-1
DBTEST_cA[DBTEST_cA[,17]<=3,17] <-0
DBTEST_cA[DBTEST_cA[,17]>=5,17] <-1

#MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_cA<-DB_cA[,-c(16)]
#DBTEST_cA<-DBTEST_cA[,-c(16)]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., DB_cA)[,-1]
y <- factor(DB_cA$TRG) # 0->42, 1->94 , Proportion = 0.69
x.test <- model.matrix(TRG ~., DBTEST_cA)[,-1]
y.test <- factor(DBTEST_cA$TRG)

# default alpha = 1 -> LASSO, alpha = 0 -> RIDGE, alpha (0,1) -> ELASTIC NET
# Set lambda via cross validation
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "response")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
predicted <- ifelse(probabilities > 0.69, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 

#ROC curve
#roc_pred <- prediction(predictions = probabilities , labels = y.test)
#roc_perf <- performance(roc_pred , "tpr" , "fpr")
#plot(roc_perf,colorize = TRUE,
#     print.cutoffs.at= seq(0,1,0.05),
#     text.adj=c(-0.2,1.7))
#auc_ROCR <- performance(roc_pred, measure = "auc")
#auc_ROCR@y.values[[1]]
                       

################CASO Dicotomico B 
DB_cB <- DB_train[,-c(18:71)]
DBTEST_cB <- DB_test[,-c(18:71)]
################TRG ={0,1} (1,3+5) 
DB_cB[DB_cB[,17]<=1,17] <-0
DB_cB[DB_cB[,17]>=3,17] <-1
DBTEST_cB[DBTEST_cB[,17]<=1,17] <-0
DBTEST_cB[DBTEST_cB[,17]>=3,17] <-1

##MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_cB<-DB_cB[,-c(16)]
#DBTEST_cB<-DBTEST_cB[,-c(16)]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., DB_cB)[,-1]
y <- factor(DB_cB$TRG) # 0->21, 1->115 , Proportion = 0.84
x.test <- model.matrix(TRG ~., DBTEST_cB)[,-1]
y.test <- factor(DBTEST_cB$TRG)

# Set lambda via cross validation
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "response")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
predicted <- ifelse(probabilities > 0.84, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test)


################CASO multinomiale
DB_cMUL <- DB_train[,-c(18:71)]
DBTEST_cMUL <- DB_test[,-c(18:71)]
####TRG ={1,3,5}   1->21, 3->21, 5->94 

#type.multinomial() allows the usage of a grouped lasso penalty
#(q=2->is a grouped-lasso penalty on all the K coefficients for a particular variable, which makes them 
#all be zero or nonzero together) if type.multinomial = "grouped". The default is type.multinomial = "ungrouped" 
#(q=1 -> lasso penalty on each of the parameters, called partial Newton algorithm)

#MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_cMUL<-DB[,-c(16)]
#DBTEST_cMUL<-DBTEST_cMUL[,-c(16)]

#MODELLO
# Divide predicotors and target
x <- model.matrix(TRG ~., DB_cMUL)[,-1]
y <- factor(DB_cMUL$TRG) 
x.test <- model.matrix(TRG ~., DBTEST_cMUL)[,-1]
y.test <- factor(DBTEST_cMUL$TRG)

# Set lambda via cross validation
#type measure:deviance(default), class, mae,mse
#type.multinomial = grouped -> results are worse or similar
cv.lasso <- cv.glmnet(x, y, family = "multinomial")#,type.measure = "class")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
#probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
#predicted <- ifelse(probabilities > ?, "1","0")
predicted <-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")

#Model Perfomance
mean(predicted == y.test)

#################################################################
#################CASO Dicotomico A (CLINCHE + RADIOMICHE_CORE)
DB_ccA <- DB_train[,-c(45:71)]
DBTEST_ccA <- DB_test[,-c(45:71)]
################TRG ={0,1} (1+3,5)
DB_ccA[DB_ccA[,17]<=3,17] <-0
DB_ccA[DB_ccA[,17]>=5,17] <-1
DBTEST_ccA[DBTEST_ccA[,17]<=3,17] <-0
DBTEST_ccA[DBTEST_ccA[,17]>=5,17] <-1


#STANDARDIZZAZIONE MANUALE -> nel caso cambiare "standardize = FALSE" riga 278
#standardize Radiomics AND età, diam_mm, cea
#DB_ccA[,c(2,7,8,18:44)] <- scale(DB_ccA[,c(2,7,8,18:44)])

#MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_ccA<-DB_ccA[,-c(16)]
#DBTEST_ccA<-DBTEST_ccA[,-c(16)]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., DB_ccA)[,-1]
y <- factor(DB_ccA$TRG) # 0->42, 1->94 , Proportion = 0.69
x.test <- model.matrix(TRG ~., DBTEST_ccA)[,-1]
y.test <- factor(DBTEST_ccA$TRG)

# Set lambda via cross validation
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")#,standardize = FALSE)
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "response")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
predicted <- ifelse(probabilities > 0.69, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 

################CASO Dicotomico B
DB_ccB <- DB_train[,-c(45:71)]
DBTEST_ccB <- DB_test[,-c(45:71)]
################TRG ={0,1} (1,3+5) 
DB_ccB[DB_ccB[,17]<=1,17] <-0
DB_ccB[DB_ccB[,17]>=3,17] <-1
DBTEST_ccB[DBTEST_ccB[,17]<=1,17] <-0
DBTEST_ccB[DBTEST_ccB[,17]>=3,17] <-1

#MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_ccB<-DB_ccB[,-c(16)]
#DBTEST_ccB<-DBTEST_ccB[,-c(16)]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., DB_ccB)[,-1]
y <- factor(DB_ccB$TRG) # 0->21 ,1->115 , Proportion = 0.84
x.test <- model.matrix(TRG ~., DBTEST_ccB)[,-1]
y.test <- factor(DBTEST_ccB$TRG)

# Set lambda via cross validation
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")#, standardize = FALSE)
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "response")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
predicted <- ifelse(probabilities > 0.84, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 


#############CASO multinomiale
DB_ccMUL <- DB_train[,-c(45:71)]
DBTEST_ccMUL <- DB_test[,-c(45:71)]
####TRG ={1,3,5}  
#type.multinomial() allows the usage of a grouped lasso penalty
#(q=2->is a grouped-lasso penalty on all the K coefficients for a particular variable, which makes them 
#all be zero or nonzero together) if type.multinomial = "grouped". The default is type.multinomial = "ungrouped" 
#(q=1 -> lasso penalty on each of the parameters, called partial Newton algorithm)

#MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_ccMUL<-DB_ccMUL[,-c(16)]
#DBTEST_ccMUL<-DBTEST_ccMUL[,-c(16)]

#MODELLO
# Divide predicotors and target
x <- model.matrix(TRG ~., DB_ccMUL)[,-1]
y <- factor(DB_ccMUL$TRG) 
x.test <- model.matrix(TRG ~., DBTEST_ccMUL)[,-1]
y.test <- factor(DBTEST_ccMUL$TRG)

# Set lambda via cross validation
cv.lasso <- cv.glmnet(x, y, family = "multinomial",type.measure = "mse")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
#probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
#predicted <- ifelse(probabilities > 0.86, "1","0")
predicted <-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")

#Model Perfomance
mean(predicted == y.test)

#################################################################
#################CASO Dicotomico A (CLINCHE + CORE + MARGIN)
DB_ccmA <- DB_train
DBTEST_ccmA <- DB_test
################TRG ={0,1} (1+3,5)
DB_ccmA[DB_ccmA[,17]<=3,17] <-0
DB_ccmA[DB_ccmA[,17]>=5,17] <-1
DBTEST_ccmA[DBTEST_ccmA[,17]<=3,17] <-0
DBTEST_ccmA[DBTEST_ccmA[,17]>=5,17] <-1

#standardize Radiomics AND età,diam_mm,cea
#DB_ccmA[,c(2,7,8,18:44)] <- scale(DB_ccmA[,c(2,7,8,18:44)])

#MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_ccmA<-DB_ccmA[,-c(16)]
#DBTEST_ccmA<-DBTEST_ccmA[,-c(16)]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., DB_ccmA)[,-1]
y <- factor(DB_ccmA$TRG) # 0->42, 1->94 , Proportion = 0.69
x.test <- model.matrix(TRG ~., DBTEST_ccmA)[,-1]
y.test <- factor(DBTEST_ccmA$TRG)

# Set lambda via cross validation
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")#,standardize = FALSE)
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "response")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
predicted <- ifelse(probabilities > 0.69, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 

################CASO Dicotomico B 
DB_ccmB <- DB_train
DBTEST_ccmB <- DB_test
################TRG ={0,1} (1,3+5) 
DB_ccmB[DB_ccmB[,17]<=1,17] <-0
DB_ccmB[DB_ccmB[,17]>=3,17] <-1
DBTEST_ccmB[DBTEST_ccmB[,17]<=1,17] <-0
DBTEST_ccmB[DBTEST_ccmB[,17]>=3,17] <-1


#MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_ccmB <- DB_ccmB[,-c(16)]
#DBTEST_ccmB <- DBTEST_ccmB[,-c(16)]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., DB_ccmB)[,-1]
y <- factor(DB_ccmB$TRG) # # 0->21 ,1->115 , Proportion = 0.84
x.test <- model.matrix(TRG ~., DBTEST_ccmB)[,-1]
y.test <- factor(DBTEST_ccmB$TRG)

# Set lambda via cross validation
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "response")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
predicted <- ifelse(probabilities > 0.84, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 


#############CASO multinomiale 
DB_ccmMUL <- DB_train
DBTEST_ccmMUL <- DB_test
####TRG ={1,3,5}   1->24, 3->28, 5->117  
#type.multinomial() allows the usage of a grouped lasso penalty
#(q=2->is a grouped-lasso penalty on all the K coefficients for a particular variable, which makes them 
#all be zero or nonzero together) if type.multinomial = "grouped". The default is type.multinomial = "ungrouped" 
#(q=1 -> lasso penalty on each of the parameters, called partial Newton algorithm)


#MODELLO SENZA RIPOSTA RADIOLOGICA
#DB_ccmMUL<-DB_ccmMUL[,-c(16)]
#DBTEST_ccmMUL<-DBTEST_ccmMUL[,-c(16)]

#MODELLO
# Divide predicotors and target
x <- model.matrix(TRG ~., DB_ccmMUL)[,-1]
y <- factor(DB_ccmMUL$TRG) 
x.test <- model.matrix(TRG ~., DBTEST_ccmMUL)[,-1]
y.test <- factor(DBTEST_ccmMUL$TRG)

# Set lambda via cross validation
#type measure:deviance(default), class, mae,mse
#type.multinomial = grouped -> results are worse or similar
cv.lasso <- cv.glmnet(x, y, family = "multinomial",type.measure = "mse")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
#probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
#predicted <- ifelse(probabilities > ?, "1","0")
predicted <-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")

#Model Perfomance
mean(predicted == y.test)


###########################################################
############RANDOM FOREST
############caso multinomiale

# da importare i tre dataset DB_cliniche
# DB_cliniche_core, DB_cliniche_core_margin

library(randomForest)

set.seed(123)
training.samples <- DB_cliniche$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB_cliniche[training.samples, ]
test.data <- DB_cliniche[-training.samples, ]

x1 <- model.matrix(TRG ~., train.data)[,-1]
y1 <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)


modello_rf<-randomForest(x= x1,y =y1, ntree = 200, proximity =TRUE)
modello_rf

predicted2 <- predict(modello_rf, x.test)

confusionMatrix(as.factor(predicted2), as.factor(y.test))

##### cliniche + core
set.seed(123)
training.samples <- DB_cliniche_core$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB_cliniche_core[training.samples, ]
test.data <- DB_cliniche_core[-training.samples, ]

x1 <- model.matrix(TRG ~., train.data)[,-1]
y1 <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

modello_rf<-randomForest(x= x1,y =y1, ntree = 200, proximity =TRUE)
modello_rf

predicted2 <- predict(modello_rf, x.test)

confusionMatrix(as.factor(predicted2), as.factor(y.test))


###### cliniche + core + margin

set.seed(123)
training.samples <- DB_cliniche_core_margin$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB_cliniche_core_margin[training.samples, ]
test.data <- DB_cliniche_core_margin[-training.samples, ]

x1 <- model.matrix(TRG ~., train.data)[,-1]
y1 <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)


modello_rf<-randomForest(x= x1,y =y1, ntree = 200, proximity =TRUE)
modello_rf

predicted2 <- predict(modello_rf, x.test)

confusionMatrix(as.factor(predicted2), as.factor(y.test))
