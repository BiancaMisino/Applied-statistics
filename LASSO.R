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

#collinearità radiomiche margin
Database_finale<-Database_finale[,-c(81,82,83,85,87,94,95,97,98,101,103,104,107,108,109,110,116,118,120,121,123)]
#collinearità radiomiche core
Database_finale<-Database_finale[,-c(33,34,35,37,39,40,46,49,51,55,58,59,60,61,62,66,69,72,73,74,75)]

##### MODELLO _ milano+torino 
DB_cliniche <- Database_finale[,-c(29:82)]
DB_cl_core <- Database_finale[,-c(56:82)]
DB_cc_mar <- Database_finale

#levo le cliniche ridondanti
DB<-DB_cliniche[,-c(1,2,3,4,5,11,24,25)]

# per analisi successive
DB<-DB_cl_core[,-c(1,2,3,4,5,11,24,25)]
DB <- DB_cc_mar[,-c(1,2,3,4,5,11,24,25)]

#rendo numeriche e le variabili categoriche:
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

#M_CONVENTIONAL_HU_MIN (radiomica_margin)
DB[,48]<-as.numeric(unlist(DB[,48]))

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

#DBpul <- DB[-45,] (outlier CEA)

#elimino le variabili con NA
DB<-DB[,-c(10,11,12)]

####################################################################################
################CASO Dicotomico A (CLINICHE)
################TRG ={0,1} (1+3,5) 

DB[DB[,17]<=3,17] <-0
DB[DB[,17]>=5,17] <-1

Database_validazione[Database_validazione_[17]<=3,17] <-0
Database_validazione[Database_validazione[,17]>=5,17] <-1

#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) # 0->43, 1->93 , Proportion = 0.684
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

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
predicted <- ifelse(probabilities > 0.68, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 

#ROC curve
roc_pred <- prediction(predictions = probabilities , labels = y.test)
#roc_perf <- performance(roc_pred , "tpr" , "fpr")
#plot(roc_perf,colorize = TRUE,
#     print.cutoffs.at= seq(0,1,0.05),
#     text.adj=c(-0.2,1.7))
auc_ROCR <- performance(roc_pred, measure = "auc")
auc_ROCR@y.values[[1]]
                       


################CASO Dicotomico B 
################TRG ={0,1} (1,3+5) ->  0->24;  1->145
DB[DB[,17]<=1,17] <-0
DB[DB[,17]>=3,17] <-1

Database_validazione[Database_validazione_[,28]>=3,28] <-1
Database_validazione[Database_validazione[,28]==1,28] <-0

#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) # 0->19, 1->117 , Proportion = 0.8602941
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

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
predicted <- ifelse(probabilities > 0.86, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test)


####CASO multinomiale  
####TRG ={1,3,5}   1->24, 3->28, 5->117  
#type.multinomial() allows the usage of a grouped lasso penalty
#(q=2->is a grouped-lasso penalty on all the K coefficients for a particular variable, which makes them 
#all be zero or nonzero together) if type.multinomial = "grouped". The default is type.multinomial = "ungrouped" 
#(q=1 -> lasso penalty on each of the parameters, called partial Newton algorithm)

#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predicotors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

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
#predicted <- ifelse(probabilities > 0.86, "1","0")
predicted <-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")

#Model Perfomance
mean(predicted == y.test)

#################################################################
#################CASO Dicotomico A (CLINCHE + RADIOMICHE_CORE)
################TRG ={0,1} (1+3,5)
DB[DB[,17]<=3,17] <-0
DB[DB[,17]>=5,17] <-1

Database_validazione_lesione_gr[Database_validazione_lesione_gr[,28]<=3,28] <-0
Database_validazione_lesione_gr[Database_validazione_lesione_gr[,28]>=5,28] <-1
#DBpul <- DB[-45,]

#standardize Radiomics AND età,diam_mm,cea
#DB[,c(2,7,8,18:44)] <- scale(DB[,c(2,7,8,18:44)])

#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) # 0->43, 1->93 , Proportion = 0.684
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

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
predicted <- ifelse(probabilities > 0.68, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 

################CASO Dicotomico B 
################TRG ={0,1} (1,3+5) ->  0->24;  1->145
DB[DB[,17]<=1,17] <-0
DB[DB[,17]>=3,17] <-1

Database_validazione[Database_validazione_[,28]>=3,28] <-1
Database_validazione[Database_validazione[,28]==1,28] <-0


#standardize Radiomics AND età,diam_mm,cea
#DB[,c(2,7,8,18:44)] <- scale(DB[,c(2,7,8,18:44)])

#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) # 0->19 ,1->117 ,  Proportion = 0.8602941
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

# Set lambda via cross validation, standardize only the radiomics (see above)
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")#, standardize = FALSE)
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "response")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
predicted <- ifelse(probabilities > 0.86, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 


####CASO multinomiale  
####TRG ={1,3,5}   1->24, 3->28, 5->117  
#type.multinomial() allows the usage of a grouped lasso penalty
#(q=2->is a grouped-lasso penalty on all the K coefficients for a particular variable, which makes them 
#all be zero or nonzero together) if type.multinomial = "grouped". The default is type.multinomial = "ungrouped" 
#(q=1 -> lasso penalty on each of the parameters, called partial Newton algorithm)

#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predicotors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

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
#predicted <- ifelse(probabilities > 0.86, "1","0")
predicted <-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")

#Model Perfomance
mean(predicted == y.test)

#################################################################
#################CASO Dicotomico A (CLINCHE + CORE + MARGIN)
################TRG ={0,1} (1+3,5)
DB[DB[,17]<=3,17] <-0
DB[DB[,17]>=5,17] <-1

Database_validazione_lesione_gr[Database_validazione_lesione_gr[,28]<=3,28] <-0
Database_validazione_lesione_gr[Database_validazione_lesione_gr[,28]>=5,28] <-1
#DBpul <- DB[-45,]

#standardize Radiomics AND età,diam_mm,cea
#DB[,c(2,7,8,18:44)] <- scale(DB[,c(2,7,8,18:44)])

#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) # 0->43, 1->93 , Proportion = 0.684
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

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
predicted <- ifelse(probabilities > 0.68, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 

################CASO Dicotomico B 
################TRG ={0,1} (1,3+5) ->  0->24;  1->145
DB[DB[,17]<=1,17] <-0
DB[DB[,17]>=3,17] <-1

Database_validazione[Database_validazione_[,28]>=3,28] <-1
Database_validazione[Database_validazione[,28]==1,28] <-0

#standardize Radiomics AND età,diam_mm,cea
#DB[,c(2,7,8,18:44)] <- scale(DB[,c(2,7,8,18:44)])

#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predictors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) # 0->19 ,1->117 ,  Proportion = 0.8602941
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

# Set lambda via cross validation, standardize only the radiomics (see above)
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")#, standardize = FALSE)
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min
#cv.lasso$cvm #measure

# Display coefficients AND Make predictions on the test data
probabilities <- predict(cv.lasso, newx = x.test, s = "lambda.min", type = "response")
coeff<-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "coefficients")
coeff
predicted <- ifelse(probabilities > 0.86, "1","0")

#Model Perfomance
confusionMatrix(as.factor(predicted), as.factor(y.test)) #mode = "prec_recall"
AUC(predicted,y.test) 


####CASO multinomiale  
####TRG ={1,3,5}   1->24, 3->28, 5->117  
#type.multinomial() allows the usage of a grouped lasso penalty
#(q=2->is a grouped-lasso penalty on all the K coefficients for a particular variable, which makes them 
#all be zero or nonzero together) if type.multinomial = "grouped". The default is type.multinomial = "ungrouped" 
#(q=1 -> lasso penalty on each of the parameters, called partial Newton algorithm)


#Cliniche senza risposta_radiologica
DB<-DB[,-c(16)]

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predicotors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- factor(train.data$TRG) 
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- factor(test.data$TRG)

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
#predicted <- ifelse(probabilities > 0.86, "1","0")
predicted <-predict(cv.lasso, newx = x.test, s = "lambda.min", type = "class")

#Model Perfomance
mean(predicted == y.test)
