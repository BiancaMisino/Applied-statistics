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

#table(DB$TRG)
#1   3   5 
#24  28 117 

##### MODELLO CLINICHE _ milano+torino 
DB_cliniche <- Database_finale[,-c(29:124)]


#levo le cliniche ridondanti

DB<-DB_cliniche[,-c(1,2,3,4,5,11,24,25)]

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

#M_HU_CONVENTIONAL (radiomica)
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

#collinearità radiomiche
DB_coll<-DB_rad[,-c(6,7,8,10,12,13,19,22,24,28,31,32,33,34,35,39,42,45,46,47,48)]
#30 -> 0.86
#15->0.85
#40 -> 0.88
#DBpul <- DB[-45,] (outlier CEA)

#elimino le variabili con NA
DB<-DB[,-c(10,11,12)]

################CASO Dicotomico A 
################TRG ={0,1} (5,1+3)

DB[DB[,17]<=3,20] <-0
DB[DB[,17]>=5,20] <-1

Database_validazione[Database_validazione_[,28]<=3,28] <-0
Database_validazione[Database_validazione[,28]>=5,28] <-1


#SELECTION
x <- model.matrix(TRG ~., DB)[,-1]
y <- DB$TRG


#x <- as.matrix(DB[, -c(which(colnames(DB)=='TRG'))])
#standardize/normalize the variables
#x <- data.Normalization (x,type="n4",normalization="column") #n1 -> standardization #n4 ->(0,1)
#x <-scale(x) -> non serve perchè glmnet() lo fa di default 

# default alpha = 1 -> LASSO, alpha = 0 -> RIDGE, alpha (0,1) -> ELASTIC NET
cv.lasso <- cv.glmnet(x, y, family = "binomial", type.measure = "class")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min

# Display coefficients
model <- glmnet(x, y, family = "binomial",type.measure = "class",lambda = cv.lasso$lambda.min)
coef(model)

#parametri(lambda)
fit <- glmnet(x,y)
x11()
plot(fit,xvar='lambda',label=TRUE, col = rainbow(dim(x)[2]))
legend('topright', dimnames(x)[[2]], col =  rainbow(dim(x)[2]), lty=1, cex=1)

# Split the data into training and test set
set.seed(123)
training.samples <- DB$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB[training.samples, ]
test.data <- DB[-training.samples, ]

#MODELLO
# Divide predicotors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- train.data$TRG
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- test.data$TRG

# set lambda via cross validation

cv.lasso <- cv.glmnet(x, y, family = "binomial", type.measure = "auc")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min

# Display coefficients

model <- glmnet(x, y, family = "binomial",type.measure = "auc",lambda = cv.lasso$lambda.min)
coef(model)

# Make predictions on the test data
#probabilities <- predict(model,newx= x.test, type = "class")
#predicted.classes <- ifelse(probabilities > 0.5, "1", "0")

predicted.classes <- predict(model,newx= x.test, type = "class")

# Model accuracy
mean(predicted.classes == y.test)

 
################CASO Dicotomico A 
################TRG ={0,1} (5,1+3)
DB[DB[,28]<=3,28] <-0
DB[DB[,28]>=5,28] <-1

Database_validazione_lesione_gr[Database_validazione_lesione_gr[,28]<=3,28] <-0
Database_validazione_lesione_gr[Database_validazione_lesione_gr[,28]>=5,28] <-1
#DBpul <- DB[-45,]

DB_rad<-subset(DB[,c(28:76)])
DB_rad_val<-subset(DB[,c(28:76)])

#collinearità 
DB_coll<-DB_rad[,-c(6,7,8,10,12,13,19,22,24,28,31,32,33,34,35,39,42,45,46,47,48)]
#30 -> 0.86
#15->0.85
#40 -> 0.88

#SELECTION
x <- as.matrix(DB_rad[, -c(which(colnames(DB_rad)=='TRG'))])
y <- DB_rad$TRG

#standardize/normalize the variables
#x <- data.Normalization (x,type="n4",normalization="column") #n1 -> standardization #n4 ->(0,1)
x <-scale(x)

# default alpha = 1 -> LASSO, alpha = 0 -> RIDGE, alpha (0,1) -> ELASTIC NET
cv.lasso <- cv.glmnet(x, y, family = "binomial", type.measure = "auc")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min

# Display coefficients
model <- glmnet(x, y, family = "binomial",type.measure = "auc",lambda = cv.lasso$lambda.min)
coef(model)

#parametri(lambda)
fit <- glmnet(x,y)
x11()
plot(fit,xvar='lambda',label=TRUE, col = rainbow(dim(x)[2]))
legend('topright', dimnames(x)[[2]], col =  rainbow(dim(x)[2]), lty=1, cex=1)

# Split the data into training and test set
set.seed(123)
training.samples <- DB_rad$TRG %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- DB_rad[training.samples, ]
test.data <- DB_rad[-training.samples, ]

#MODELLO
# Divide predicotors and target
x <- model.matrix(TRG ~., train.data)[,-1]
y <- train.data$TRG
x.test <- model.matrix(TRG ~., test.data)[,-1]
y.test <- test.data$TRG

#standardization
x <-scale(x)
x.test <-scale(x.test)

# set lambda via cross validation

cv.lasso <- cv.glmnet(x, y, family = "binomial", type.measure = "auc")
plot(cv.lasso)
abline(v=log(cv.lasso$lambda.min), lty=1)
cv.lasso$lambda.min

# Display coefficients

model <- glmnet(x, y, family = "binomial",type.measure = "auc",lambda = cv.lasso$lambda.min)
coef(model)

# Make predictions on the test data
#probabilities <- predict(model,newx= x.test, type = "class")
#predicted.classes <- ifelse(probabilities > 0.5, "1", "0")

predicted.classes <- predict(model,newx= x.test, type = "class")

# Model accuracy
mean(predicted.classes == y.test)

################CASO Dicotomico B 
################TRG ={0,1} (5+3,1)


####CASO multinomiale  
####TRG ={1,3,5}