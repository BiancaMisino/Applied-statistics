setwd("~/martina/università/magistrale/I anno/II sem/applied statistics-secchi/Project/dati_09_05")

##Togliamo ID, codice_TC,...3,Interval_CT_Surgery,Numero_metastasi, Cicli_di_chemio e cicli_last_line
#Database_finale<-Database_finale1

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
DB_cliniche<-Database_finale[,c(1:3,5,7:20,117:123)]


##rimuovo gli NA che non sappiamo trattare
DB_cliniche<-DB_cliniche[,-c(8,9,10)]

DB_cliniche[DB_cliniche[,15]<=3,15] <-0
DB_cliniche[DB_cliniche[,15]>=5,15] <-1
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
#IL CORE CHECK
summary(DB_core)

#IL MARGIN
DB_margin<-Database_finale[,c(69:116)]

###tolgo le correlate del core
DB_margin<-DB_margin[,-c(5,6,7,9,11,18,19,21,22,25,27,28,31,32,33,34,40,42,44,45,47)]
#IL MARGIN CHECK
summary(DB_margin)

#####################################################################################################################
################################ SOLO VARIABILI CLINICHE ######################################################
################################### TRG 0 1###################################################################


####primo caso dicotomico TRG ={0,1} (1+3,5)
library(olsrr)
model <- lm(TRG ~ ., data = DB_cliniche)
k<-ols_step_forward_p(model,progress = 'TRUE')
plot(k)
k$model

#######ottengo lo stesso risultato
library(leaps)
library(ISLR)

regfit.step <- regsubsets(TRG~.,data=DB_cliniche,nvmax=22,method="forward")
summary(regfit.step)

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
which.max(reg.summary$adjr2)
coef(regfit.step,10)
max(summary(regfit.step)$adjr2)

###mi bastano 12 variabili per avere il miglior R^2 adj
# estimation on the full dataset
reg.best <- regsubsets(TRG~.,data=DB_cliniche, nvmax=20)
coef(reg.best,6)

#########################DA FARE MEGLIO#########################
#######################PROCEDO CON CROSS- VALIDATION###################################
### Choosing among models using the k-fold cross-validation approach
### (exhaustive search)

k <- 10

set.seed(1)
folds <- sample(1:k,nrow(DB_cliniche),replace=TRUE)
folds
table(folds)

# function that performs the prediction for regsubsets
predict.regsubsets <- function(object,newdata,id){
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form,newdata)
  coefi <- coef(object,id=id)
  xvars <- names(coefi)
  mat[,xvars]%*%coefi
}

cv.errors <- matrix(NA,k,19, dimnames=list(NULL, paste(1:19)))
for(j in 1:k){
  best.fit <- regsubsets(TRG~.,data=DB_cliniche[folds!=j,],nvmax=19)
  for(i in 1:19){
    pred <- predict(best.fit,DB_cliniche[folds==j,],id=i)
    predicted<-ifelse(pred>0.7,"1","0")  
    cv.errors[j,i] <- mean( (DB_cliniche$TRG[folds==j]-as.numeric(predicted))^2 )
  }
}
cv.errors
root.mean.cv.errors <- sqrt(apply(cv.errors,2,mean)) # average over the columns
root.mean.cv.errors

plot(root.mean.cv.errors,type='b')


which.min(root.mean.cv.errors)
points(which.min(root.mean.cv.errors),root.mean.cv.errors[which.min(root.mean.cv.errors)], col='red',pch=19)
# estimation on the full dataset
reg.best <- regsubsets(TRG~.,data=DB_cliniche, nvmax=19)
coef(reg.best,10)


################################################provo a togliere la risposta radiologica###################################
DB_cliniche_senzarisprad<-DB_cliniche[,-c(14)]

regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_senzarisprad,nvmax=22,method="forward")
summary(regfit.step)

x11(height=7,width=14)
par(mfrow=c(1,3))
plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

reg.summary <- summary(regfit.step_senza_risp_rad)
which.max(reg.summary$adjr2)
coef(regfit.step_senza_risp_rad,11)
max(summary(regfit.step_senza_risp_rad)$adjr2)


#####################################################################################################################
################################ SOLO VARIABILI CLINICHE+CORE ######################################################
################################### TRG 0 1###################################################################

DB_cliniche_core<-data.frame(DB_cliniche,DB_core)

regfit.step <- regsubsets(TRG~.,data=DB_cliniche_core,nvmax=30,method="forward")
summary(regfit.step)

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
which.max(reg.summary$adjr2)
coef(regfit.step,14)
max(summary(regfit.step)$adjr2)

###mi bastano 12 variabili per avere il miglior R^2 adj
# estimation on the full dataset
reg.best <- regsubsets(TRG~.,data=DB_cliniche_core, nvmax=30)
coef(reg.best,6)

#########################DA FARE MEGLIO#########################
#######################PROCEDO CON CROSS- VALIDATION###################################
### Choosing among models using the k-fold cross-validation approach
### (exhaustive search)

k <- 10

set.seed(1)
folds <- sample(1:k,nrow(DB_cliniche_core),replace=TRUE)
folds
table(folds)

# function that performs the prediction for regsubsets
predict.regsubsets <- function(object,newdata,id){
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form,newdata)
  coefi <- coef(object,id=id)
  xvars <- names(coefi)
  mat[,xvars]%*%coefi
}

cv.errors <- matrix(NA,k,46, dimnames=list(NULL, paste(1:46)))
for(j in 1:k){
  best.fit <- regsubsets(TRG~.,data=DB_cliniche_core[folds!=j,],nvmax=46)
  for(i in 1:46){
    pred <- predict(best.fit,DB_cliniche_core[folds==j,],id=i)
    predicted<-ifelse(pred>0.7,"1","0")  
    cv.errors[j,i] <- mean( (DB_cliniche_core$TRG[folds==j]-as.numeric(predicted))^2 )
  }
}
cv.errors
root.mean.cv.errors <- sqrt(apply(cv.errors,2,mean)) # average over the columns
root.mean.cv.errors

plot(root.mean.cv.errors,type='b')


which.min(root.mean.cv.errors)
points(which.min(root.mean.cv.errors),root.mean.cv.errors[which.min(root.mean.cv.errors)], col='red',pch=19)
# estimation on the full dataset
reg.best <- regsubsets(TRG~.,data=DB_cliniche, nvmax=19)
coef(reg.best,10)

################################################provo a togliere la risposta radiologica###################################
DB_cliniche_core_senzarisprad<-DB_cliniche_core[,-c(14)]

regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_core_senzarisprad,nvmax=30,method="forward")
summary(regfit.step)

x11(height=7,width=14)
par(mfrow=c(1,3))
plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

reg.summary <- summary(regfit.step_senza_risp_rad)
which.max(reg.summary$adjr2)
coef(regfit.step_senza_risp_rad,11)
max(summary(regfit.step_senza_risp_rad)$adjr2)


#####################################################################################################################
################################ SOLO VARIABILI CLINICHE+CORE+MARGIN ################################################
################################### TRG 0 1#########################################################################



DB_cliniche_core_margin<-data.frame(DB_cliniche,DB_core,DB_margin)

regfit.step <- regsubsets(TRG~.,data=DB_cliniche_core_margin,nvmax=50,really.big=TRUE,method="forward")
summary(regfit.step)

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
which.max(reg.summary$adjr2)
coef(regfit.step,19)
max(summary(regfit.step)$adjr2)

###mi bastano 12 variabili per avere il miglior R^2 adj
# estimation on the full dataset
reg.best <- regsubsets(TRG~.,data=DB_cliniche_core_margin, nvmax=50)
coef(reg.best,6)
################################################provo a togliere la risposta radiologica###################################
DB_cliniche_core_margin_senzarisprad<-DB_cliniche_core_margin[,-c(14)]

regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_core_margin_senzarisprad,nvmax=50,method="forward")
summary(regfit.step)

x11(height=7,width=14)
par(mfrow=c(1,3))
plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

reg.summary <- summary(regfit.step_senza_risp_rad)
which.max(reg.summary$adjr2)
coef(regfit.step_senza_risp_rad,11)
max(summary(regfit.step_senza_risp_rad)$adjr2)
