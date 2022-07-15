library(stats)
library(olsrr)
library(leaps)
library(ISLR)
library(caret)
library(cvAUC)
library( rms )
library(arm)
library(ResourceSelection)
library(pROC)
library(PRROC)
library(mice)
library(VIM)
library(lattice)
library(ggplot2)
##Togliamo ID, codice_TC,...3,Interval_CT_Surgery,Numero_metastasi, Cicli_di_chemio e cicli_last_line
Database_finale<-Database_finale1

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
colonne_con_na<-DB_cliniche[,c(8)]
DB_cliniche<-DB_cliniche[,-c(8,9,10)]

###tolgo quelle lin dipendenti-> nelle dummy le avevo aggiunte ma non aveva senso -> le ritolgo
DB_cliniche<-DB_cliniche[,-c(18,22)]



DB_cliniche<-data.frame(DB_cliniche,colonne_con_na)
#LE CLINICHE SONO PRONTE
summary(DB_cliniche)

#marginplot(DB_cliniche[c(21,22)])
#Here is an explanation of the parameters used:

#m  - Refers to 5 imputed data sets
#maxit - Refers to no. of iterations taken to impute missing values
#method - Refers to method used in imputation. we used predictive mean matching.

mice_plot <- aggr(colonne_con_na, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(colonne_con_na), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))

imputed_Data <- mice(DB_cliniche, m=5, maxit = 50, method = 'pmm', seed = 500)
summary(imputed_Data) 
##ho fatto 5 inputazioni e ne ho scelta una delle 5 per proseguire con le analisi
DB_cliniche <- complete(imputed_Data,2)

DB_cliniche[,2]<-scale(DB_cliniche[,2],center=FALSE,scale=TRUE)
DB_cliniche[,5]<-scale(DB_cliniche[,5],center=FALSE,scale=TRUE)
DB_cliniche[,6]<-scale(DB_cliniche[,6],center=FALSE,scale=TRUE)

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

#####################################################################################################################
################################ SOLO VARIABILI CLINICHE ######################################################
################################### TRG 0 1###################################################################

DB_cliniche_casoa=DB_cliniche
DB_cliniche_casoa[DB_cliniche_casoa[,15]<=3,15] <-0
DB_cliniche_casoa[DB_cliniche_casoa[,15]>3,15] <-1

#variables selection
reg0<-glm(TRG~1,data=DB_cliniche_casoa,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_casoa,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients
summary(best.fit)
##interval greater than 30-> lo levo


#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casoa[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casoa[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~ Risposta_radiologica + ANTI_VEGF + CEA
             + Linee_di_chemioterapia+OXALIPLATINO  ,
             data=DB_cliniche_casoa[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_casoa[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~ Risposta_radiologica + ANTI_VEGF + CEA
           + Linee_di_chemioterapia+OXALIPLATINO  ,
           data=DB_cliniche_casoa, family=binomial)

odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI

fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_casoa$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

################################################provo a togliere la risposta radiologica###################################
DB_cliniche_senzarispradA<-DB_cliniche_casoa[,-c(14)]

#regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_senzarisprad,nvmax=22,method="forward")
#summary(regfit.step)

#x11(height=7,width=14)
#par(mfrow=c(1,3))
#plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
#plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
#plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

#reg.summary <- summary(regfit.step_senza_risp_rad)
#which.max(reg.summary$adjr2)
#coef(regfit.step_senza_risp_rad,11)
#max(summary(regfit.step_senza_risp_rad)$adjr2)
reg0<-glm(TRG~1,data=DB_cliniche_senzarispradA,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_senzarispradA,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients
summary(best.fit)

##

#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_senzarispradA[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_senzarispradA[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~ CEA + colonne_con_na + sync + ANTI_VEGF + 
               primary_tumor_site_left  ,
             data=DB_cliniche_senzarispradA[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_senzarispradA[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_senzarispradA[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_senzarispradA[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~ CEA + colonne_con_na + sync + ANTI_VEGF + 
             primary_tumor_site_left  ,
           data=DB_cliniche_senzarispradA, family=binomial)


odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))
CI<-confint(model,level=0.95)
odds.ratio
CI

fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_senzarispradA$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)
#######################################################CASO B#####################################################

####secondo caso dicotomico TRG ={0,1} (1,3+5)
DB_cliniche_casob=DB_cliniche
DB_cliniche_casob[DB_cliniche_casob[,15]<3,15] <-0
DB_cliniche_casob[DB_cliniche_casob[,15]>=3,15] <-1


#variables selection
reg0<-glm(TRG~1,data=DB_cliniche_casob,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_casob,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients

summary(best.fit)
### RISP RAD HA PVALUE ALTISSIMO , LINEE DI CHEMIO
## LINEE DI CHEMIO TOLTE

#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casob[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casob[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~ Risposta_radiologica + primary_tumor_site_left + 
               ANTI_VEGF + Linee_di_chemioterapia   
               ,
             data=DB_cliniche_casob[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casob[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.86, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casob[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_casob[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~  Risposta_radiologica + primary_tumor_site_left + 
             ANTI_VEGF + Linee_di_chemioterapia
             ,
           data=DB_cliniche_casob, family=binomial)

odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))
CI<-confint(model,level=0.95)
odds.ratio
CI


################################################provo a togliere la risposta radiologica###################################
DB_cliniche_senzarispradB<-DB_cliniche_casob[,-c(14)]

#regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_senzarisprad,nvmax=22,method="forward")
#summary(regfit.step)

#x11(height=7,width=14)
#par(mfrow=c(1,3))
#plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
#plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
#plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

#reg.summary <- summary(regfit.step_senza_risp_rad)
#which.max(reg.summary$adjr2)
#coef(regfit.step_senza_risp_rad,11)
#max(summary(regfit.step_senza_risp_rad)$adjr2)
reg0<-glm(TRG~1,data=DB_cliniche_senzarispradB,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_senzarispradB,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients

summary(best.fit)
##va tolto solo diam in mm 

#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_senzarispradB[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_senzarispradB[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~  primary_tumor_site_left + Diametro_in_mm + 
               sync + ANTI_VEGF + Linee_di_chemioterapia + bilater + greater_6_cicli  ,
             data=DB_cliniche_senzarispradB[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_senzarispradB[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.86, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_senzarispradB[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_senzarispradB[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~ primary_tumor_site_left + Diametro_in_mm + 
             sync + ANTI_VEGF + Linee_di_chemioterapia + bilater + greater_6_cicli  ,
           data=DB_cliniche_senzarispradB, family=binomial)

odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))
CI<-confint(model,level=0.95)
odds.ratio
CI


#####################################################################################################################
################################ SOLO VARIABILI CLINICHE+CORE ######################################################
################################### TRG 0 1###################################################################

DB_cliniche_core<-data.frame(DB_cliniche,DB_core)

DB_cliniche_casoa=DB_cliniche_core
DB_cliniche_casoa[DB_cliniche_casoa[,15]<=3,15] <-0
DB_cliniche_casoa[DB_cliniche_casoa[,15]>3,15] <-1

#variables selection
reg0<-glm(TRG~1,data=DB_cliniche_casoa,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_casoa,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients

summary(best.fit)
##TOLGO OXALIPLATINO
#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casoa[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casoa[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA  
  model<-glm(TRG ~ Risposta_radiologica + ANTI_VEGF + CEA + 
               SHAPE_Volume..mL. + Diametro_in_mm + sync + GLRLM_GLNU + 
               primary_tumor_site_left + GLRLM_SRLGE   
                 ,
             data=DB_cliniche_casoa[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_casoa[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~ Risposta_radiologica + ANTI_VEGF + CEA + 
             SHAPE_Volume..mL. + Diametro_in_mm + sync + GLRLM_GLNU + 
             primary_tumor_site_left + GLRLM_SRLGE 
             ,data=DB_cliniche_casoa,family = binomial)
#n=9
#st.dev<-summary(model)$coefficients[2:10,2]
#beta_estimate<-summary(model)$coefficients[2:10,1]
#C<-diag(1,9)
#cfr.t<-qt(1-0.05/2,n-1)
#CI<-cbind(inf = C%*%beta_estimate-cfr.t*st.dev, center= C%*%beta_estimate, sup= C%*%beta_estimate+cfr.t*st.dev)

#colnames(CI)<-c("Inf","Mean","Sup")
#rownames(CI)<-c("Risposta_radiologica","ANTI_VEGF","CEA ","SHAPE_Volume..mL.","Diametro_in_mm ","sync ","GLRLM_GLNU ","primary_tumor_site_left","GLRLM_SRLGE")
#CI
odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI

fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_casoa$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

################################################provo a togliere la risposta radiologica###################################
DB_cliniche_senzarispradA<-DB_cliniche_casoa[,-c(14)]

#regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_senzarisprad,nvmax=22,method="forward")
#summary(regfit.step)

#x11(height=7,width=14)
#par(mfrow=c(1,3))
#plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
#plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
#plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

#reg.summary <- summary(regfit.step_senza_risp_rad)
#which.max(reg.summary$adjr2)
#coef(regfit.step_senza_risp_rad,11)
#max(summary(regfit.step_senza_risp_rad)$adjr2)
reg0<-glm(TRG~1,data=DB_cliniche_senzarispradA,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_senzarispradA,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients

summary(best.fit)
##TOLGO PRIMARY TUMOR SITE LEFT

#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_senzarispradA[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_senzarispradA[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~ CEA + colonne_con_na + GLRLM_SRLGE + sync + 
               Interval._greater_30.days + primary_tumor_site_left + NGLDM_Busyness + 
               Diametro_in_mm + ANTI_VEGF + HISTO_ExcessKurtosis + GLZLM_LZE   ,
             data=DB_cliniche_senzarispradA[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_senzarispradA[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_senzarispradA[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_senzarispradA[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~ CEA + colonne_con_na + GLRLM_SRLGE + sync + 
             Interval._greater_30.days + primary_tumor_site_left + NGLDM_Busyness + 
             Diametro_in_mm + ANTI_VEGF + HISTO_ExcessKurtosis + GLZLM_LZE  ,
           data=DB_cliniche_casoa, family=binomial)

odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI


fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_senzarispradA$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)
#######################################################CASO B#####################################################

####secondo caso dicotomico TRG ={0,1} (1,3+5)
DB_cliniche_casob=DB_cliniche_core
DB_cliniche_casob[DB_cliniche_casob[,15]<3,15] <-0
DB_cliniche_casob[DB_cliniche_casob[,15]>=3,15] <-1


#variables selection
reg0<-glm(TRG~1,data=DB_cliniche_casob,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_casob,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients

summary(best.fit)

##TOLGO RISP RAD , SYNC LINEE DI CHEMIO
#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casob[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casob[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~ Risposta_radiologica + GLZLM_LGZE + primary_tumor_site_left + 
               GLZLM_ZP + sync + Linee_di_chemioterapia  ,
             data=DB_cliniche_casob[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casob[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.86, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casob[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_casob[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~  Risposta_radiologica + GLZLM_LGZE + primary_tumor_site_left + 
             GLZLM_ZP + sync + Linee_di_chemioterapia ,
           data=DB_cliniche_casob, family=binomial)


odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI


################################################provo a togliere la risposta radiologica###################################
DB_cliniche_senzarispradB<-DB_cliniche_casob[,-c(14)]

#regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_senzarisprad,nvmax=22,method="forward")
#summary(regfit.step)

#x11(height=7,width=14)
#par(mfrow=c(1,3))
#plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
#plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
#plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

#reg.summary <- summary(regfit.step_senza_risp_rad)
#which.max(reg.summary$adjr2)
#coef(regfit.step_senza_risp_rad,11)
#max(summary(regfit.step_senza_risp_rad)$adjr2)
reg0<-glm(TRG~1,data=DB_cliniche_senzarispradB,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_senzarispradB,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients
vif(best.fit)
summary(best.fit)
##TOLGO IRINOTECAN

#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_senzarispradB[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_senzarispradB[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~ GLRLM_SRLGE + primary_tumor_site_left + Linee_di_chemioterapia + 
               greater_6_cicli + OXALIPLATINO + GLCM_Homogeneity...Inverse.difference. + 
               IRINOTECAN + Diametro_in_mm  ,
             data=DB_cliniche_senzarispradB[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_senzarispradB[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.86, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_senzarispradB[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_senzarispradB[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
###CAMBIA MODELLO OGNI VOLTA
best.fit$coefficients
model<-glm(TRG ~ GLRLM_SRLGE + primary_tumor_site_left + Linee_di_chemioterapia + 
             greater_6_cicli + OXALIPLATINO + GLCM_Homogeneity...Inverse.difference. + 
             IRINOTECAN + Diametro_in_mm   ,
           data=DB_cliniche_senzarispradB, family=binomial)

odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI

#####################################################################################################################
################################ SOLO VARIABILI CLINICHE+CORE+MARGIN ################################################
################################### TRG 0 1#########################################################################



DB_cliniche_core_margin<-data.frame(DB_cliniche,DB_core,DB_margin)

DB_cliniche_casoa=DB_cliniche_core_margin
DB_cliniche_casoa[DB_cliniche_casoa[,15]<=3,15] <-0
DB_cliniche_casoa[DB_cliniche_casoa[,15]>3,15] <-1

#variables selection
reg0<-glm(TRG~1,data=DB_cliniche_casoa,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_casoa,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients
vif(best.fit)
summary(best.fit)
##TOLGO OXALIPLATINO, 

#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casoa[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casoa[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~ Risposta_radiologica + ANTI_VEGF + CEA + 
               M_SHAPE_Sphericity..only.for.3D.ROI..nz.1. + M_HISTO_Entropy_log2 + 
               GLRLM_SRLGE + M_GLRLM_LGRE + M_CONVENTIONAL_HUmean + OXALIPLATINO + 
               M_GLZLM_SZE + M_GLZLM_GLNU + Diametro_in_mm + numerometastasi_2to3 + 
               SHAPE_Compacity.only.for.3D.ROI..nz.1.  ,
             data=DB_cliniche_casoa[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_casoa[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~ Risposta_radiologica + ANTI_VEGF + CEA + 
             M_SHAPE_Sphericity..only.for.3D.ROI..nz.1. + M_HISTO_Entropy_log2 + 
             GLRLM_SRLGE + M_GLRLM_LGRE + M_CONVENTIONAL_HUmean + OXALIPLATINO + 
             M_GLZLM_SZE + M_GLZLM_GLNU + Diametro_in_mm + numerometastasi_2to3 + 
             SHAPE_Compacity.only.for.3D.ROI..nz.1. ,
           data=DB_cliniche_casoa, family=binomial)


odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI

soglia = 0.68

valori.reali  = DB_cliniche_casoa$TRG
valori.predetti = as.numeric( model$fitted.values > soglia )
# 1 se > soglia, 0 se < = soglia
valori.predetti

tab = table( valori.reali, valori.predetti )

tab

# La tabella riportata è detta matrice di confusione, e riporta le osservazioni dette 
# Veri Positivi (True Positive o TP, osservazioni 1 classificate come 1), 
# Veri Negativi (True Negative o TN, osservazioni 0 classificate come 0), 
# Falsi Positivi (False Positive o FP, osservazioni 0 classificati come 1), 
# Falsi Negativi (Falsi Negativi o FN, osservazioni 1 classificati come 0). 

# Ci sono numerose metriche che permettono di valutare le performance del modello, a seconda delle esigenze:
# Accuracy, Sensitivity, Specificity

# Accuracy

# % di casi classificati correttamente:
round( sum( diag( tab ) ) / sum( tab ), 2 )

# % di casi misclassificati:
round( ( tab [ 1, 2 ] + tab [ 2, 1 ] ) / sum( tab ), 2 )

# Sensitivity
sensitivita =  tab [ 2, 2 ] /( tab [ 2, 1 ] + tab [ 2, 2 ] ) 
sensitivita

# Specificity 
specificita = tab[ 1, 1 ] /( tab [ 1, 2 ] + tab [ 1, 1 ] )
specificita



## 3. Curva ROC

# Costruire la Curva ROC a partire dai valori predetti per la risposta dal modello `mod.low2` 
# dell'analisi della variabile LOWBT.

# Le curve ROC (Receiver Operating Characteristic, anche note come Relative Operating Characteristic) 
# sono degli schemi grafici per un classificatore binario.

# Una curva ROC è il grafico dell'insieme delle coppie (FP, TP) al variare di un parametro del classificatore. 

fit2 = model$fitted


#media campionaria della prob di sopravvivenza nel campione

soglia_roc  = seq( 0, 1, length.out = 2e2 )
lens = length( soglia_roc )-1
ascissa_roc  = rep( NA, lens )
ordinata_roc = rep( NA, lens )

for ( k in 1 : lens )
{
  soglia = soglia_roc [ k ]
  
  classification = as.numeric( sapply( fit2, function( x ) ifelse( x < soglia, 0, 1 ) ) )
  
  #  ATTENZIONE, voglio sulle righe il vero e sulle colonne il predetto
  # t.misc = table( lw$LOW, classification )
  
  ordinata_roc[ k ] = sum( classification[ which( DB_cliniche_casoa$TRG == 1 ) ] == 1 ) /
    length( which( DB_cliniche_casoa$TRG == 1 ) )
  
  ascissa_roc[ k ] = sum( classification[ which( DB_cliniche_casoa$TRG == 0 ) ] == 1 ) /
    length( which( DB_cliniche_casoa$TRG == 0 ) )
  
  # ordinata_roc [ k ]  = t.misc [ 1, 1 ] /( t.misc [ 1, 1 ] + t.misc [ 1, 2 ] )
  #
  # ascissa_roc [ k ]  = t.misc [ 2, 1 ] /( t.misc [ 2, 1 ] + t.misc [ 2, 2 ] )
}


# Visualizziamo la curva ROC.

plot( ascissa_roc, ordinata_roc, type = "l", xlab = "1 - Specificity", ylab = "Sensitivity",
      main = "Curva ROC", lwd = 2, col = 'darkblue', ylim = c( 0, 1 ), xlim = c( 0, 1 ) )
abline( h = c( 0, 1 ), v = c( 0, 1 ), lwd = 1, lty = 2, col = 'red' )
abline( a = 0, b = 1, lty = 2, col = 'black' )

# qual era il nostro punto?
abline( v = 1 - specificita,  h = sensitivita, lty = 3, col = 'blue' )
points( 1 - specificita, sensitivita, pch = 4, lwd = 3, cex = 1.5, col = 'blue' )

# Le linee tratteggiate corrispondono alle due metriche calcolate con la threshold = 0.5 che abbiamo scelto. 

# Attraverso l'analisi delle curve ROC si valuta la capacità del classificatore calcolando 
# l'area sottesa alla curva ROC (Area Under Curve, AUC). 

# R fa tutto in automatico --> calcola AUC e p0ottimale 

PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_casoa$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)


################################################provo a togliere la risposta radiologica###################################
DB_cliniche_senzarispradA<-DB_cliniche_casoa[,-c(14)]

#regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_senzarisprad,nvmax=22,method="forward")
#summary(regfit.step)

#x11(height=7,width=14)
#par(mfrow=c(1,3))
#plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
#plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
#plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

#reg.summary <- summary(regfit.step_senza_risp_rad)
#which.max(reg.summary$adjr2)
#coef(regfit.step_senza_risp_rad,11)
#max(summary(regfit.step_senza_risp_rad)$adjr2)
reg0<-glm(TRG~1,data=DB_cliniche_senzarispradA,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_senzarispradA,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients
vif(best.fit)
summary(best.fit)

#model selection in 5-fold cv con variabili selezionate nel best fit 
##TOLGO MOLTE VARIABILI
k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_senzarispradA[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_senzarispradA[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~ CEA + colonne_con_na + M_GLRLM_LRLGE + GLRLM_SRLGE + 
               M_GLZLM_HGZE + OXALIPLATINO + GLCM_Contrast...Variance. + 
               M_GLZLM_GLNU + HISTO_Energy...Uniformity. + HISTO_Entropy_log2 + 
               Interval._greater_30.days + Diametro_in_mm + IRINOTECAN + 
               Linee_di_chemioterapia + GLZLM_ZP   ,
             data=DB_cliniche_senzarispradA[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_senzarispradA[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.68, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_senzarispradA[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_senzarispradA[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~ CEA + colonne_con_na + M_GLRLM_LRLGE + GLRLM_SRLGE + 
             M_GLZLM_HGZE + OXALIPLATINO + GLCM_Contrast...Variance. + 
             M_GLZLM_GLNU + HISTO_Energy...Uniformity. + HISTO_Entropy_log2 + 
             Interval._greater_30.days + Diametro_in_mm + IRINOTECAN + 
             Linee_di_chemioterapia + GLZLM_ZP ,
           data=DB_cliniche_senzarispradA, family=binomial)

odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI

soglia = 0.68

valori.reali  = DB_cliniche_senzarispradA$TRG
valori.predetti = as.numeric( model$fitted.values > soglia )
# 1 se > soglia, 0 se < = soglia
valori.predetti

tab = table( valori.reali, valori.predetti )

tab

# La tabella riportata è detta matrice di confusione, e riporta le osservazioni dette 
# Veri Positivi (True Positive o TP, osservazioni 1 classificate come 1), 
# Veri Negativi (True Negative o TN, osservazioni 0 classificate come 0), 
# Falsi Positivi (False Positive o FP, osservazioni 0 classificati come 1), 
# Falsi Negativi (Falsi Negativi o FN, osservazioni 1 classificati come 0). 

# Ci sono numerose metriche che permettono di valutare le performance del modello, a seconda delle esigenze:
# Accuracy, Sensitivity, Specificity

# Accuracy

# % di casi classificati correttamente:
round( sum( diag( tab ) ) / sum( tab ), 2 )

# % di casi misclassificati:
round( ( tab [ 1, 2 ] + tab [ 2, 1 ] ) / sum( tab ), 2 )

# Sensitivity
sensitivita =  tab [ 2, 2 ] /( tab [ 2, 1 ] + tab [ 2, 2 ] ) 
sensitivita

# Specificity 
specificita = tab[ 1, 1 ] /( tab [ 1, 2 ] + tab [ 1, 1 ] )
specificita



## 3. Curva ROC

# Costruire la Curva ROC a partire dai valori predetti per la risposta dal modello `mod.low2` 
# dell'analisi della variabile LOWBT.

# Le curve ROC (Receiver Operating Characteristic, anche note come Relative Operating Characteristic) 
# sono degli schemi grafici per un classificatore binario.

# Una curva ROC è il grafico dell'insieme delle coppie (FP, TP) al variare di un parametro del classificatore. 

fit2 = model$fitted


#media campionaria della prob di sopravvivenza nel campione

soglia_roc  = seq( 0, 1, length.out = 2e2 )
lens = length( soglia_roc )-1
ascissa_roc  = rep( NA, lens )
ordinata_roc = rep( NA, lens )

for ( k in 1 : lens )
{
  soglia = soglia_roc [ k ]
  
  classification = as.numeric( sapply( fit2, function( x ) ifelse( x < soglia, 0, 1 ) ) )
  
  #  ATTENZIONE, voglio sulle righe il vero e sulle colonne il predetto
  # t.misc = table( lw$LOW, classification )
  
  ordinata_roc[ k ] = sum( classification[ which( DB_cliniche_senzarispradA$TRG == 1 ) ] == 1 ) /
    length( which( DB_cliniche_senzarispradA$TRG == 1 ) )
  
  ascissa_roc[ k ] = sum( classification[ which( DB_cliniche_senzarispradA$TRG == 0 ) ] == 1 ) /
    length( which( DB_cliniche_senzarispradA$TRG == 0 ) )
  
  # ordinata_roc [ k ]  = t.misc [ 1, 1 ] /( t.misc [ 1, 1 ] + t.misc [ 1, 2 ] )
  #
  # ascissa_roc [ k ]  = t.misc [ 2, 1 ] /( t.misc [ 2, 1 ] + t.misc [ 2, 2 ] )
}


# Visualizziamo la curva ROC.

plot( ascissa_roc, ordinata_roc, type = "l", xlab = "1 - Specificity", ylab = "Sensitivity",
      main = "Curva ROC", lwd = 2, col = 'darkblue', ylim = c( 0, 1 ), xlim = c( 0, 1 ) )
abline( h = c( 0, 1 ), v = c( 0, 1 ), lwd = 1, lty = 2, col = 'red' )
abline( a = 0, b = 1, lty = 2, col = 'black' )

# qual era il nostro punto?
abline( v = 1 - specificita,  h = sensitivita, lty = 3, col = 'blue' )
points( 1 - specificita, sensitivita, pch = 4, lwd = 3, cex = 1.5, col = 'blue' )

# Le linee tratteggiate corrispondono alle due metriche calcolate con la threshold = 0.5 che abbiamo scelto. 

# Attraverso l'analisi delle curve ROC si valuta la capacità del classificatore calcolando 
# l'area sottesa alla curva ROC (Area Under Curve, AUC). 

# R fa tutto in automatico --> calcola AUC e p0ottimale 

PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_senzarispradA$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

#######################################################CASO B#####################################################

####secondo caso dicotomico TRG ={0,1} (1,3+5)
DB_cliniche_casob=DB_cliniche_core_margin
DB_cliniche_casob[DB_cliniche_casob[,15]<3,15] <-0
DB_cliniche_casob[DB_cliniche_casob[,15]>=3,15] <-1


#variables selection
reg0<-glm(TRG~1,data=DB_cliniche_casob,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_casob,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients
vif(best.fit)
summary(best.fit)
##TOLGO RISP RAD (ENORME),LINEE DI CHEMIO E SYNC

#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casob[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casob[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~  Risposta_radiologica + GLZLM_LGZE + primary_tumor_site_left + 
               M_CONVENTIONAL_HUstd + CONVENTIONAL_HUstd + Linee_di_chemioterapia + 
               M_NGLDM_Contrast
                 ,
             data=DB_cliniche_casob[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casob[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.86, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casob[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_casob[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~  Risposta_radiologica + GLZLM_LGZE + primary_tumor_site_left + 
             M_CONVENTIONAL_HUstd + CONVENTIONAL_HUstd + Linee_di_chemioterapia + 
             M_NGLDM_Contrast ,
           data=DB_cliniche_casob, family=binomial)

odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI

soglia = 0.86

valori.reali  = DB_cliniche_casob$TRG
valori.predetti = as.numeric( model$fitted.values > soglia )
# 1 se > soglia, 0 se < = soglia
valori.predetti

tab = table( valori.reali, valori.predetti )

tab

# La tabella riportata è detta matrice di confusione, e riporta le osservazioni dette 
# Veri Positivi (True Positive o TP, osservazioni 1 classificate come 1), 
# Veri Negativi (True Negative o TN, osservazioni 0 classificate come 0), 
# Falsi Positivi (False Positive o FP, osservazioni 0 classificati come 1), 
# Falsi Negativi (Falsi Negativi o FN, osservazioni 1 classificati come 0). 

# Ci sono numerose metriche che permettono di valutare le performance del modello, a seconda delle esigenze:
# Accuracy, Sensitivity, Specificity

# Accuracy

# % di casi classificati correttamente:
round( sum( diag( tab ) ) / sum( tab ), 2 )

# % di casi misclassificati:
round( ( tab [ 1, 2 ] + tab [ 2, 1 ] ) / sum( tab ), 2 )

# Sensitivity
sensitivita =  tab [ 2, 2 ] /( tab [ 2, 1 ] + tab [ 2, 2 ] ) 
sensitivita

# Specificity 
specificita = tab[ 1, 1 ] /( tab [ 1, 2 ] + tab [ 1, 1 ] )
specificita



## 3. Curva ROC

# Costruire la Curva ROC a partire dai valori predetti per la risposta dal modello `mod.low2` 
# dell'analisi della variabile LOWBT.

# Le curve ROC (Receiver Operating Characteristic, anche note come Relative Operating Characteristic) 
# sono degli schemi grafici per un classificatore binario.

# Una curva ROC è il grafico dell'insieme delle coppie (FP, TP) al variare di un parametro del classificatore. 

fit2 = model$fitted


#media campionaria della prob di sopravvivenza nel campione

soglia_roc  = seq( 0, 1, length.out = 2e2 )
lens = length( soglia_roc )-1
ascissa_roc  = rep( NA, lens )
ordinata_roc = rep( NA, lens )

for ( k in 1 : lens )
{
  soglia = soglia_roc [ k ]
  
  classification = as.numeric( sapply( fit2, function( x ) ifelse( x < soglia, 0, 1 ) ) )
  
  #  ATTENZIONE, voglio sulle righe il vero e sulle colonne il predetto
  # t.misc = table( lw$LOW, classification )
  
  ordinata_roc[ k ] = sum( classification[ which( DB_cliniche_casob$TRG == 1 ) ] == 1 ) /
    length( which( DB_cliniche_casob$TRG == 1 ) )
  
  ascissa_roc[ k ] = sum( classification[ which( DB_cliniche_casob$TRG == 0 ) ] == 1 ) /
    length( which( DB_cliniche_casob$TRG == 0 ) )
  
  # ordinata_roc [ k ]  = t.misc [ 1, 1 ] /( t.misc [ 1, 1 ] + t.misc [ 1, 2 ] )
  #
  # ascissa_roc [ k ]  = t.misc [ 2, 1 ] /( t.misc [ 2, 1 ] + t.misc [ 2, 2 ] )
}


# Visualizziamo la curva ROC.

plot( ascissa_roc, ordinata_roc, type = "l", xlab = "1 - Specificity", ylab = "Sensitivity",
      main = "Curva ROC", lwd = 2, col = 'darkblue', ylim = c( 0, 1 ), xlim = c( 0, 1 ) )
abline( h = c( 0, 1 ), v = c( 0, 1 ), lwd = 1, lty = 2, col = 'red' )
abline( a = 0, b = 1, lty = 2, col = 'black' )

# qual era il nostro punto?
abline( v = 1 - specificita,  h = sensitivita, lty = 3, col = 'blue' )
points( 1 - specificita, sensitivita, pch = 4, lwd = 3, cex = 1.5, col = 'blue' )

# Le linee tratteggiate corrispondono alle due metriche calcolate con la threshold = 0.5 che abbiamo scelto. 

# Attraverso l'analisi delle curve ROC si valuta la capacità del classificatore calcolando 
# l'area sottesa alla curva ROC (Area Under Curve, AUC). 

# R fa tutto in automatico --> calcola AUC e p0ottimale 

PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_casob$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)
################################################provo a togliere la risposta radiologica###################################
DB_cliniche_senzarispradB<-DB_cliniche_casob[,-c(14)]

#regfit.step_senza_risp_rad <- regsubsets(TRG~.,data=DB_cliniche_senzarisprad,nvmax=22,method="forward")
#summary(regfit.step)

#x11(height=7,width=14)
#par(mfrow=c(1,3))
#plot(summary(regfit.step_senza_risp_rad)$rsq,xlab="Number of Variables",ylab="R-squared",type="b")
#plot(summary(regfit.step_senza_risp_rad)$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="b")
#plot(summary(regfit.step_senza_risp_rad)$rss,xlab="Number of Variables",ylab="RSS",type="b")

#reg.summary <- summary(regfit.step_senza_risp_rad)
#which.max(reg.summary$adjr2)
#coef(regfit.step_senza_risp_rad,11)
#max(summary(regfit.step_senza_risp_rad)$adjr2)
reg0<-glm(TRG~1,data=DB_cliniche_senzarispradB,family = binomial)
reg1<-glm(TRG~.,data=DB_cliniche_senzarispradB,family = binomial)
best.fit<-step(reg0,scope=formula(reg1), direction="forward",k=2)
best.fit$coefficients
vif(best.fit)
summary(best.fit)
##TOLGO IRINOTECAN

#model selection in 5-fold cv con variabili selezionate nel best fit 

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_senzarispradB[1:127,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_senzarispradB[128:169,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,169)
folds[1:127]<-folds_mi
folds[128:169]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:5) {
  ###CAMBIA MODELLO OGNI VOLTA
  model<-glm(TRG ~  GLRLM_SRLGE + primary_tumor_site_left + Linee_di_chemioterapia + 
               greater_6_cicli + M_GLZLM_SZE + GLZLM_LGZE + sync + Diametro_in_mm  ,
             data=DB_cliniche_senzarispradB[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_senzarispradB[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.86, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_senzarispradB[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_senzarispradB[folds==j,]$TRG)
}
accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)
best.fit$coefficients
###CAMBIA MODELLO OGNI VOLTA
model<-glm(TRG ~  GLRLM_SRLGE + primary_tumor_site_left + Linee_di_chemioterapia + 
             greater_6_cicli + M_GLZLM_SZE + GLZLM_LGZE + sync + Diametro_in_mm    ,
           data=DB_cliniche_senzarispradB, family=binomial)

odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))

CI<-confint(model,level=0.95)
odds.ratio
CI

soglia = 0.86

valori.reali  = DB_cliniche_senzarispradB$TRG
valori.predetti = as.numeric( model$fitted.values > soglia )
# 1 se > soglia, 0 se < = soglia
valori.predetti

tab = table( valori.reali, valori.predetti )

tab

# La tabella riportata è detta matrice di confusione, e riporta le osservazioni dette 
# Veri Positivi (True Positive o TP, osservazioni 1 classificate come 1), 
# Veri Negativi (True Negative o TN, osservazioni 0 classificate come 0), 
# Falsi Positivi (False Positive o FP, osservazioni 0 classificati come 1), 
# Falsi Negativi (Falsi Negativi o FN, osservazioni 1 classificati come 0). 

# Ci sono numerose metriche che permettono di valutare le performance del modello, a seconda delle esigenze:
# Accuracy, Sensitivity, Specificity

# Accuracy

# % di casi classificati correttamente:
round( sum( diag( tab ) ) / sum( tab ), 2 )

# % di casi misclassificati:
round( ( tab [ 1, 2 ] + tab [ 2, 1 ] ) / sum( tab ), 2 )

# Sensitivity
sensitivita =  tab [ 2, 2 ] /( tab [ 2, 1 ] + tab [ 2, 2 ] ) 
sensitivita

# Specificity 
specificita = tab[ 1, 1 ] /( tab [ 1, 2 ] + tab [ 1, 1 ] )
specificita



## 3. Curva ROC

# Costruire la Curva ROC a partire dai valori predetti per la risposta dal modello `mod.low2` 
# dell'analisi della variabile LOWBT.

# Le curve ROC (Receiver Operating Characteristic, anche note come Relative Operating Characteristic) 
# sono degli schemi grafici per un classificatore binario.

# Una curva ROC è il grafico dell'insieme delle coppie (FP, TP) al variare di un parametro del classificatore. 

fit2 = model$fitted


#media campionaria della prob di sopravvivenza nel campione

soglia_roc  = seq( 0, 1, length.out = 2e2 )
lens = length( soglia_roc )-1
ascissa_roc  = rep( NA, lens )
ordinata_roc = rep( NA, lens )

for ( k in 1 : lens )
{
  soglia = soglia_roc [ k ]
  
  classification = as.numeric( sapply( fit2, function( x ) ifelse( x < soglia, 0, 1 ) ) )
  
  #  ATTENZIONE, voglio sulle righe il vero e sulle colonne il predetto
  # t.misc = table( lw$LOW, classification )
  
  ordinata_roc[ k ] = sum( classification[ which( DB_cliniche_senzarispradB$TRG == 1 ) ] == 1 ) /
    length( which( DB_cliniche_senzarispradB$TRG == 1 ) )
  
  ascissa_roc[ k ] = sum( classification[ which( DB_cliniche_senzarispradB$TRG == 0 ) ] == 1 ) /
    length( which( DB_cliniche_senzarispradB$TRG == 0 ) )
  
  # ordinata_roc [ k ]  = t.misc [ 1, 1 ] /( t.misc [ 1, 1 ] + t.misc [ 1, 2 ] )
  #
  # ascissa_roc [ k ]  = t.misc [ 2, 1 ] /( t.misc [ 2, 1 ] + t.misc [ 2, 2 ] )
}


# Visualizziamo la curva ROC.

plot( ascissa_roc, ordinata_roc, type = "l", xlab = "1 - Specificity", ylab = "Sensitivity",
      main = "Curva ROC", lwd = 2, col = 'darkblue', ylim = c( 0, 1 ), xlim = c( 0, 1 ) )
abline( h = c( 0, 1 ), v = c( 0, 1 ), lwd = 1, lty = 2, col = 'red' )
abline( a = 0, b = 1, lty = 2, col = 'black' )

# qual era il nostro punto?
abline( v = 1 - specificita,  h = sensitivita, lty = 3, col = 'blue' )
points( 1 - specificita, sensitivita, pch = 4, lwd = 3, cex = 1.5, col = 'blue' )

# Le linee tratteggiate corrispondono alle due metriche calcolate con la threshold = 0.5 che abbiamo scelto. 

# Attraverso l'analisi delle curve ROC si valuta la capacità del classificatore calcolando 
# l'area sottesa alla curva ROC (Area Under Curve, AUC). 

# R fa tutto in automatico --> calcola AUC e p0ottimale 

PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_senzarispradB$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)