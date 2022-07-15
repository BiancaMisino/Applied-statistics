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
library(glmnet)

library(leaps)
library(ISLR)



Database_finale<-subset(DATABASE_PER_LESION[,-c(1,2,3,5,11,24,25)])


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



#####################################################################################################################
################################ SOLO VARIABILI CORE ######################################################
################################### TRG 0 1 ###################################################################

TRG<-DB_cliniche[,c(15)]
DB_cliniche_core<-data.frame(TRG,DB_core)
rm(TRG)

DB_cliniche_casoa=DB_cliniche_core
DB_cliniche_casoa[DB_cliniche_casoa[,1]<=3,1] <-0
DB_cliniche_casoa[DB_cliniche_casoa[,1]>3,1] <-1
DB_cliniche_casoa$TRG<-factor(DB_cliniche_casoa$TRG,c(0,1),c(0,1))

#variables selection
x <- model.matrix(TRG ~., DB_cliniche_casoa)[,-1]
y <- factor(DB_cliniche_casoa$TRG) 


set.seed(12)
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model<-glm(TRG ~ CONVENTIONAL_HUmin+CONVENTIONAL_HUmean+CONVENTIONAL_HUmax+HISTO_ExcessKurtosis +
             HISTO_Entropy_log2+GLCM_Correlation +GLRLM_GLNU+ NGLDM_Contrast+ GLZLM_LZE +GLZLM_SZLGE +                            
             GLZLM_SZHGE + GLZLM_ZP+ HISTO_Energy...Uniformity.+ SHAPE_Sphericity..only.for.3D.ROI..nz.1.
           +SHAPE_Compacity.only.for.3D.ROI..nz.1.+GLCM_Contrast...Variance. , data=DB_cliniche_casoa,family = binomial)



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

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casoa[1:319,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casoa[320:378,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,378)
folds[1:319]<-folds_mi
folds[320:378]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG ~ CONVENTIONAL_HUmin+CONVENTIONAL_HUmean+CONVENTIONAL_HUmax+HISTO_ExcessKurtosis +
               HISTO_Entropy_log2+GLCM_Correlation +GLRLM_GLNU+ NGLDM_Contrast+ GLZLM_LZE +GLZLM_SZLGE +                            
               GLZLM_SZHGE + GLZLM_ZP+ HISTO_Energy...Uniformity.+ SHAPE_Sphericity..only.for.3D.ROI..nz.1.
             +SHAPE_Compacity.only.for.3D.ROI..nz.1.+GLCM_Contrast...Variance., data=DB_cliniche_casoa[folds!=j,],family = binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.72, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  ##NON FUNZIA AUC
  #auc.cv[j]<-AUC(predicted,DB_cliniche_casoa[folds==j,]$TRG)
}


accuracy.cv
auc.cv
model$coefficients
mean(accuracy.cv)
mean(auc.cv)

#######################################################CASO B#####################################################

####secondo caso dicotomico TRG ={0,1} (1,3+5)
DB_cliniche_casob=DB_cliniche_core
DB_cliniche_casob[DB_cliniche_casob[,11]<3,1] <-0
DB_cliniche_casob[DB_cliniche_casob[,11]>=3,1] <-1


DB_cliniche_casoa$TRG<-factor(DB_cliniche_casoa$TRG,c(0,1),c(0,1))

#variables selection
x <- model.matrix(TRG ~., DB_cliniche_casoa)[,-1]
y <- factor(DB_cliniche_casoa$TRG) 


set.seed(12)
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model<-glm(TRG ~ CONVENTIONAL_HUmin+CONVENTIONAL_HUmean+CONVENTIONAL_HUmax+HISTO_ExcessKurtosis +
             HISTO_Entropy_log2+GLCM_Correlation +GLRLM_GLNU+ NGLDM_Contrast+ GLZLM_LZE +GLZLM_SZLGE +                            
             GLZLM_SZHGE + GLZLM_ZP+ HISTO_Energy...Uniformity.+ SHAPE_Sphericity..only.for.3D.ROI..nz.1.
           +SHAPE_Compacity.only.for.3D.ROI..nz.1.+GLCM_Contrast...Variance. , data=DB_cliniche_casoa,family = binomial)

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

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casoa[1:319,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casoa[320:378,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,378)
folds[1:319]<-folds_mi
folds[320:378]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG ~ CONVENTIONAL_HUmin+CONVENTIONAL_HUmean+CONVENTIONAL_HUmax+HISTO_ExcessKurtosis +
               HISTO_Entropy_log2+GLCM_Correlation +GLRLM_GLNU+ NGLDM_Contrast+ GLZLM_LZE +GLZLM_SZLGE +                            
               GLZLM_SZHGE + GLZLM_ZP+ HISTO_Energy...Uniformity.+ SHAPE_Sphericity..only.for.3D.ROI..nz.1.
             +SHAPE_Compacity.only.for.3D.ROI..nz.1.+GLCM_Contrast...Variance., data=DB_cliniche_casoa[folds!=j,],family = binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.84, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  ##NON FUNZIA AUC
  #auc.cv[j]<-AUC(predicted,DB_cliniche_casoa[folds==j,]$TRG)
}


accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)



#####################################################################################################################
################################ SOLO VARIABILI CORE+MARGIN ################################################
################################### TRG 0 1#########################################################################

TRG<-DB_cliniche[,c(15)]

DB_cliniche_core_margin<-data.frame(TRG,DB_core,DB_margin)

DB_cliniche_casoa=DB_cliniche_core_margin
DB_cliniche_casoa[DB_cliniche_casoa[,1]<=3,1] <-0
DB_cliniche_casoa[DB_cliniche_casoa[,1]>3,1] <-1
DB_cliniche_casoa$TRG<-factor(DB_cliniche_casoa$TRG,c(0,1),c(0,1))


rm(TRG)



#variables selection
x <- model.matrix(TRG ~., DB_cliniche_casoa)[,-1]
y <- factor(DB_cliniche_casoa$TRG) 


set.seed(12)
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model<-glm(TRG ~ CONVENTIONAL_HUmean + CONVENTIONAL_HUmax + HISTO_ExcessKurtosis + HISTO_Entropy_log2 +
             HISTO_Energy...Uniformity.+ SHAPE_Volume..mL. + SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
             SHAPE_Compacity.only.for.3D.ROI..nz.1.+ GLCM_Entropy_log2...Joint.entropy.+ GLRLM_LRE +
             GLRLM_GLNU + NGLDM_Contrast + GLZLM_LZE + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_ZP + M_CONVENTIONAL_HUmin+
             M_CONVENTIONAL_HUmean + M_CONVENTIONAL_HUstd +M_CONVENTIONAL_HUmax + M_HISTO_Energy...Uniformity.+
             M_SHAPE_Sphericity..only.for.3D.ROI..nz.1. + M_GLCM_Homogeneity...Inverse.difference.+
             M_GLCM_Correlation + M_GLRLM_SRE + M_GLRLM_LGRE + M_NGLDM_Coarseness + M_NGLDM_Contrast +
             M_GLZLM_SZE + M_GLZLM_LZE + M_GLZLM_HGZE + M_GLZLM_SZHGE + M_GLZLM_GLNU + M_GLZLM_ZP
           ,data=DB_cliniche_casoa,family = binomial)
CI<-confint(model, level = 0.95)
odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))
model$coefficients
odds.ratio
CI
fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_casoa$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casoa[1:319,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casoa[320:378,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,378)
folds[1:319]<-folds_mi
folds[320:378]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG ~ CONVENTIONAL_HUmean + CONVENTIONAL_HUmax + HISTO_ExcessKurtosis + HISTO_Entropy_log2 +
               HISTO_Energy...Uniformity.+ SHAPE_Volume..mL. + SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
               SHAPE_Compacity.only.for.3D.ROI..nz.1.+ GLCM_Entropy_log2...Joint.entropy.+ GLRLM_LRE +
               GLRLM_GLNU + NGLDM_Contrast + GLZLM_LZE + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_ZP + M_CONVENTIONAL_HUmin+
               M_CONVENTIONAL_HUmean + M_CONVENTIONAL_HUstd +M_CONVENTIONAL_HUmax + M_HISTO_Energy...Uniformity.+
               M_SHAPE_Sphericity..only.for.3D.ROI..nz.1. + M_GLCM_Homogeneity...Inverse.difference.+
               M_GLCM_Correlation + M_GLRLM_SRE + M_GLRLM_LGRE + M_NGLDM_Coarseness + M_NGLDM_Contrast +
               M_GLZLM_SZE + M_GLZLM_LZE + M_GLZLM_HGZE + M_GLZLM_SZHGE + M_GLZLM_GLNU + M_GLZLM_ZP
             ,data=DB_cliniche_casoa,family = binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.72, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  #auc.cv[j]<-AUC(predicted,DB_cliniche_casoa[folds==j,]$TRG)
}


accuracy.cv
auc.cv
model$coefficients
mean(accuracy.cv)
mean(auc.cv)

#######################################################CASO B#####################################################

####secondo caso dicotomico TRG ={0,1} (1,3+5)
DB_cliniche_casob=DB_cliniche_core_margin
DB_cliniche_casob[DB_cliniche_casob[,1]<3,1] <-0
DB_cliniche_casob[DB_cliniche_casob[,1]>=3,1] <-1


DB_cliniche_casoa$TRG<-factor(DB_cliniche_casoa$TRG,c(0,1),c(0,1))


#variables selection
x <- model.matrix(TRG ~., DB_cliniche_casoa)[,-1]
y <- factor(DB_cliniche_casoa$TRG) 


set.seed(12)
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables

model<-glm(TRG ~ CONVENTIONAL_HUmean + CONVENTIONAL_HUmax + HISTO_ExcessKurtosis + HISTO_Entropy_log2 +
             HISTO_Energy...Uniformity.+ SHAPE_Volume..mL. + SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
             SHAPE_Compacity.only.for.3D.ROI..nz.1.+ GLCM_Entropy_log2...Joint.entropy.+ GLRLM_LRE +
             GLRLM_GLNU + NGLDM_Contrast + GLZLM_LZE + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_ZP + M_CONVENTIONAL_HUmin+
             M_CONVENTIONAL_HUmean + M_CONVENTIONAL_HUstd +M_CONVENTIONAL_HUmax + M_HISTO_Energy...Uniformity.+
             M_SHAPE_Sphericity..only.for.3D.ROI..nz.1. + M_GLCM_Homogeneity...Inverse.difference.+
             M_GLCM_Correlation + M_GLRLM_SRE + M_GLRLM_LGRE + M_NGLDM_Coarseness + M_NGLDM_Contrast +
             M_GLZLM_SZE + M_GLZLM_LZE + M_GLZLM_HGZE + M_GLZLM_SZHGE + M_GLZLM_GLNU + M_GLZLM_ZP
           ,data=DB_cliniche_casoa,family = binomial)
CI<-confint(model, level = 0.95)
odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))
model$coefficients
odds.ratio
CI
fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_casoa$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

k <- 5

set.seed(1)
folds_mi <- sample(1:k,nrow(DB_cliniche_casoa[1:319,]),replace=TRUE)
folds_mi
folds_to <- sample(1:k,nrow(DB_cliniche_casoa[320:378,]),replace=TRUE)
folds_to
table(folds_to)
folds<-matrix(NA,1,378)
folds[1:319]<-folds_mi
folds[320:378]<-folds_to

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

accuracy.cv<-matrix(NA,1,k)
auc.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG ~ CONVENTIONAL_HUmean + CONVENTIONAL_HUmax + HISTO_ExcessKurtosis + HISTO_Entropy_log2 +
               HISTO_Energy...Uniformity.+ SHAPE_Volume..mL. + SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
               SHAPE_Compacity.only.for.3D.ROI..nz.1.+ GLCM_Entropy_log2...Joint.entropy.+ GLRLM_LRE +
               GLRLM_GLNU + NGLDM_Contrast + GLZLM_LZE + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_ZP + M_CONVENTIONAL_HUmin+
               M_CONVENTIONAL_HUmean + M_CONVENTIONAL_HUstd +M_CONVENTIONAL_HUmax + M_HISTO_Energy...Uniformity.+
               M_SHAPE_Sphericity..only.for.3D.ROI..nz.1. + M_GLCM_Homogeneity...Inverse.difference.+
               M_GLCM_Correlation + M_GLRLM_SRE + M_GLRLM_LGRE + M_NGLDM_Coarseness + M_NGLDM_Contrast +
               M_GLZLM_SZE + M_GLZLM_LZE + M_GLZLM_HGZE + M_GLZLM_SZHGE + M_GLZLM_GLNU + M_GLZLM_ZP
             ,data=DB_cliniche_casoa,family = binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  #andrebbe messo la frequenza di trg=1
  predicted <- ifelse(probabilities > 0.84, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  #auc.cv[j]<-AUC(predicted,DB_cliniche_casoa[folds==j,]$TRG)
}


accuracy.cv
auc.cv

mean(accuracy.cv)
mean(auc.cv)



#####################################################################################################################
################################ SOLO VARIABILI CLINICHE+CORE+MARGIN ################################################
################################### TRG 0 1#########################################################################

#con risposta radiologica
DB_core<-DB_core[,c(2,4,6,7,8,9,10,11,15,17,20,21,23,24,26,27)]
DB_margin<-DB_margin[,c(1,2,3,4,8,10,12,13,15,16,19,20,22,23,24,25,26,27)]


##creo il mio database aggiungendo le clinche dopo la selezione attraverso modello lasso
DB<-data.frame(DB_cliniche,DB_core,DB_margin)

DB_cliniche_casoa=DB
DB_cliniche_casoa[DB_cliniche_casoa[,15]<=3,15] <-0
DB_cliniche_casoa[DB_cliniche_casoa[,15]>3,15] <-1
DB_cliniche_casoa$TRG<-factor(DB_cliniche_casoa$TRG,c(0,1),c(0,1))

#variables selection
x <- model.matrix(TRG ~., DB_cliniche_casoa)[,-1]
y <- factor(DB_cliniche_casoa$TRG) 


set.seed(12)
cv.lasso <- cv.glmnet(x, y, family = "binomial",type.measure = "auc")
c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables
##sono decisamente troppe ... concludiamo passando da 124 a 52

model<-glm(TRG ~ Interval._greater_30.days+ Age + Sex + sync+ Diametro_in_mm+ CEA+ bilater+ OXALIPLATINO+ IRINOTECAN + 
             ANTI_VEGF+ANTI_EGFR+ Linee_di_chemioterapia+ greater_6_cicli+Risposta_radiologica+primary_tumor_site_left+primary_tumor_site_right+
             numerometastasi_1+ numerometastasi_2to3+ numerometastasi_4to9+ CONVENTIONAL_HUmean + CONVENTIONAL_HUmax + HISTO_ExcessKurtosis + HISTO_Entropy_log2 +
             HISTO_Energy...Uniformity.+ SHAPE_Volume..mL. + SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
             SHAPE_Compacity.only.for.3D.ROI..nz.1.+ GLCM_Entropy_log2...Joint.entropy.+ GLRLM_LRE +
             GLRLM_GLNU + NGLDM_Contrast + GLZLM_LZE + GLZLM_LGZE + GLZLM_SZHGE + GLZLM_ZP + M_CONVENTIONAL_HUmin+
             M_CONVENTIONAL_HUmean + M_CONVENTIONAL_HUstd +M_CONVENTIONAL_HUmax + M_HISTO_Energy...Uniformity.+
             M_SHAPE_Sphericity..only.for.3D.ROI..nz.1. + M_GLCM_Homogeneity...Inverse.difference.+
             M_GLCM_Correlation + M_GLRLM_SRE + M_GLRLM_LGRE + M_NGLDM_Coarseness + M_NGLDM_Contrast +
             M_GLZLM_SZE + M_GLZLM_LZE + M_GLZLM_HGZE + M_GLZLM_SZHGE + M_GLZLM_GLNU + M_GLZLM_ZP
           ,data=DB_cliniche_casoa,family = binomial)
CI<-confint(model, level = 0.95)
odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio<-matrix(NA,1,model$rank - 1)
odds.ratio=exp(coef(model))
model$coefficients
odds.ratio
CI
fit2<-model$fitted
PRROC_obj <- roc.curve(scores.class0 = fit2, weights.class0=as.numeric(paste(DB_cliniche_casoa$TRG)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

vif(model)
