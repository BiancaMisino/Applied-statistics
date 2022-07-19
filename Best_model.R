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




###### MODELLO
DB_cliniche_core_margin<-data.frame(DB_cliniche,DB_core,DB_margin)

DB_cliniche_casoa=DB_cliniche_core_margin
DB_cliniche_casoa[DB_cliniche_casoa[,15]<=3,15] <-0
DB_cliniche_casoa[DB_cliniche_casoa[,15]>3,15] <-1




## ############## LASSO FINO MARGIN CASO A CON RISPOSTA ############################

#variables selection

model<-glm(TRG ~ sync+Diametro_in_mm+OXALIPLATINO+ANTI_VEGF+ANTI_EGFR+Risposta_radiologica+
             HISTO_ExcessKurtosis+GLRLM_SRLGE+M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
             M_SHAPE_Compacity.only.for.3D.ROI..nz.1.+M_GLRLM_LRLGE+M_GLZLM_SZE+M_GLZLM_HGZE+
             primary_tumor_site_left+numerometastasi_1,
           data=DB_cliniche_casoa,family = binomial)

summary(model)

#intervalli di confidenza dei beta
CI<-confint(fit.lasso, level = 0.95)

#Odds ratio
odds.ratio<-matrix(NA,1,model$rank-1)
odds.ratio<-exp(model$coefficients)

model$coefficients
CI
odds.ratio



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
accuracy_train.cv<-matrix(NA,1,k)
auc_train.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model<-glm(TRG~ sync+Diametro_in_mm+OXALIPLATINO+ANTI_VEGF+ANTI_EGFR+Risposta_radiologica+
               HISTO_ExcessKurtosis+GLRLM_SRLGE+M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
               M_SHAPE_Compacity.only.for.3D.ROI..nz.1.+M_GLRLM_LRLGE+M_GLZLM_SZE+M_GLZLM_HGZE+
               primary_tumor_site_left+numerometastasi_1,
             data=DB_cliniche_casoa[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > 0.69, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  accuracy.cv[j]<-cm$overall[1]
  auc.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_casoa[folds==j,]$TRG)
  
}

accuracy.cv
auc.cv

mean(accuracy.cv)  #0.73
mean(auc.cv)        #0.73
mean(accuracy_train.cv)
mean(auc_train.cv)

model_2<-glm(TRG ~sync+ Diametro_in_mm+ANTI_VEGF+Risposta_radiologica+
            GLRLM_SRLGE+M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
            M_GLZLM_SZE+primary_tumor_site_left+numerometastasi_1,
           data=DB_cliniche_casoa,family = binomial)
summary(model_2)

accuracy2.cv<-matrix(NA,1,k)
auc2.cv<-matrix(NA,1,k)

for (j in 1:k) {
  
  model2<-glm(TRG~ sync+Diametro_in_mm+ANTI_VEGF+Risposta_radiologica+
              GLRLM_SRLGE+M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.+
               M_GLZLM_HGZE+
              primary_tumor_site_left+numerometastasi_1,
             data=DB_cliniche_casoa[folds!=j,], family=binomial)
  probabilities <- predict(model,newdata=DB_cliniche_casoa[folds==j,],type = 'response')
  predicted <- ifelse(probabilities > 0.69, "1","0")
  cm<-confusionMatrix(as.factor(predicted),as.factor(DB_cliniche_casoa[folds==j,]$TRG))
  accuracy2.cv[j]<-cm$overall[1]
  auc2.cv[j]<-AUC(as.numeric(predicted),DB_cliniche_casoa[folds==j,]$TRG)
  
}



mean(accuracy2.cv)    #0.74
mean(auc2.cv)        #0.75

#levo volta volta le variabili con pvalue alto:
#MODELLO FINALE model_2 :     Null deviance: 208.63  on 168  degrees of freedom
#                              Residual deviance: 155.02  on 159  degrees of freedom
#                              AIC: 175.02
#n. variabili=9
#   
#  VS 
#MODELLO INIZIALE model:       Null deviance: 208.63  on 168  degrees of freedom
#                              Residual deviance: 151.56  on 153  degrees of freedom
#                              AIC: 183.56
# n.variabili=15   

#pensiamo vada bene model_2 perchè riduciamo la dimensione dell'input space 
#mantenendo circa la stessa bontà del modello

#Odds ratio
odds.ratio<-matrix(NA,1,model$rank-1)
odds.ratio<-exp(model_2$coefficients)

model_2$coefficients


#Se il valore dell'OR è uguale a 1, significa che l'odds di esposizione nei sani è uguale all'odds di esposizione nei malati, cioè il fattore in esame è ininfluente sulla comparsa della malattia.
#Se il valore dell'OR è maggiore di 1, il fattore in esame può essere implicato nella comparsa della malattia (fattore di rischio).
#Se il valore dell'OR è minore di 1 il fattore in esame è una difesa contro la malattia (fattore protettivo).

# le variabili Risposta Radiologica, SHAPE_Sphericity e Primary tumor site left 
#hanno valore dell'odds ratio alto




########################################## MANOVA #####################################

#1) factor TRG:
i0 <- which(DB_cliniche_casoa$TRG==0)
i1 <- which(DB_cliniche_casoa$TRG==1)

data<-DB_cliniche_casoa[,c(15,39,57,69)]
Ps <- c(mcshapiro.test(data[ i0, 2:4])$p,
        mcshapiro.test(data[ i1, 2:4])$p)
Ps    ### P-value:  0.92 0.016

S1 <-  cov(data[i0,2:4])
S2 <-  cov(data[i1,2:4])
x11(width=21)
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
#questo ci piace abbastanza(forse)

#procedo con manova:
fit <- manova(as.matrix(data[,2:4]) ~ data$TRG)
summary.manova(fit,test="Wilks")  #test non significativo
summary.aov(fit)
#in particolare vediamo che vi è una differenza significativa al 10 % nelle variabili GLRLM e M_SHAPE_Sphericity,
#calcolo intervalli di Bonferroni

k <-3 # number of intervals I want to compute (set in advance)
n1<-52
n2<-117
alpha<-0.1
cfr.t <- qt(1-alpha/(2*k),n1+n2-2)
x.mean1<-sapply(data[i1,2:4],mean)
x.mean2<-sapply(data[i0,2:4],mean)
Spooled<-((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
Bf1 <- cbind(inf = (x.mean1[1]-x.mean2[1]) - cfr.t*sqrt(Spooled[1,1]*(1/n1+1/n2)),
             sup = (x.mean1[1]-x.mean2[1]) + cfr.t*sqrt(Spooled[1,1]*(1/n1+1/n2)))
Bf2 <- cbind(inf = (x.mean1[2]-x.mean2[2]) - cfr.t*sqrt(Spooled[2,2]*(1/n1+1/n2)),
             sup = (x.mean1[2]-x.mean2[2]) + cfr.t*sqrt(Spooled[2,2]*(1/n1+1/n2)))
Bf3 <- cbind(inf = (x.mean1[3]-x.mean2[3]) - cfr.t*sqrt(Spooled[3,3]*(1/n1+1/n2)),
             sup = (x.mean1[3]-x.mean2[3]) + cfr.t*sqrt(Spooled[3,3]*(1/n1+1/n2)))
Bf <- rbind(Bf1, Bf2,Bf3)
dimnames(Bf)[[2]] <- c('inf','sup')    
Bf
# INTERVALLI
#                                                inf          sup
#GLRLM_SRLGE                                -0.03589553 -0.006258329    ->non contiene lo zero, torna con manova
#M_SHAPE_Sphericity..only.for.3D.ROI..nz.1. -0.02598625  0.206376770
#M_GLZLM_SZE                                -0.01003626  0.037747411

#2)factor greater_6_cicli
i0 <- which(DB_cliniche_casoa$primary_tumor_site_left==0)
i1 <- which(DB_cliniche_casoa$primary_tumor_site_left==1)

data<-DB_cliniche_casoa[,c(16,39,57)]  #solo core
Ps <- c(mcshapiro.test(data[ i0, 2:3])$p,
        mcshapiro.test(data[ i1, 2:3])$p)
Ps    ### P-value:  0.39 0.73


S1 <-  cov(data[i0,2:3])
S2 <-  cov(data[i1,2:3])
x11(width=21)
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))

#procedo con manova:
fit <- manova(as.matrix(data[,2:3]) ~ data$primary_tumor_site_left)
summary.manova(fit,test="Wilks")  #test non significativo
summary.aov(fit)
#non viene rilevante
#ho provato con interval_greater_30, cicli_greater_6, primary tumor site left e non viene niente



##ANALISI ESPLORATIVA CON STARPLOT SU RADIOMICHE I ORDINE(CONVENTIONAL E HISTO)
x1<-DB_cliniche_casoa[1:169,21:28]
x11()
stars(x1, full = TRUE, scale = TRUE, radius = TRUE,
      labels = dimnames(x)[[1]], locations = NULL,
      nrow = NULL, ncol = NULL, len = 1,
      key.loc = NULL, key.labels = dimnames(x)[[2]],
      key.xpd = TRUE,
      xlim = NULL, ylim = NULL, flip.labels = NULL,
      draw.segments = FALSE,
      col.segments = 1:n.seg, col.stars = NA, col.lines = NA,
      axes = FALSE,
      main = NULL, sub = NULL, xlab = "", ylab = "",
      cex = 0.8, lwd = 0.25, lty = par("lty"), xpd = FALSE,
      add = FALSE, plot = TRUE)
x2<-DB_cliniche_casoa[36:70,21:28]
x11()
stars(x2, full = TRUE, scale = TRUE, radius = TRUE,
      labels = dimnames(x)[[1]], locations = NULL,
      nrow = NULL, ncol = NULL, len = 1,
      key.loc = NULL, key.labels = dimnames(x)[[2]],
      key.xpd = TRUE,
      xlim = NULL, ylim = NULL, flip.labels = NULL,
      draw.segments = FALSE,
      col.segments = 1:n.seg, col.stars = NA, col.lines = NA,
      axes = FALSE,
      main = NULL, sub = NULL, xlab = "", ylab = "",
      cex = 0.8, lwd = 0.25, lty = par("lty"), xpd = FALSE,
      add = FALSE, plot = TRUE)
x3<-DB_cliniche_casoa[71:105,21:28]
x11()
stars(x3, full = TRUE, scale = TRUE, radius = TRUE,
      labels = dimnames(x)[[1]], locations = NULL,
      nrow = NULL, ncol = NULL, len = 1,
      key.loc = NULL, key.labels = dimnames(x)[[2]],
      key.xpd = TRUE,
      xlim = NULL, ylim = NULL, flip.labels = NULL,
      draw.segments = FALSE,
      col.segments = 1:n.seg, col.stars = NA, col.lines = NA,
      axes = FALSE,
      main = NULL, sub = NULL, xlab = "", ylab = "",
      cex = 0.8, lwd = 0.25, lty = par("lty"), xpd = FALSE,
      add = FALSE, plot = TRUE)
x4<-DB_cliniche_casoa[106:140,21:28]
x11()
stars(x4, full = TRUE, scale = TRUE, radius = TRUE,
      labels = dimnames(x)[[1]], locations = NULL,
      nrow = NULL, ncol = NULL, len = 1,
      key.loc = NULL, key.labels = dimnames(x)[[2]],
      key.xpd = TRUE,
      xlim = NULL, ylim = NULL, flip.labels = NULL,
      draw.segments = FALSE,
      col.segments = 1:n.seg, col.stars = NA, col.lines = NA,
      axes = FALSE,
      main = NULL, sub = NULL, xlab = "", ylab = "",
      cex = 0.8, lwd = 0.25, lty = par("lty"), xpd = FALSE,
      add = FALSE, plot = TRUE)
x5<-DB_cliniche_casoa[141:169,21:28]
x11()
stars(x5, full = TRUE, scale = TRUE, radius = TRUE,
      labels = dimnames(x)[[1]], locations = NULL,
      nrow = NULL, ncol = NULL, len = 1,
      key.loc = NULL, key.labels = dimnames(x)[[2]],
      key.xpd = TRUE,
      xlim = NULL, ylim = NULL, flip.labels = NULL,
      draw.segments = FALSE,
      col.segments = 1:n.seg, col.stars = NA, col.lines = NA,
      axes = FALSE,
      main = NULL, sub = NULL, xlab = "", ylab = "",
      cex = 0.8, lwd = 0.25, lty = par("lty"), xpd = FALSE,
      add = FALSE, plot = TRUE)

#per vedere se ci sono clusters simili per forma -> qualcosa si trova


########################### CLUSTERING  #############################################


#SU RADIOMICHE model_2:
data<-DB_cliniche_casoa[,c(39,57,71)]   
summary(data)
 #distance
data.e<-dist(data,method = 'euclidean')
# hierarchical clustering:
HC<- hclust(data.e, method='ward.D2')  #ward.D2
x11()
plot(HC, main='euclidean-Ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='') 

rect.hclust(HC, k=2) 
#cut-tree
cluster.es <- cutree(HC, k=2) 
utilities <- cbind(data, member = cluster.es)
table(utilities$member)
by(data = utilities, INDICES = utilities$member, FUN = summary)
util_summ <- group_by(utilities, member) %>%
summarise(GLRLM_SRLGE=mean(GLRLM_SRLGE),
            M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.=mean(M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.),
            M_GLZLM_HGZE=mean(M_GLZLM_HGZE))
palette(rainbow(3))
# Standardizza i valori in modo che il massimo sia 1
to.draw <- apply(util_summ[, -1], 2, function(x) x/max(x))
# Genera il grafico
x11()
stars(to.draw, draw.segments = TRUE, scale = FALSE, key.loc = c(4.6, 2.3),
      labels = c("CLUSTER 1", "CLUSTER 2"),
      main = "Dati Utility (profili cluster)", cex = .75,
      flip.labels = TRUE)            
cluster.es[which(cluster.es==2)]<-0

table(label.true = t(DB_cliniche_casoa$TRG), label.cluster = cluster.es)


#SU RADIOMICHE model_2 (2 variabili->levo MGLZLM_HGZE):
data<-DB_cliniche_casoa[,c(39,57)]   
summary(data)
#distance
data.e<-dist(data,method = 'euclidean')
# hierarchical clustering:
HC<- hclust(data.e, method='ward.D2')  #ward.D2
x11()
plot(HC, main='euclidean-Ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='') 

rect.hclust(HC, k=2) 
#cut-tree
cluster.es <- cutree(HC, k=2) 
utilities <- cbind(data, member = cluster.es)
table(utilities$member)
by(data = utilities, INDICES = utilities$member, FUN = summary)
util_summ <- group_by(utilities, member) %>%
summarise(GLRLM_SRLGE=mean(GLRLM_SRLGE),
         M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.=mean(M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.))
            
palette(rainbow(2))
# Standardizza i valori in modo che il massimo sia 1
to.draw <- apply(util_summ[, -1], 2, function(x) x/max(x))
# Genera il grafico
x11()
stars(to.draw, draw.segments = TRUE, scale = FALSE, key.loc = c(4.6, 2.3),
      labels = c("CLUSTER 1", "CLUSTER 2"),
      main = "Dati Utility (profili cluster)", cex = .75,
      flip.labels = TRUE)            
cluster.es[which(cluster.es==2)]<-0

table(label.true = t(DB_cliniche_casoa$TRG), label.cluster = cluster.es)
#error rate:

#accuracy:

#precision:
18/(18+37)

#sensitivity


x11()
plot(data , col=cluster.es+1, asp=1, pch=16, lwd=2)
cl<-c(1,2)
legend('bottomleft', legend=cl,fill=c('red','light blue') ,cex=.7)

##SU RADIOMICHE model_2:
data<-DB_cliniche_casoa[,c(39,57,71)]   
summary(data)
#distance
data.e<-dist(data,method = 'euclidean')
# hierarchical clustering:
HC<- hclust(data.e, method='ward.D2')  #ward.D2
x11()
plot(HC, main='euclidean-Ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='') 

rect.hclust(HC, k=2) 
#cut-tree
cluster.es <- cutree(HC, k=2) 
utilities <- cbind(data, member = cluster.es)
table(utilities$member)
by(data = utilities, INDICES = utilities$member, FUN = summary)
util_summ <- group_by(utilities, member) %>%
  summarise(GLRLM_SRLGE=mean(GLRLM_SRLGE),
            M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.=mean(M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.),
            M_GLZLM_HGZE=mean(M_GLZLM_HGZE))
palette(rainbow(3))
# Standardizza i valori in modo che il massimo sia 1
to.draw <- apply(util_summ[, -1], 2, function(x) x/max(x))
# Genera il grafico
x11()
stars(to.draw, draw.segments = TRUE, scale = FALSE, key.loc = c(4.6, 2.3),
      labels = c("CLUSTER 1", "CLUSTER 2"),
      main = "Dati Utility (profili cluster)", cex = .75,
      flip.labels = TRUE)            
cluster.es[which(cluster.es==2)]<-0

table(label.true = t(DB_cliniche_casoa$TRG), label.cluster = cluster.es)
x11()
plot(data , col=cluster.es+1, asp=1, pch=16, lwd=2)
cl<-c(1,2)
legend('bottomleft', legend=cl,fill=c('red','green') ,cex=.7)


#SU RADIOMICHE modello marti (2 variabili->le 2 shape):
data<-DB_cliniche_casoa[,c(30,57)]   
summary(data)
#distance
data.e<-dist(data,method = 'canberra')
# hierarchical clustering:
HC<- hclust(data.e, method='ward.D2')  #ward.D2
x11()
plot(HC, main='euclidean-Ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='') 

rect.hclust(HC, k=2) 
#cut-tree
cluster.es <- cutree(HC, k=2) 
utilities <- cbind(data, member = cluster.es)
table(utilities$member)
by(data = utilities, INDICES = utilities$member, FUN = summary)
util_summ <- group_by(utilities, member) %>%
summarise(SHAPE_Sphericity..only.for.3D.ROI..nz.1.=mean(SHAPE_Sphericity..only.for.3D.ROI..nz.1.),
            M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.=mean(M_SHAPE_Sphericity..only.for.3D.ROI..nz.1.))

palette(rainbow(2))
# Standardizza i valori in modo che il massimo sia 1
to.draw <- apply(util_summ[, -1], 2, function(x) x/max(x))
# Genera il grafico
x11()
stars(to.draw, draw.segments = TRUE, scale = FALSE, key.loc = c(4.6, 2.3),
      labels = c("CLUSTER 1", "CLUSTER 2"),
      main = "Dati Utility (profili cluster)", cex = .75,
      flip.labels = TRUE)            
cluster.es[which(cluster.es==2)]<-0

table(label.true = t(DB_cliniche_casoa$TRG), label.cluster = cluster.es)
#accuracy:
108/169


GLRLM_SRLGE + M_GLZLM_SZE + SHAPE_Volume..mL. + 
  GLRLM_GLNU 

#SU RADIOMICHE modello marti (2 variabili->le 2 shape):

data<-DB_cliniche_casoa[,c(29,39,40,71)]   
summary(data)
#distance
data.e<-dist(data,method = 'canberra')
# hierarchical clustering:
HC<- hclust(data.e, method='ward.D2')  #ward.D2
x11()
plot(HC, main='euclidean-Ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='') 

rect.hclust(HC, k=2) 
#cut-tree
cluster.es <- cutree(HC, k=2) 
utilities <- cbind(data, member = cluster.es)
table(utilities$member)
by(data = utilities, INDICES = utilities$member, FUN = summary)
util_summ <- group_by(utilities, member) %>%
summarise(SHAPE_Volume..mL.=mean(SHAPE_Volume..mL.),
            GLRLM_SRLGE=mean(GLRLM_SRLGE),
            GLRLM_GLNU =mean(GLRLM_GLNU),
            M_GLZLM_HGZE=mean( M_GLZLM_HGZE))

palette(rainbow(4))
# Standardizza i valori in modo che il massimo sia 1
to.draw <- apply(util_summ[, -1], 2, function(x) x/max(x))
# Genera il grafico
x11()
stars(to.draw, draw.segments = TRUE, scale = FALSE, key.loc = c(4.6, 2.3),
      labels = c("CLUSTER 1", "CLUSTER 2"),
      main = "Dati Utility (profili cluster)", cex = .75,
      flip.labels = TRUE)            
cluster.es[which(cluster.es==2)]<-0

table(label.true = t(DB_cliniche_casoa$TRG), label.cluster = cluster.es)
#accuracy:
80/169




