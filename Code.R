
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



#importo dataset

#View(Database_TRG_anonymized_06_03_22)
summary(Database_TRG_anonymized_06_03_22) 
dim(Database_TRG_anonymized_06_03_22)


#solo variabili cliniche fino al TRG (incluso)
varcli <- subset(Database_TRG_anonymized_06_03_22[,1:28])

########################### CASO A ####################################
#considero solo la lesione più grande per ogni paziente:
#essendo già in ordine decrescente elimino tutti gli NA sotto l'id del paziente,
#tengo solo la prima riga per ogni paziente (check VOI == 1)

DB_cliniche <- varcli[complete.cases(varcli[,1]),]

summary(DB_cliniche)
dim(DB_cliniche)

#n <- dim(varcliniche)[1]
#p <- dim(varcliniche)[2]
#DB_cliniche<-na.omit(varcliniche) #elimina tutte le righe con almeno un NA

#solo variabili cliniche numeriche (escludo anche codice e id dei pazienti)

DB_cliniche_num<-DB_cliniche[,-c(1,2,3,4,6,8,9,10,12,15:23,25:28)]

#DATA CLEANING
#exploration
x11()
pairs(DB_cliniche_num)

#da questo primo grafico si può notare un comportamento strano della 
#variabile CEA analizzo più in dettaglio con un boxplot

x11()
boxplot(DB_cliniche_num)

#x11()
#boxplot(scale(x=DB_cliniche_num,center = F, scale=F), las=2, col='gold',outline=FALSE)

# come ipoteizzato -> CEA presenta un valore non ammissibile (>4000)
# quindi va eliminata (riga 45)

NDB <- DB_cliniche_num[-45,]


x11()
boxplot(scale(x=NDB,center = F, scale=F), las=2, col='gold')

# noto la presenza di altri ouliers (CEA) ma non li elimino poichè
# in totale rappresentano più del 20% dei pazienti (troppo)

# eliminare outlier tramite IQR
#Q <- quantile(DB_cliniche_num$CEA, probs=c(.25, .75), na.rm = FALSE)
#iqr <- IQR(DB_cliniche_num$CEA)
#up <-  Q[2]+1.5*iqr # Upper Range  
#low<- Q[1]-1.5*iqr # Lower Range

#eliminated<- subset(DB_cliniche_num, DB_cliniche_num$CEA > (Q[1] - 1.5*iqr) & DB_cliniche_num$CEA < (Q[2]+1.5*iqr)
#DB_cliniche_num(eliminated, Interval CT-Surgery,days breaks, outlier.tagging = TRUE) 

# PCA on NDB dataset

picia <- princomp(NDB,cor=TRUE,scores = TRUE)
summary(picia)

# le prime 4 componenti spiegano l'80% della variabilità


# loadings ( coefficients of the linear combination of the original 
#           variables that defines each principal component)
load.tour <- picia$loadings
load.tour

load.tour[,1:6]

# graphical representation of the loadings of the first six principal components

x11()
par(mfcol = c(2,3))
for(i in 1:6) barplot(load.tour[,i], ylim = c(-1, 1), main=paste("PC",i))


##Explained variance
#standardizzo
NDB.sd <- scale(NDB)
NDB.sd <- data.frame(NDB.sd)

x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(picia, las=2, main='Principal Components', ylim=c(0,3))
abline(h=1, col='blue')
barplot(sapply(NDB.sd,sd)^2, las=2, main='Original Variables', ylim=c(0,3), ylab='Variances')
plot(cumsum(picia$sde^2)/sum(picia$sde^2), type='b', axes=F, xlab='Number of components', ylab='Contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(NDB.sd),labels=1:ncol(NDB.sd),las=2)

#scores NDB
scores.NDB <- picia$scores

x11()
plot(scores.NDB[,1:4])
abline(h=0, v=0, lty=2, col='grey')

x11()
layout(matrix(c(1,2),2))
boxplot(NDB.sd, las=2, col='gold', main='Standardized variables')
scores.NDB <- data.frame(scores.NDB)
boxplot(scores.NDB, las=2, col='gold', main='Principal components')

x11()
biplot(picia)

graphics.off()


############ ANALISI VARIABILI RADIOMICHE #############

#come per le clinche considero solamente la prima lesione (più grande)
DB_lesmag<-subset(Database_TRG_anonymized_06_03_22[,c(1,29:76)])
DB_rad <-na.omit(DB_lesmag)
DB_rad <- DB_rad[-45,]# CEA 4216 elimino
summary(DB_rad)

#elimino ID pazienti
DB_rad<-DB_rad[,-1]
dim(DB_rad)

# esplorazione grafica di tutte le radiomiche * gruppo indicatore
# per analizzare le correlazioni nel gruppo 

#CONVENTIONAL
x11()
pairs(DB_rad[,1:7]) 
#tolgo tutti i quartili 5,6,7 (tengo la media)

#HISTO
x11()
pairs(DB_rad[,8:13]) 
# tolgo curtosis 9 (tengo excesscurtosis)
# tolgo entropy_log10 11 (tengo log2)
# energy correlata neg con entropy, togliere una delle due?

#SHAPE
x11()
pairs(DB_rad[,14:16])  
#controllare relazione tra shape_volume e compattezza

#GLCM
x11()
pairs(DB_rad[,17:23])
#tolgo entropy_log10 21 (tengo log2)

#GLRLM
x11()
pairs(DB_rad[,24:34]) 
#tolgo LRLGE e LRLGHE  30/31 (tengo LRE)
#tolgo 33 (tengo GLNU)
#tolgo 34 (tengo SRE)

#NGLDM
x11()
pairs(DB_rad[,35:37]) 

#GLZLM
x11()
pairs(DB_rad[,38:48]) 
#tolgo LZLGE e LZHGE 44/45 (tengo LZE)
#tolgo ZLNU 47 (tengo GLNU)

#considero solo le variabili selezionate

DB_rad<-DB_rad[,-c(5,6,7,9,11,21,30,31,33,34,44,45,47)]


#PCA radiomiche
rad_picia <- princomp(DB_rad,cor=TRUE,scores = TRUE)
summary(rad_picia)  
#con 5/6/7 principal components ho 82/88/91% della variability spiegato
#-> PCA molto efficace : da 35 variabili a 5/6/7!!

#standardizzo (unità di misura)
DB_rad.sd <- scale(DB_rad)
DB_rad.sd <- data.frame(DB_rad.sd)

x11()
layout(matrix(c(1,2,1,2),2,byrow=T))
plot(rad_picia, las=2, main='Principal Components', ylim=c(0,12))
abline(h=1, col='blue')
plot(cumsum(rad_picia$sde^2)/sum(rad_picia$sde^2), type='b', axes=F, xlab='Number of components', ylab='Contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(DB_rad.sd),labels=1:ncol(DB_rad.sd),las=2)

load.tour <- rad_picia$loadings
load.tour

# graphical representation of the loadings of the first seven principal components
x11()
barplot(load.tour[,1], ylim = c(-1, 1), main=paste("PC",1))## bars proportional to loadings 
x11()
barplot(load.tour[,2], ylim = c(-1, 1), main=paste("PC",2))## bars proportional to loadings 
x11()
barplot(load.tour[,3], ylim = c(-1, 1), main=paste("PC",3))## bars proportional to loadings 
x11()
barplot(load.tour[,4], ylim = c(-1, 1), main=paste("PC",4))## bars proportional to loadings 
x11()
barplot(load.tour[,5], ylim = c(-1, 1), main=paste("PC",5))## bars proportional to loadings 
x11()
barplot(load.tour[,6], ylim = c(-1, 1), main=paste("PC",6))## bars proportional to loadings 
x11()
barplot(load.tour[,7], ylim = c(-1, 1), main=paste("PC",7))## bars proportional to loadings 


