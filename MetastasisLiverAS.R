
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

##############################  MANOVA ETA' ##############################################
##proviamo a fare MANOVA dividendo in due gruppi per età:0-65 vs 65-82 anni con variabili p variabili radiomiche dei vari tipi

DB_lesmag_age<-subset(Database_TRG_anonymized_06_03_22[,c(1,7,29:76)])
DB_rad_age <-na.omit(DB_lesmag_age)
DB_rad_age<-DB_rad_age[,-c(7,8,9,11,13,23,32,33,35,36,46,47,49)]

###CONVENTIONAL

DB_lesmag_age<-subset(Database_TRG_anonymized_06_03_22[,c(1,7,29:35)])
DB_rad_age <-na.omit(DB_lesmag_age)

summary(DB_rad_age)

var.risp <-DB_rad_age [,3:9]
group.names <- DB_rad_age[,2]
dim(var.risp)[2]
levels(group.names)

p<-7
g<-2

i1 <- which(group.names<65)
i2 <- which(group.names>=65)

ng <- c(length(i1),length(i2))
ng
N <- sum(ng)

# verifichiamo Gaussianità e omogeneità della matrice covarianza
# N.B. da caricare funzione mcshapiro da lab5 sennò non parte nulla!!

#1)Gaussianità
Ps <- c(mcshapiro.test(var.risp[ i1, ])$p,
        mcshapiro.test(var.risp[ i2, ])$p)
Ps
  ###non gaussiane infatti vediamo istogrammi e qqplot non carini

x11()
par(mfrow=c(2,4))

hist(as.matrix(var.risp[,1]), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.matrix(var.risp[,1])),sd(as.matrix(var.risp[,1]))), col='blue', lty=2)

hist(as.matrix(var.risp[,2]), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.matrix(var.risp[,2])),sd(as.matrix(var.risp[,2]))), col='blue', lty=2)

hist(as.matrix(var.risp[,3]), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.matrix(var.risp[,3])),sd(as.matrix(var.risp[,3]))), col='blue', lty=2)

hist(as.matrix(var.risp[,4]), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.matrix(var.risp[,4])),sd(as.matrix(var.risp[,4]))), col='blue', lty=2)

hist(as.matrix(var.risp[,5]), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.matrix(var.risp[,5])),sd(as.matrix(var.risp[,5]))), col='blue', lty=2)

hist(as.matrix(var.risp[,6]), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.matrix(var.risp[,6])),sd(as.matrix(var.risp[,6]))), col='blue', lty=2)

hist(as.matrix(var.risp[,7]), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.matrix(var.risp[,7])),sd(as.matrix(var.risp[,7]))), col='blue', lty=2)

x11()
par(mfrow=c(2,4))

qqnorm(as.matrix(var.risp[,1]), main='QQplot of x')
qqline(as.matrix(var.risp[,1]))

qqnorm(as.matrix(var.risp[,2]), main='QQplot of y')
qqline(as.matrix(var.risp[,2]))

qqnorm(as.matrix(var.risp[,3]), main='QQplot of y')
qqline(as.matrix(var.risp[,3]))

qqnorm(as.matrix(var.risp[,4]), main='QQplot of y')
qqline(as.matrix(var.risp[,4]))

qqnorm(as.matrix(var.risp[,5]), main='QQplot of y')
qqline(as.matrix(var.risp[,5]))

qqnorm(as.matrix(var.risp[,6]), main='QQplot of y')
qqline(as.matrix(var.risp[,6]))

qqnorm(as.matrix(var.risp[,7]), main='QQplot of y')
qqline(as.matrix(var.risp[,7]))

##alcune fanno cahà (variabili 1,3,4) quindi proviamo boxcox trasformation ma non prende il comando powerTransform(var.risp[,1:7])
##le variabili 1 e 5 di var.risp hanno dati negativi.
##il comando Powertrasformation funziona solo con dati positivi.
##quindi abbiamo provato trasformazioni logaritmiche sulle variabili 3,4 ma gli shapiro vengo sempre bassissimi:

#variabile 3 ->trasf log
y_3<-matrix(0,127,1)
for (i in 1:127)
  {if(var.risp[i,3]<0)
  {
    y_3[i]<- -log(-var.risp[i,3]+1)
  }
  else
  {
    y_3[i]<-log(var.risp[i,3]+1)
  }
}
x11()
par(mfrow=c(1,2))
hist(as.numeric(y_3), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.numeric(y_3)),sd(as.numeric(y_3))), col='blue', lty=2) 
#caruccio
qqnorm(as.numeric(y_3), main='QQplot of y')
qqline(as.numeric(y_3))

shapiro.test(as.numeric(y_3))


##var4 ->trasf log
y_4<-matrix(0,127,1)
for (i in 1:127)
{if(var.risp[i,4]<0)
{
  y_4[i]<- -log(-var.risp[i,4]+1)
}
  else
  {
    y_4[i]<-log(var.risp[i,4]+1)
  }
}
x11()
par(mfrow=c(1,2))
hist(as.numeric(y_4), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.numeric(y_4)),sd(as.numeric(y_4))), col='blue', lty=2)
#caruccio
qqnorm(as.numeric(y_4), main='QQplot of y')
qqline(as.numeric(y_4))

shapiro.test(as.numeric(y_4))

##visto che non si risolveva niente abbiamo provato a usare power trasf sulle variabili con dati positivi(quindi senza var 1 e 5)

lambda.mult <- powerTransform(var.risp[,c(2,3,4,6,7)])    
lambda.mult

#vediamo che solo lambda[2] è lontano da 1 -->trasformiamo solo la variabile CONVENTIONAL_HUstd (var.risp[2])
BC.var3 <- bcPower(var.risp[,3],  lambda.mult$lambda[2]) 
x11()
par(mfrow=c(1,2))
hist(as.matrix(BC.var3), prob=T,col='grey85')
lines(0:1000 / 100, dnorm(0:1000 / 100,mean(as.numeric(BC.var3)),sd(as.numeric(BC.var3))), col='blue', lty=2)
qqnorm(as.numeric(BC.var3), main='QQplot of y')
qqline(as.numeric(BC.var3))

shapiro.test(as.matrix(BC.var3))  ##shapiro ancora bruttino
##sostituiamo la var 3 con quella trasformata:

var.risp[,3]<-BC.var3

#proviamo a ricalcolare lo shapiro multivariato ma non cambia nulla e p_value sempre (0,0)

Ps <- c(mcshapiro.test(var.risp[ i1, c(2,3,4,6,7) ])$p,
        mcshapiro.test(var.risp[ i2, c(2,3,4,6,7)])$p)
Ps
###no-->abbiamo mandato mail a arnone  





#2)Omogeneità di varianza
S  <-  cov(var.risp)
S1 <-  cov(var.risp[i1,])
S2 <-  cov(var.risp[i2,])

# Qualitatively:
round(S1,digits=1)
round(S2,digits=1)


x11(width=21)
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
##questo ci piace abbastanza(forse)

#test MANOVA

fit <- manova(as.matrix(var.risp) ~ as.matrix(group.names))
summary.manova(fit,test="Wilks")

##pvalue alto -->sembra che l'età non abbia effetto sulla media delle radiomiche


####DA ANDARE AVANTI CON MANOVA PER OGNI GRUPPO DI RADIOMICHE




##############################  MANOVA CICLI DI CHEMIO ##############################################
############################## DA RIFARE  #######################################

DB_lesmag_chemio<-subset(Database_TRG_anonymized_06_03_22[,c(1,24,29:76)])
DB_rad_chemio <-na.omit(DB_lesmag_chemio)

DB_rad_chemio<-DB_rad_chemio[,-c(7,8,9,11,13,23,32,33,35,36,46,47,49)]


summary(DB_rad_chemio)

var.risp <-DB_rad_chemio [,3:37]
group.names <- DB_rad_chemio[,2]
dim(var.risp)[2]
levels(group.names)

p<-35
g<-2

i1 <- which(group.names<9)
i2 <- which(group.names>=9)

ng <- c(length(i1),length(i2))
ng
N <- sum(ng)

#per ora tutto va bene menomale, c'era ansia, andiamo avanti e verifichiamo Gaussianità e omogeneità della matrice covarianza
#N.B. da caricare funzione mcshapiro da lab5 sennò non parte nulla!!

#1)Gaussianità
Ps <- c(mcshapiro.test(var.risp[ i1, ])$p,
        mcshapiro.test(var.risp[ i2, ])$p)
Ps
### problema: non funziona lo shapiro perchè la matrice cov è singolare perchè sono correlate le features

#2)Omogeneità di varianza
S  <-  cov(var.risp)
S1 <-  cov(var.risp[i1,])
S2 <-  cov(var.risp[i2,])

# Qualitatively:
round(S1,digits=1)
round(S2,digits=1)


x11(width=21)
par(mfrow=c(1,2))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
##questo ci piace abbastanza(forse)

#test MANOVA

fit <- manova(as.matrix(var.risp) ~ as.matrix(group.names))
summary.manova(fit,test="Wilks")





















