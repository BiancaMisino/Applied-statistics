setwd("C:/Users/sarar/Desktop/POLI/APPLIED STATISTIC/Project/Applied-statistics-main")

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
library(readxl)
library(class)


################################################################################
################################################################################
##                      SUPERVISED CLASSIFICATION



##-------------------------------------------------------------------------------
## CLINICHE - TRG (0,1)
##     TRG=0 se era 1,3
##     TRG=1 se era 5

#importo dataset e funzione shapiro
Database_TRG_anonymized_06_03_22 <- read_excel("Database TRG anonymized 06_03_22.xlsx")
load("C:/Users/sarar/Desktop/POLI/APPLIED STATISTIC/Labs/Lab 5/mcshapiro.test.RData")

#solo variabili cliniche fino al TRG (incluso)-lesione più grande
varcli <- subset(Database_TRG_anonymized_06_03_22[,1:28])
DB_cliniche <- varcli[complete.cases(varcli[,1]),]

#database con variabili cliniche numeriche e il TRG -elimino l'outlier
DBC<-DB_cliniche[-45,-c(1,2,3,4,6,8,9,10,12,15:23,25:27)]
head(DBC)
dim(DBC)

#modifico il dataset in modo da avere TRG= 1 o 0
DBC[DBC[,7]<=3, 7]<-0
DBC[DBC[,7]>=5, 7]<-1

attach(DBC)

trg.val <- factor(TRG, labels=c('0','1'))
trg.val
g=2 

trg1 <- which(trg.val== '0')
trg2 <- which(trg.val=='1')
trg1
trg2

n1 <- length(trg1)
n2 <- length(trg2)
n <- n1+n2

detach(DBC)


#SCELGO SOLO DUE VARIABILI PER LA CLASSIFICAZIONE

# plot the data
x11()
plot(DBC, col=trg.val)


# tengo: 1 intervallo e  2 età
DBC2 <- DBC[,c(4,5)] #12 
head(DBC2)

# plot the data
x11()
plot(DBC2, main='VAR1 e VAR2', xlab='VAr1', ylab='VAR2', pch=20)
points(DBC2[trg1,], col='red', pch=20)
points(DBC2[trg2,], col='green', pch=20)
legend("topright", legend=levels(trg.val), fill=c('red','green'))

m <-  colMeans(DBC2)
m1 <- colMeans(DBC2[trg1,])
m2 <- colMeans(DBC2[trg2,])

S1 <- cov(DBC2[trg1,])
S2 <- cov(DBC2[trg2,])
Sp  <- ((n1-1)*S1+(n2-1)*S2)/(n-g)

# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) sigma.A=sigma.B (SOLO PER LDA)
# 3) c(A|B)=c(B|A) (equal misclassification costs)

# verify assumptions: normality within the groups 

mcshapiro.test(DBC2[trg1,])
mcshapiro.test(DBC2[trg2,])

##------------------------------------------------
# per verificare le ipotesi rendo normali i dati 
# separatamente in base all'etichetta
# secondo me questa cosa non ha senso

lambda.mult1 <- powerTransform(DBC2[trg1,])    
lambda.mult1

tr.11 <- bcPower(DBC2[trg1,1],  lambda.mult1$lambda[1]) 
tr.12 <- bcPower(DBC2[trg1,2],  lambda.mult1$lambda[2]) 
DBC2[trg1,1]<-tr.11
DBC2[trg1,2]<-tr.12
mcshapiro.test(DBC2[trg1,])

lambda.mult2 <- powerTransform(DBC2[trg2,])    
lambda.mult2

tr.21 <- bcPower(DBC2[trg2,1],  lambda.mult2$lambda[1]) 
tr.22 <- bcPower(DBC2[trg2,2],  lambda.mult2$lambda[2]) 
DBC2[trg2,1]<-tr.21
DBC2[trg2,2]<-tr.22
mcshapiro.test(DBC2[trg2,])


#------------------
#dato che quello fatto prima per me non ha senso lo faccio

lambda.mult <- powerTransform(DBC2[,])    
lambda.mult

tr.1 <- bcPower(DBC2[,1],  lambda.mult$lambda[1]) 
tr.2 <- bcPower(DBC2[,2],  lambda.mult$lambda[2]) 
DBC2[,1]<-tr.1
DBC2[,2]<-tr.2
mcshapiro.test(DBC2[trg1,])
mcshapiro.test(DBC2[trg2,])


### Quadratic Discriminand Analysis (QDA)
###---------------------------------------

qda.dbc2 <- qda(DBC2, trg.val)
qda.dbc2

Qda.dbc2 <- predict(qda.dbc2, DBC2)

# compute the APER
Qda.dbc2$class
trg.val
table(class.true=trg.val, class.assigned=Qda.dbc2$class)

errorsq <- (Qda.dbc2$class != trg.val)
errorsq

APERq   <- sum(errorsq)/length(trg.val)
APERq


# Compute the estimate of the AER by leave-out-out cross-validation 
QdaCV.dbc2 <- qda(DBC2, trg.val, CV=T)
QdaCV.dbc2$class
trg.val
table(class.true=trg.val, class.assignedCV=QdaCV.dbc2$class)

errorsqCV <- (QdaCV.dbc2$class != trg.val)
errorsqCV

AERqCV   <- sum(errorsqCV)/length(trg.val)
AERqCV


# Plot the partition induced by QDA
x11()
plot(DBC2, main='VAR1 e VAR2', xlab='VAr1', ylab='VAR2', pch=20)
points(DBC2[trg1,], col='red', pch=20)
points(DBC2[trg2,], col='green', pch=20)
legend("topright", legend=levels(trg.val), fill=c('red','green'))
points(qda.dbc2$means, col=c('red','green'), pch=4, lwd=2, cex=1.5)
x  <- seq(min(DBC2[,1]), max(DBC2[,1]), length=200)
y  <- seq(min(DBC2[,2]), max(DBC2[,2]), length=200)
xy <- expand.grid(Diametro_in_mm =x, CEA=y)   ###############modificare#######
z  <- predict(qda.dbc2, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]     
contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)


### knn-classifier
###----------------
# Plot the partition induced by knn

k <- 20

x11()
plot(DBC2, main='PLOT VAR1 e VAR2', xlab='VA1', ylab='VAR2', pch=20)
points(DBC2[trg1,], col=2, pch=20)
points(DBC2[trg2,], col=3, pch=20)
legend("topright", legend=levels(trg.val), fill=c(2,3,4))

x  <- seq(min(DBC2[,1]), max(DBC2[,1]), length=200)
y  <- seq(min(DBC2[,2]), max(DBC2[,2]), length=200)
xy <- expand.grid((Interval_CT_Surgery_days =x, Age=y)      #################modificre
              
dbc2.knn <- knn(train = DBC2, test = xy, cl = DBC$TRG, k = k)
                  
z  <- as.numeric(dbc2.knn)
                  
contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)
                  
graphics.off()



###----------------------------------------------------------------------------
### CLINICHE - TRG (1,3,5)

##  --> Ho fatto un po di prove ma mi sembra tutto troppo confuso. 
##      Meglio TRG(0,1) alla sezione precedente

#importo dataset
Database_TRG_anonymized_06_03_22 <- read_excel("Database TRG anonymized 06_03_22.xlsx")
load("C:/Users/sarar/Desktop/POLI/APPLIED STATISTIC/Labs/Lab 5/mcshapiro.test.RData")

#solo variabili cliniche fino al TRG (incluso)
varcli <- subset(Database_TRG_anonymized_06_03_22[,1:28])
DB_cliniche <- varcli[complete.cases(varcli[,1]),]


#creo un database che contenga tutte le variabili cliniche numeriche e il TRG
# e elimino l'outlier
DBC<-DB_cliniche[-45,-c(1,2,3,4,6,8,9,10,12,15:23,25:27)]
head(DBC)
dim(DBC)

attach(DBC)

trg.val <- factor(TRG, labels=c('1','3','5'))
trg.val
g=3 

trg1 <- which(trg.val=='1')
trg2 <- which(trg.val=='3')
trg3 <- which(trg.val=='5')

n1 <- length(trg1)
n2 <- length(trg2)
n3 <- length(trg3)
n <- n1+n2+n3

detach(DBC)


#SCELGO SOLO DUE VARIABILI PER LA CLASSIFICAZIONE

# plot all pairs
x11()
plot(DBC[,c(1,2,5)], col=trg.val)

# tengo: 1 intervallo e  2 età
DBC2 <- DBC[,c(1,2)] #12 
head(DBC2)

# plot the data
x11()
plot(DBC2, main='PLOT VAR1 E VAR2', xlab='VAR1', ylab='VAR2', pch=19)
points(DBC2[trg1,], col='red', pch=19)
points(DBC2[trg2,], col='green', pch=19)
points(DBC2[trg3,], col='blue', pch=19)
legend("topright", legend=levels(trg.val), fill=c('red','green','blue'))

m <-  colMeans(DBC2)
m1 <- colMeans(DBC2[trg1,])
m2 <- colMeans(DBC2[trg2,])
m3 <- colMeans(DBC2[trg3,])

S1 <- cov(DBC2[trg1,])
S2 <- cov(DBC2[trg2,])
S3 <- cov(DBC2[trg3,])
Sp  <- ((n1-1)*S1+(n2-1)*S2+(n3-1)*S3)/(n-g)

# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) sigma.A=sigma.B (SOLO PER LDA)
# 3) c(A|B)=c(B|A) (equal misclassification costs)

# verify assumptions:
# 1) normality within the groups 
mcshapiro.test(DBC2[trg1,])$pvalue
mcshapiro.test(DBC2[trg2,])$pvalue
mcshapiro.test(DBC2[trg3,])$pvalue

lambda.mult <- powerTransform(DBC2[,])    
lambda.mult

tr.1 <- bcPower(DBC2[,1],  lambda.mult$lambda[1]) 
tr.2 <- bcPower(DBC2[,2],  lambda.mult$lambda[2]) 
DBC2[,1]<-tr.1
DBC2[,2]<-tr.2
mcshapiro.test(DBC2[trg1,])
mcshapiro.test(DBC2[trg2,])
mcshapiro.test(DBC2[trg3,])



### Quadratic Discriminand Analysis (QDA)
###---------------------------------------

qda.dbc2 <- qda(DBC2, trg.val)
qda.dbc2

Qda.dbc2 <- predict(qda.dbc2, DBC2)

# compute the APER
Qda.dbc2$class
trg.val
table(class.true=trg.val, class.assigned=Qda.dbc2$class)

errorsq <- (Qda.dbc2$class != trg.val)
errorsq

APERq   <- sum(errorsq)/length(trg.val)
APERq


# Compute the estimate of the AER by leave-out-out cross-validation 
QdaCV.dbc2 <- qda(DBC2, trg.val, CV=T)
QdaCV.dbc2$class
trg.val
table(class.true=trg.val, class.assignedCV=QdaCV.dbc2$class)

errorsqCV <- (QdaCV.dbc2$class != trg.val)
errorsqCV

AERqCV   <- sum(errorsqCV)/length(trg.val)
AERqCV


# Plot the partition induced by QDA
x11()
plot(DBC2, main='VAR1 e VAR2', xlab='VAr1', ylab='VAR2', pch=20)
points(DBC2[trg1,], col='red', pch=20)
points(DBC2[trg2,], col='green', pch=20)
points(DBC2[trg3,], col='blue', pch=20)
legend("topright", legend=levels(trg.val), fill=c('red','green','blue'))

points(qda.dbc2$means, col=c('red','green','blue'), pch=4, lwd=2, cex=1.5)

x  <- seq(min(DBC2[,1]), max(DBC2[,1]), length=200)
y  <- seq(min(DBC2[,2]), max(DBC2[,2]), length=200)
xy <- expand.grid(Interval_CT_Surgery_days=x, Age=y)   ###############modificare

z  <- predict(qda.dbc2, xy)$post  
z1 <- z[,1] - pmax(z[,2], z[,3])    
z2 <- z[,2] - pmax(z[,1], z[,3])    
z3 <- z[,3] - pmax(z[,1], z[,2])

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z3, 200), levels=0, drawlabels=F, add=T)


### knn-classifier
###----------------
# Plot the partition induced by knn

k <- 7

x11()
plot(DBC2, main='PLOT VAR1 e VAR2', xlab='VA1', ylab='VAR2', pch=20)
points(DBC2[trg1,], col=2, pch=20)
points(DBC2[trg2,], col=3, pch=20)
points(DBC2[trg3,], col=4, pch=20)
legend("topright", legend=levels(trg.val), fill=c(2,3,4))

x  <- seq(min(DBC2[,1]), max(DBC2[,1]), length=200)
y  <- seq(min(DBC2[,2]), max(DBC2[,2]), length=200)
xy <- expand.grid(Interval_CT_Surgery_days=x, Age=y)      #################modificare

dbc2.knn <- knn(train = DBC2, test = xy, cl = DBC$TRG, k = k)

z  <- as.numeric(dbc2.knn)

contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)

graphics.off()

                  
##---------------------------------------------------------------------
# vARIABILI RADIOLOGICHE SIGNIFICATIVE - RISPOSTA RADIOLOGICA
                  
#importo dataset
Database_TRG_anonymized_06_03_22 <- read_excel("Database TRG anonymized 06_03_22.xlsx")
load("C:/Users/sarar/Desktop/POLI/APPLIED STATISTIC/Labs/Lab 5/mcshapiro.test.RData")
                  
#come per le clinche considero solamente la prima lesione (più grande)
DB<-subset(Database_TRG_anonymized_06_03_22[,c(1,27,29:76)])
DB_rad <-na.omit(DB)
DB_rad <- DB_rad[-45,] # CEA 4216 elimino
DB_rad<-DB_rad[,-1]    # elimino ID pazienti
                  
                  
#creo un database che contenga tutte le variabili radiomiche 
# e la RISPOSTA RADIOLOGICA
                  
DBR<-DB_rad[,-c(5,6,7,9,11,21,30,31,33,34,44,45,47)]
                  
head(DBR)
dim(DBR)
                  
attach(DBR)
                  
rr.val <- factor(Risposta_radiologica, labels=c('PR','SD'))
rr.val
g=2 
                
rr1 <- which(rr.val== 'PR')
rr2 <- which(rr.val=='SD')
rr1
rr2
                  
n1 <- length(rr1)
n2 <- length(rr2)
n <- n1+n2
                  
detach(DBR)
                  
                  
#SCELGO SOLO DUE VARIABILI PER LA CLASSIFICAZIONE
                  
# plot the data
x11()
plot(DBR[,10:15], col=rr.val)
                  
# tengo: 
DBR2 <- DBR[,c(3,8)]  
head(DBR2)
                  
# plot the data
x11()
plot(DBR2, main='VAR1 e VAR2', xlab='VAr1', ylab='VAR2', pch=20)
points(DBR2[rr1,], col='red', pch=20)
points(DBR2[rr2,], col='green', pch=20)
legend("topright", legend=levels(rr.val), fill=c('red','green'))                  
                  
m <-  colMeans(DBR2)
m1 <- colMeans(DBR2[rr1,])
m2 <- colMeans(DBR2[rr2,])                  
                  
S1 <- cov(DBR2[rr1,])
S2 <- cov(DBR2[rr2,])
Sp  <- ((n1-1)*S1+(n2-1)*S2)/(n-g)

# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) sigma.A=sigma.B (SOLO PER LDA)
# 3) c(A|B)=c(B|A) (equal misclassification costs)

# verify assumptions:
# 1) normality within the groups 
mcshapiro.test(DBR2[rr1,])
mcshapiro.test(DBR2[rr2,])                  

### Quadratic Discriminand Analysis (QDA)
###---------------------------------------

qda.dbr2 <- qda(DBR2, rr.val)
qda.dbr2

Qda.dbr2 <- predict(qda.dbr2, DBR2)

# compute the APER
Qda.dbr2$class
rr.val
table(class.true=rr.val, class.assigned=Qda.dbr2$class)

errorsq <- (Qda.dbr2$class != rr.val)
errorsq

APERq   <- sum(errorsq)/length(rr.val)
APERq


# Compute the estimate of the AER by leave-out-out cross-validation 
QdaCV.dbr2 <- qda(DBR2, rr.val, CV=T)
QdaCV.dbr2$class
rr.val
table(class.true=rr.val, class.assignedCV=QdaCV.dbr2$class)

errorsqCV <- (QdaCV.dbr2$class != rr.val)
errorsqCV

AERqCV   <- sum(errorsqCV)/length(rr.val)
AERqCV


# Plot the partition induced by QDA
x11()
plot(DBR2, main='VAR1 e VAR2', xlab='VAr1', ylab='VAR2', pch=20)
points(DBR2[rr1,], col='red', pch=20)
points(DBR2[rr2,], col='green', pch=20)
legend("topright", legend=levels(rr.val), fill=c('red','green'))

points(qda.dbr2$means, col=c('red','green'), pch=4, lwd=2, cex=1.5)

x  <- seq(min(DBR2[,1]), max(DBR2[,1]), length=200)
y  <- seq(min(DBR2[,2]), max(DBR2[,2]), length=200)
xy <- expand.grid(CONVENTIONAL_HUmean=x, HISTO_Entropy_log2=y)   ###############modificare

z  <- predict(qda.dbr2, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]     

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)


open3d()
points3d(DBR2[rr1,1], DBR2[rr1,2], 0, col='red', pch=15)
points3d(DBR2[rr2,1], DBR2[rr2,2], 0, col='green', pch=15)
surface3d(x,y,matrix(dmvnorm(xy, m1, S1) / 3, 50), alpha=0.4, color='red')
surface3d(x,y,matrix(dmvnorm(xy, m2, S2) / 3, 50), alpha=0.4, color='green', add=T)
box3d()
##?????

### knn-classifier
###----------------

k <- 7

x11()
plot(DBR2, main='PLOT VAR1 e VAR2', xlab='VA1', ylab='VAR2', pch=20)
points(DBR2[rr1,], col=2, pch=20)
points(DBR2[rr2,], col=3, pch=20)
legend("topright", legend=levels(rr.val), fill=c(2,3,4))

x  <- seq(min(DBR2[,1]), max(DBR2[,1]), length=200)
y  <- seq(min(DBR2[,2]), max(DBR2[,2]), length=200)
xy <- expand.grid(CONVENTIONAL_HUmean=x, HISTO_Entropy_log2=y)     #################modificare


dbr2.knn <- knn(train = DBR2, test = xy, cl = DBR$Risposta_radiologica, k = k)

z  <- as.numeric(dbr2.knn)

contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)

graphics.off()



##---------------------------------------------------------------------
## VARIABILI RADIOLOGICHE SIGNIFICATIVE CON TRG

#importo dataset
Database_TRG_anonymized_06_03_22 <- read_excel("Database TRG anonymized 06_03_22.xlsx")
load("C:/Users/sarar/Desktop/POLI/APPLIED STATISTIC/Labs/Lab 5/mcshapiro.test.RData")

#come per le clinche considero solamente la prima lesione (più grande)
DB<-subset(Database_TRG_anonymized_06_03_22[,c(1,28,29:76)])
DB_rad <-na.omit(DB)
DB_rad <- DB_rad[-45,] # CEA 4216 elimino
DB_rad<-DB_rad[,-1]    # elimino ID pazienti

#creo un database che contenga tutte le variabili radiomiche e la  TRG

DBR<-DB_rad[,-c(6,7,9,11,12,13,19,20,21,30,31,33,34,45,47)]
head(DBR)
dim(DBR)


#modifico il dataset in modo da avere TRG= 1 o 0
DBR[DBR[,1]<=3, 1]<-0
DBR[DBR[,1]>=5, 1]<-1

attach(DBR)

trg.val <- factor(TRG, labels=c('0','1'))
trg.val
g=2 

trg1 <- which(trg.val== '0')
trg2 <- which(trg.val=='1')
trg1
trg2

n1 <- length(trg1)
n2 <- length(trg2)
n <- n1+n2

detach(DBR)

#SCELGO SOLO DUE VARIABILI PER LA CLASSIFICAZIONE

# plot the data
x11()
plot(DBR[,c(5,31)], col=trg.val)

# tengo: 1 intervallo e  2 età
DBR2 <- DBR[,c(5,31)] #12 
head(DBR2)

# plot the data
x11()
plot(DBR2, main='VAR1 e VAR2', xlab='VAr1', ylab='VAR2', pch=20)
points(DBR2[trg1,], col='red', pch=20)
points(DBR2[trg2,], col='green', pch=20)
legend("topright", legend=levels(trg.val), fill=c('red','green'))


m <-  colMeans(DBR2)
m1 <- colMeans(DBR2[trg1,])
m2 <- colMeans(DBR2[trg2,])


S1 <- cov(DBR2[trg1,])
S2 <- cov(DBR2[trg2,])
Sp  <- ((n1-1)*S1+(n2-1)*S2)/(n-g)

# Assumptions:
# 1) if L=i, X.i ~ N(mu.i, sigma.i^2), i=A,B
# 2) sigma.A=sigma.B (SOLO PER LDA)
# 3) c(A|B)=c(B|A) (equal misclassification costs)

# verify assumptions:
# 1) normality within the groups 
mcshapiro.test(DBR2[trg1,])
mcshapiro.test(DBR2[trg2,])

### Quadratic Discriminand Analysis (QDA)
###---------------------------------------

qda.dbr2 <- qda(DBR2, trg.val)
qda.dbr2

Qda.dbr2 <- predict(qda.dbr2, DBR2)

# compute the APER
Qda.dbr2$class
rr.val
table(class.true=rr.val, class.assigned=Qda.dbr2$class)

errorsq <- (Qda.dbr2$class != rr.val)
errorsq

APERq   <- sum(errorsq)/length(rr.val)
APERq


# Compute the estimate of the AER by leave-out-out cross-validation 
QdaCV.dbr2 <- qda(DBR2, rr.val, CV=T)
QdaCV.dbr2$class
rr.val
table(class.true=rr.val, class.assignedCV=QdaCV.dbr2$class)

errorsqCV <- (QdaCV.dbr2$class != rr.val)
errorsqCV

AERqCV   <- sum(errorsqCV)/length(rr.val)
AERqCV


# Plot the partition induced by QDA
x11()
plot(DBR2, main='VAR1 e VAR2', xlab='VAr1', ylab='VAR2', pch=20)
points(DBR2[trg1,], col='red', pch=20)
points(DBR2[trg2,], col='green', pch=20)
legend("topright", legend=levels(trg.val), fill=c('red','green'))

points(qda.dbr2$means, col=c('red','green'), pch=4, lwd=2, cex=1.5)

x  <- seq(min(DBR2[,1]), max(DBR2[,1]), length=200)
y  <- seq(min(DBR2[,2]), max(DBR2[,2]), length=200)
xy <- expand.grid(CONVENTIONAL_HUmean=x, HISTO_Entropy_log2=y)   ###############modificare

z  <- predict(qda.dbr2, xy)$post  
z1 <- z[,1] - z[,2] 
z2 <- z[,2] - z[,1]     

contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)
contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)


open3d()
points3d(DBR2[trg1,1], DBR2[trg1,2], 0, col='red', pch=15)
points3d(DBR2[trg2,1], DBR2[trg2,2], 0, col='green', pch=15)
surface3d(x,y,matrix(dmvnorm(xy, m1, S1) / 3, 50), alpha=0.4, color='red')
surface3d(x,y,matrix(dmvnorm(xy, m2, S2) / 3, 50), alpha=0.4, color='green', add=T)
box3d()
##?????

### knn-classifier
###----------------

k <- 3

x11()
plot(DBR2, main='PLOT VAR1 e VAR2', xlab='VA1', ylab='VAR2', pch=20)
points(DBR2[trg1,], col=2, pch=20)
points(DBR2[trg2,], col=3, pch=20)
legend("topright", legend=levels(trg.val), fill=c(2,3,4))

x  <- seq(min(DBR2[,1]), max(DBR2[,1]), length=200)
y  <- seq(min(DBR2[,2]), max(DBR2[,2]), length=200)
xy <- expand.grid(CONVENTIONAL_HUmax =x, GLZLM_SZHGE=y)     #################modificare


dbr2.knn <- knn(train = DBR2, test = xy, cl = DBR$TRG, k = k)

z  <- as.numeric(dbr2.knn)

contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)

graphics.off()

                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  