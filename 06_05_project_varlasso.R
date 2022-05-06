##TRG 01 
Database_TRG_anonymized_06_03_22[Database_TRG_anonymized_06_03_22[,28]<=3,28]<-1
Database_TRG_anonymized_06_03_22[Database_TRG_anonymized_06_03_22[,28]>=5,28]<-0
DB_lesmag_age<-subset(Database_TRG_anonymized_06_03_22[,c(1,28,32,34,35,36,37,41,43,45,46,48,51,62,64,65,71,75)])
DB_rad_age <-na.omit(DB_lesmag_age)

species.name <- DB_rad_age[,2]
var.rad        <- DB_rad_age[,3:18]

##Matrice con distanza euclidea
var.rade <- dist(var.rad, method='euclidean')

x11()
image(1:127,1:127,as.matrix(var.rade), main='metrics: Euclidean', asp=1, xlab='i', ylab='j')

var.radm <- dist(var.rad, method='manhattan')
var.radc <- dist(var.rad, method='canberra')


x11()
par(mfrow=c(1,3))
image(1:127,1:127,as.matrix(var.rade), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )
image(1:127,1:127,as.matrix(var.radc), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )
image(1:127,1:127,as.matrix(var.radm), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

###abbiamo scelto la distanza di canberra-> centrino più carino
var.radms <- hclust(var.radc, method='single')
var.radma <- hclust(var.radc, method='average')
var.radmc <- hclust(var.radc, method='complete')

##disegno il dendogram
x11()
par(mfrow=c(1,3))
plot(var.radms, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(var.radmc, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(var.radma, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
###scelgo di suddividere in 2 cluster-> con complete linkage
x11()
par(mfrow=c(1,3))
plot(var.radms, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(var.radms, k=2)
plot(var.radmc, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(var.radmc, k=2)
plot(var.radma, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(var.radma, k=2)

# Fix k=2 clusters:
cluster.ec <- cutree(var.radmc, k=2) # euclidean-complete:
cluster.ec

table(label.true = t(species.name), label.cluster = cluster.ec)

##viene il 55% -> si può solo fare di meglio

# ward linkage

clustw <- hclust(var.radc, method='ward.D2')
x11()
plot(clustw, hang=-0.1, labels=FALSE, main='ward', xlab='', sub='')
rect.hclust(clustw, k=2)
clusterw <- cutree(clustw, 2)
table(label.true = t(species.name), label.cluster = clusterw)
##ancora peggio 53%


##K_MEANS

result.k <- kmeans(var.rad, centers=2) # Centers: fixed number of clusters

names(result.k)

x11()
plot(var.rad[,c(1,15)], col = result.k$cluster+1)



# RIPETO LE ANALISI CON TRG SULLE VARIABILI RADIOLOGICHE SIGNIFICATIVE

attach(DB_rad_age)

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

detach(DB_rad_age)

#SCELGO SOLO DUE VARIABILI PER LA CLASSIFICAZIONE

# plot the data
x11()
plot(DB_rad_age[,c(5,10)], col=trg.val)

# tengo: 1 intervallo e  2 età
DBR2 <- DB_rad_age[,c(5,10)] #12 
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
library(MASS)
qda.dbr2 <- qda(DBR2, trg.val)
qda.dbr2

Qda.dbr2 <- predict(qda.dbr2, DBR2)

# compute the APER
Qda.dbr2$class

table(class.true=trg.val, class.assigned=Qda.dbr2$class)

errorsq <- (Qda.dbr2$class != trg.val)
errorsq

APERq   <- sum(errorsq)/length(trg.val)
APERq


# Compute the estimate of the AER by leave-out-out cross-validation 
QdaCV.dbr2 <- qda(DBR2, trg.val, CV=T)
QdaCV.dbr2$class

table(class.true=trg.val, class.assignedCV=QdaCV.dbr2$class)

errorsqCV <- (QdaCV.dbr2$class != trg.val)
errorsqCV

AERqCV   <- sum(errorsqCV)/length(trg.val)
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



### knn-classifier
###----------------
library(class)
k <- 3

x11()
plot(DBR2, main='PLOT VAR1 e VAR2', xlab='VAR1', ylab='VAR2', pch=20)
points(DBR2[trg1,], col=2, pch=20)
points(DBR2[trg2,], col=3, pch=20)
legend("topright", legend=levels(trg.val), fill=c(2,3,4))

x  <- seq(min(DBR2[,1]), max(DBR2[,1]), length=200)
y  <- seq(min(DBR2[,2]), max(DBR2[,2]), length=200)
xy <- expand.grid(CONVENTIONAL_HUmax =x, GLZLM_SZHGE=y)     #################modificare


dbr2.knn <- knn(train = DBR2, test = xy, cl = DB_rad_age$TRG, k = k)

z  <- as.numeric(dbr2.knn)

contour(x, y, matrix(z, 200), levels=c(1.5, 2.5), drawlabels=F, add=T)

graphics.off()