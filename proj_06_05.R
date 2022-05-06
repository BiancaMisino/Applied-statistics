
##TRG 01 
Database_TRG_anonymized_06_03_22[Database_TRG_anonymized_06_03_22[,28]<=3,28]<-1
Database_TRG_anonymized_06_03_22[Database_TRG_anonymized_06_03_22[,28]>=5,28]<-0
DB_lesmag_age<-subset(Database_TRG_anonymized_06_03_22[,c(1,28,29:76)])
DB_rad_age <-na.omit(DB_lesmag_age)

species.name <- DB_rad_age[,2]
var.rad        <- DB_rad_age[,3:50]

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

result.k <- kmeans(var.rad, centers=3) # Centers: fixed number of clusters

names(result.k)

x11()
plot(var.rad[,c(6,8,9,10,11)], col = result.k$cluster+1)
