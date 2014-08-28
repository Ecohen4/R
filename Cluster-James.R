library(maps)
library(mapdata)
library(mclust)

## Cluster the average Western U.S. winter precipitation. Compute the average winter precipitation for each climate division in the Western U.S – i.e., the data set will be a single column of one average precipitation value per division.

## x is the data matrix to cluster
## columns are the attributes to cluster on (lat-long-precip)
## rows are records/observations to cluster (years)
## "rainw" = winter precip for 81 climdivs in western U.S. (1986-2013)
x<-apply(X=rainw, MARGIN=2, FUN=mean) # 'climatology' for each climdiv

## Part (i) identify the number of clusters, K
* Select a desired number of clusters, say, j; 
* cluster the data and compute the WSS:
* repeat for j=1:10; 
* plot j versus WSS:
* select number of clusters K, where the WSS starts to saturate.
# (ii) cluster the data into K clusters and display them.
# (iii) Repeat (i)-(ii) by including latitude and longitude – the data set will now be a matrix of 3 columns – precipitation, latitude and longitude.

## K-means cluster.
nk = 7  # number of clusters to consider
wss=1:nk
x<-as.data.frame(x)
wss <- (nrow(x)-1)*sum(apply(x,2,var))
for (i in 2:nk) wss[i] <- sum(kmeans(x,centers=i)$withinss)
plot(wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

## select the K corresponding to the the minimum WSS or WSS beyond which the drop is small. Similar to the Eigen spectrum.
kbest = 3

### Use BIC to get best number of clusters..
library(mclust)
d_clust <- Mclust(as.matrix(x), G=1:10)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
# 4 clusters
# par(mfrow=c(2,2), mar=c(2,2,4,1)+0.1)
plot(d_clust)


## Repeat using lat-long-precip instead of precip alone.
x<-data.frame(x=px, y=py, rain=x)

  
  
zclust = kmeans(x,centers=kbest)
plot(lon,lat,col=c(zclust$cluster),xlab="",ylab="")
US(add=TRUE, col="grey", lwd=2,xlim=range(-125,-100))




######################## 
# Compute average winter precipitation for each division in the western US
dt = data.frame(p = colMeans(wpcpy), x = px, y = py)
dt = dt[order(dt$p),]

# Identify number of clusters (K)
clusters = lapply(1:10, function(k) kmeans(dt, k, iter.max = 10000))
wss = sapply(clusters, function(cluster) cluster$tot.withinss)
par(baseparTitle)
plot(wss, type="b", xlab="Number of clusters (K)", ylab="Within-group sum of squares", main = "Cluster WSS for mean winter precipitation (K = 1:10)")
#... using mclustBIC
K = which.min(abs(mclustBIC(dt)[,2]))
abline(v = K, col = 'grey')
points(K, wss[K], pch = 20)
axis(1, at = K)

# Plot clusters and cluster centers
par(baseparTitle)
plot(dt$x, dt$y, type = 'n', xlab = "Longitude", ylab = "Latitude)", main = "Clusters of mean winter precipitation (K = 7)")
detach("package:mclust", unload=TRUE) # mclust breaks US map...
US(add = TRUE, col = "lightgrey", lwd = 2)
for (i in 1:K) {
  plotConvexHull(dt$x[clusters[[K]]$cluster == i], dt$y[clusters[[K]]$cluster == i], lty = 2)
}
points(clusters[[K]]$centers[,2], clusters[[K]]$centers[,3], pch = '+')
text(dt$x, dt$y, clusters[[K]]$cluster)

#######################
## Partition using mediods - similar to K-means but using the cluster median..
## example using ggplot2 - feel free to replace it with your data.

library(cluster)
library(ggplot2)

x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)))

max_k <- 5
sil <- numeric(length(2:max_k))

for(i in 2:max_k) {
  p <- pam(x, i, stand=TRUE)
  sil[i-1] <- mean(silhouette(p))
}
qplot(2:max_k,sil,geom='line')+theme_bw()

# Then once you know the optimal number of clusters:
  
  k <- 2
clusters <- pam(x, k, stand=TRUE, cluster.only=T)
qplot(x[,1],x[,2],color=factor(clusters))+theme_bw()

##############