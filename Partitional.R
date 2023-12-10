library(fpc)
library(cluster)
library(ggplot2)
library(factoextra)

data <- read.csv("data.csv")
x <- data[, -c(1:9)]
distanceMatrix <- dist(x, method = "binary")
xMatrix <- as.matrix(x)

# K-medoids
pam.res <- pam(distanceMatrix, 3)
print(pam.res)
pam.res$clustering <- as.factor(pam.res$clustering)

pca <- PCA(x)
pcaScores <- data.frame(pca$scores)
ggplot(data = pcaScores, aes(x= PC.1, y = PC.2, col = pam.res$clustering)) +
  geom_point(size = 3)

fviz_nbclust(x, hcut, diss = distanceMatrix, k.max = 6, method = "silhouette", hc_method = "ward.D2")
fviz_partitional_ch(x)


# De acuerdo con Calinski y Silhouette, tomamos 3 clusters
# K-medoids
clusters <- 3
pam.res <- pam(distanceMatrix, clusters)
print(pam.res)
pam.res$clustering <- as.factor(pam.res$clustering)

# Evaluación de su estabilidad
clus.boot <- clusterboot(distanceMatrix,
                         distances = TRUE,
                         B=1000, # Number of bootstrap resamples
                         clustermethod=claraCBI, # for hierarchical clustering 
                         k=clusters, 
                         count=T) # Show progress on screen?
# clus.boot
AvgJaccard <- clus.boot$bootmean
Instability <- clus.boot$bootbrd/1000
Clusters <- c(1:clusters)
Eval <- cbind(Clusters, AvgJaccard, Instability)
Eval
#Very stable!

pca <- PCA(x)
pcaScores <- data.frame(pca$scores)
ggplot(data = pcaScores, aes(x= PC.1, y = PC.2, col = pam.res$clustering)) +
  geom_point(size = 3)


fviz_partitional_ch <- function(data) {
  ch <- c()
  distanceMatrix <- dist(data, method = "binary")
  for (i in 2:6) {
    pam.res <- pam(distanceMatrix, i)
    ch[i] <- calinhara(data, # data
                       pam.res$clustering, # cluster assignments
                       cn=max(pam.res$clustering) # total cluster number
    )
  }
  ch <-ch[2:6]
  k <- 2:6
  
  plot(k, ch, 
       lty = 1, type = "o", lwd = 1, pch = 4,
       ylab = "",
       xlab = "Number of clusters",
       main = "K-medoids", cex.main = 0.8,
       col = "dodgerblue1", cex = 0.9 ,
       bty = "l",
       las= 1, cex.axis = 0.8, tcl = -0.2)
  abline(v=which(ch==max(ch)) + 1, lwd=1, col="red", lty="dashed")
  print(max(ch))
}
