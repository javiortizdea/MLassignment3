library(fpc)
library(cluster)
library(mclust)

rawData <- read.csv("data.csv")
T6SS <- rawData[3]
x <- rawData[-c(1:9)]
# Elección de número de clusters
linkages = c("single", "complete", "average", "centroid", "ward.D", "ward.D2")
#png(filename = "AllCalinski.png", width = 600, res = 120)
par(mfrow=c(2,3))
for (linkage in linkages) {
  fviz_hierarchical_ch(x, linkage)
}
mtext("Calinksi-Harabasz index", outer = T, side = 3, line = -2)
#dev.off()

#Comprobamos con Silhouette
fviz_nbclust(x, hcut, diss = distanceMatrix, k.max = 6, method = "silhouette", hc_method = "ward.D2")


# Nos quedamos con 3 clusters, y con el de Ward D2, que parece haber logrado
#un mayor Calinksi
distanceMatrix <- dist(x, method = "binary")
clusters <- 3
clus.boot <- clusterboot(distanceMatrix,
                         distances = TRUE,
                         B=1000, # Number of bootstrap resamples
                         clustermethod=hclustCBI, # for hierarchical clustering 
                         method="ward.D2", # use what we used in "hclust"
                         k=clusters, 
                         count=T) # Show progress on screen?
# clus.boot
AvgJaccard <- clus.boot$bootmean
Instability <- clus.boot$bootbrd/1000
Clusters <- c(1:clusters)
Eval <- cbind(Clusters, AvgJaccard, Instability)
Eval
#Very stable!
#Dendrogram
par(mfrow = c(1,1))
hcWard2 <- hclust(distanceMatrix, method = "ward.D2")
plot(hcWard2, labels = F)
rect.hclust(hcWard2, k=clusters)
#Clusterplot
groups <- cutree(hcWard2, k = clusters)
groups <- as.factor(groups)
clusplot(x, groups, color = T, shade = T, labels = 4, lines = 0)

pca <- PCA(x)
pcaScores <- data.frame(pca$scores)
(pca$var[1] + pca$var[2])/pca$totalvar
ggplot(data = pcaScores, aes(x= PC.1, y = PC.2, col = groups)) +
  geom_point(size = 3)

# Hierarchical clustering
# ++++++++++++++++++++++++
# Use hcut() which compute hclust and cut the tree
hc.cut <- hcut(distanceMatrix, k = 3, hc_method = "ward.D2")
# Visualize dendrogram
fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE)
# Visualize cluster
fviz_cluster(hc.cut, x, show_labels = F)


fviz_hierarchical_ch <- function(data, linkage) {
  ch <- c()
  #linkages <- c("single", "complete", "average", "ward.D", "ward.D2", "centroid", "mcquitty")
  distanceMatrix <- dist(data, method = "binary")
  for (i in 2:6) {
    #km <- kmeans(data, i) # perform clustering
    hcClustering <- cutree(hclust(distanceMatrix, method = linkage), k = i)
    ch[i] <- calinhara(data, # data
                       hcClustering, # cluster assignments
                       cn=max(hcClustering) # total cluster number
    )
  }
  ch <-ch[2:6]
  k <- 2:6
  
  plot(k, ch, 
       lty = 1, type = "o", lwd = 1, pch = 4,
       ylab = "",
       xlab = "Number of clusters",
       main = paste(linkage, "linkage"), cex.main = 0.8,
       col = "dodgerblue1", cex = 0.9 ,
       bty = "l",
       las= 1, cex.axis = 0.8, tcl = -0.2)
  #plot(k, ch,xlab =  "Cluster number k",
  #     ylab = "Caliński - Harabasz Score",
  #     main = paste("Caliński - Harabasz Plot, linkage:", linkage), cex.main=1,
  #     col = "dodgerblue1", cex = 0.9 ,
  #     lty=1 , type="o" , lwd=1, pch=4,
  #     bty = "l",
  #     las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(ch==max(ch)) + 1, lwd=1, col="red", lty="dashed")
  print(max(ch))
  
  return(ch)
}
