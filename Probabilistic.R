library(fpc)
library(cluster)
library(ClusterR)
library(mclust)
library(factoextra)
library(ggfortify)


data <- read.csv("data.csv")
x <- data[, -c(1:9)]
distanceMatrix <- dist(x, method = "binary")



fviz_gmm_ch(x)
fviz_gmm_sil(x)

pam.res <- pam(distanceMatrix, 2)
gmm <- Mclust(x, 2, initialModel = pam.res)

xWithLabels <- cbind(x, cluster = as.factor(gmm$classification))
gmm.boot <- MclustBootstrap(gmm, nboot = 1000)


pca_res <- prcomp(x, scale. = FALSE)
autoplot(pca_res, data = xWithLabels, loadings = TRUE, colour = "cluster",
         loadings.label = T,
         loadings.label.size = 2,
         size = 3,
         scale = 0)

pca <- PCA(x)
pcaScores <- data.frame(pca$scores)
ggplot(data = pcaScores, aes(x= PC.1, y = PC.2, col = pam.res$clustering)) +
  geom_point(size = 3)

sil <- silhouette(gmm$classification, dist(x, method = "binary"))
summary(sil)

clus.boot.gmm <- clusterboot(x,
                         distances = F,
                         B=50, # Number of bootstrap resamples
                         clustermethod=gmmCBI, # for GMM (I made it!) 
                         nClusters=2,
                         count=T) # Show progress on screen?

# clus.boot
AvgJaccard <- clus.boot.gmm$bootmean
Instability <- clus.boot.gmm$bootbrd/50
Clusters <- c(1:2)
Eval <- cbind(Clusters, AvgJaccard, Instability)
Eval

sil$mean
fviz_gmm_ch <- function(data) {
  ch <- c()
  #hcWard2 <- hclust(distanceMatrix, method = "ward.D2")
  pam.res <- pam(distanceMatrix, 3)
  for (i in 2:6) {
    gmm <- Mclust(data, i, initialModel = pam.res)
    ch[i] <- calinhara(data, # data
                       gmm$classification, # cluster assignments
                       cn=max(gmm$classification) # total cluster number
    )
  }
  ch <-ch[2:6]
  k <- 2:6
  
  plot(k, ch, 
       lty = 1, type = "o", lwd = 1, pch = 4,
       ylab = "",
       xlab = "Number of clusters",
       main = "Gaussian Mixture Model", cex.main = 0.8,
       col = "dodgerblue1", cex = 0.9 ,
       bty = "l",
       las= 1, cex.axis = 0.8, tcl = -0.2)
  abline(v=which(ch==max(ch)) + 1, lwd=1, col="red", lty="dashed")
  print(max(ch))
}
fviz_gmm_sil <- function(data) {
  ch <- c()
  distanceMatrix <- dist(data, method = "binary")
  #hcWard2 <- hclust(distanceMatrix, method = "ward.D2")
  pam.res <- pam(distanceMatrix, 3)
  for (i in 2:6) {
    gmm <- Mclust(data, i, initialModel = pam.res)
    sil <- silhouette(gmm$classification, dist(data, method = "binary"))
    ch[i] <- mean(sil[,3])
  }
  ch <-ch[2:6]
  k <- 2:6
  
  plot(k, ch, 
       lty = 1, type = "o", lwd = 1, pch = 4,
       ylab = "",
       xlab = "Number of clusters",
       main = "Gaussian Mixture Model", cex.main = 0.8,
       col = "dodgerblue1", cex = 0.9 ,
       bty = "l",
       las= 1, cex.axis = 0.8, tcl = -0.2)
  abline(v=which(ch==max(ch)) + 1, lwd=1, col="red", lty="dashed")
  print(max(ch))
}

gmmCBI <- function(x, nClusters) {
  distanceMatrix <- dist(x, method = "binary")
  pam.res <- pam(distanceMatrix, nClusters)
  gmm <- Mclust(x, nClusters, initialModel = pam.res)
  res <- list("result" = gmm,
              "nc" = nClusters,
              "clusterlist" = list (
                "1" = gmm$classification == 1,
                "2" = gmm$classification == 2
              ),
              "partition" = gmm$classification,
              "clustermethod" = "Gaussian mixture model")
  return(res)
}
