# Clustering on gdc and kal data
setwd("/home/tmpereir/final")
# install.packages("factoextra")
library(factoextra)
library(stringr)
set.seed(123)
# Loading and preprocessing GDC and KAL data
gdc.data <- read.table(file = 'cb_gdc_results.txt', sep = '\t', header = TRUE)
kal.data <- read.table(file = 'cb_kal_results.txt', sep = '\t', header = TRUE)

# Only keep samples with p-value <= 0.05
gdc.data <- gdc.data[gdc.data$P.value <= 0.05,]
kal.data <- kal.data[kal.data$P.value <= 0.05,]

# Get rid of p-values, correlation, RMSE
gdc.data <- gdc.data[,1:(ncol(gdc.data)-3)]
kal.data <- kal.data[,1:(ncol(kal.data)-3)]

# Remove columns with 0 variance
gdc.data <- gdc.data[,which(apply(gdc.data, 2, var) != 0)]
kal.data <- kal.data[,which(apply(kal.data, 2, var) != 0)]

overlap.samples <- intersect(rownames(kal.data), rownames(gdc.data))
cb.gdc.val.sig.joint <- gdc.data[rownames(gdc.data) %in% overlap.samples, ]
cb.kal.val.sig.joint <- kal.data[rownames(kal.data) %in% overlap.samples, ]

# 1. Determine k using fviz_nbclust() from factoextra (Within Sum Squares)
# set.seed(123)
# print("Determining k using wss")
# fviz_nbclust(gdc.data, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2)

# Output recommends 3 clusters


# 2. Determine k using average silhouette width

library(cluster)
sil <- rep(0, 20)

# GDC:
# repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
   k1to20 <- kmeans(cb.gdc.val.sig.joint, centers = i, nstart = 25, iter.max = 20)
   ss <- silhouette(k1to20$cluster, dist(cb.gdc.val.sig.joint))
   sil[i] <- mean(ss[, 3])
}
print(sil)

# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)

# Kallisto:
sil <- rep(0, 20)

for(i in 2:20){
  k1to20 <- kmeans(cb.kal.val.sig.joint, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(cb.kal.val.sig.joint))
  sil[i] <- mean(ss[, 3])
}
print(sil)

# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)

############################################
# OUTPUT 
# figure: silhuette plot 2-20
# - for clusters sizes k 2-5
#   - gdc_2: 0.2801452
#   - gdc_3: 0.2680178
#   - gdc_4: 0.2520860
#   - gdc_5: 0.2219470
# - k recommended: 2
# 
# - kal_2: 0.2420712
# - kal_3: 0.2417190 
# - kal_4: 0.2498931
# - kal_5: 0.2295154
# - k recommended: 4

############################################

# 3. Calinsky Criterion

# install.packages("vegan")
# library(vegan)
# fit <- cascadeKM(gdc.data, 1, 20, iter = 100)
# plot(fit, sortg = TRUE, grpmts.plot = TRUE)

# calinski.best <- as.numeric(which.max(fit$results[2,]))
# cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

# Output: Calinski criterion optimal number of clusters: 3

# 4. Gap statistic

# set.seed(13)
# library(cluster)
# gap <- clusGap(gdc.data, kmeans, 20, B = 100, verbose = interactive())
# plot(gap, main = "Gap statistic")
# abline(v=which.max(gap$Tab[,3]), lty = 2)

# Output: recommends 18

#############################################

# Compute kmeans

# Compute k-means with k = 3
# set.seed(123)
km.2.gdc <- kmeans(cb.gdc.val.sig.joint, 2, nstart = 25)
km.3.gdc <- kmeans(cb.gdc.val.sig.joint, 3, nstart = 25)
km.4.gdc <- kmeans(cb.gdc.val.sig.joint, 4, nstart = 25)
km.5.gdc <- kmeans(cb.gdc.val.sig.joint, 5, nstart = 25)

km.2.kal <- kmeans(cb.kal.val.sig.joint, 2, nstart = 25)
km.3.kal <- kmeans(cb.kal.val.sig.joint, 3, nstart = 25)
km.4.kal <- kmeans(cb.kal.val.sig.joint, 4, nstart = 25)
km.5.kal <- kmeans(cb.kal.val.sig.joint, 5, nstart = 25)

# GDC Kmeans graphs
fviz_cluster(km.2.gdc, data = cb.gdc.val.sig.joint, labelsize=0, main='GDC clusters, k = 2')
fviz_cluster(km.3.gdc, data = cb.gdc.val.sig.joint, labelsize=0, main='GDC clusters, k = 3')
fviz_cluster(km.4.gdc, data = cb.gdc.val.sig.joint, labelsize=0, main='GDC clusters, k = 4')
fviz_cluster(km.5.gdc, data = cb.gdc.val.sig.joint, labelsize=0, main='GDC clusters, k = 5')

# Kallisto Kmeans graphs
fviz_cluster(km.2.kal, data = cb.kal.val.sig.joint, labelsize=0, main='KAL clusters, k = 2')
fviz_cluster(km.3.kal, data = cb.kal.val.sig.joint, labelsize=0, main='KAL clusters, k = 3')
fviz_cluster(km.4.kal, data = cb.kal.val.sig.joint, labelsize=0, main='KAL clusters, k = 4')
fviz_cluster(km.5.kal, data = cb.kal.val.sig.joint, labelsize=0, main='KAL clusters, k = 5')

# Rand index between k=2 to k=5
rand.k.2 <- rand.index(km.2.gdc$cluster[order(names(km.2.gdc$cluster))], km.2.kal$cluster[order(names(km.2.kal$cluster))])  # 0.4956
rand.k.3 <- rand.index(km.3.gdc$cluster[order(names(km.3.gdc$cluster))], km.3.kal$cluster[order(names(km.3.kal$cluster))])  # 0.5436
rand.k.4 <- rand.index(km.4.gdc$cluster[order(names(km.4.gdc$cluster))], km.4.kal$cluster[order(names(km.4.kal$cluster))])  # 0.6052
rand.k.5 <- rand.index(km.5.gdc$cluster[order(names(km.5.gdc$cluster))], km.5.kal$cluster[order(names(km.5.kal$cluster))])  # 0.6513


# Rand index scores
# k=2: 0.4956 
# k=3: 0.5436
# k=4: 0.6042
# k=5: 0.6513

saveRDS(km.4.gdc, file = "gdc_km_4.rds")
saveRDS(km.4.kal, file = "kal_km_4.rds")

# Annotations for cluster and cluster dataset
# km.4.kal$cluster
# > km.4.gdc$cluster[order(names(km.4.gdc$cluster))]
# save 2 .rds for kmeans -> push to analysis

#rand.k.4. <- rand.index(km.4.gdc$cluster[order(names(km.4.gdc$cluster))], km.4.kal$cluster[order(names(km.4.kal$cluster))])
#gdc.rownames <- order(km.4.gdc$cluster)

# gdc_counts <- read.csv("exp_gdc_recoded_aggr.tsv",sep="\t")
# 
# 
# colnames(gdc_counts) <- gsub("[.]", "-", colnames(gdc_counts))
# #for (i in 1:length(colnames(gdc_counts))) {
# #  print(gsub(".", "-", colnames(gdc_counts)[i]))
# #  #colnames(gdc_counts)[i] <- str_replace(coln, ".", "-")
# #}
# 
# cluster.1.gdc.samples <- rownames(which(clusters==1,arr.ind=T))
# cluster.2.gdc.samples <- rownames(which(clusters==2,arr.ind=T))
# cluster.3.gdc.samples <- rownames(which(clusters==3,arr.ind=T))
# 
# #write.csv(gdc_counts[cluster.2.gdc.samples,], "gdc_cluster_2.csv")
# #write.csv(gdc_counts[cluster.3.gdc.samples,], "gdc_cluster_3.csv")
# 
