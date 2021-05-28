# Clustering on gdc and kal data
setwd("/home/tmpereir/final")
# install.packages("factoextra")
library(factoextra)
library(data.table)


home.dir <- '/home/stephen/Documents/classes/bme/237/final_project'
setwd(home.dir)

# load expression and meta data
exp.gdc = readRDS(paste0(home.dir, "/data/exp_gdc_recoded_aggr.rds"))
exp.kal = readRDS(paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))
meta.data = readRDS(paste0(home.dir, "/data/meta_data.rds"))

# load cibersort results
cb.gdc <- read.table(
  file = "./analysis/cb_gdc_results.txt",
  sep = "\t",
  header = TRUE
)

cb.kal <- read.table(
  file = "./analysis/cb_kal_results.txt",
  sep = "\t",
  header = TRUE
)

cb.gdc.val <- cb.gdc[,1:(ncol(cb.gdc)-3)]
cb.kal.val <- cb.kal[,1:(ncol(cb.kal)-3)]

cb.gdc.val.sig <- cb.gdc.val[cb.gdc$P.value <= 0.05, ]


gdc_data <- cb.gdc.val.sig




# 1. Determine k using fviz_nbclust() from factoextra (Within Sum Squares)
set.seed(123)

fviz_nbclust(gdc_data, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2)

# Output recommends 3 clusters


# 2. Determine k using average silhouette width

library(cluster)
sil <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
   k1to20 <- kmeans(gdc_data, centers = i, nstart = 25, iter.max = 20)
   ss <- silhouette(k1to20$cluster, dist(gdc_data))
   sil[i] <- mean(ss[, 3])
 }

# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)

# Output: Recommends 1 clusters

# 3. Calinsky Criterion

# install.packages("vegan")
#library(vegan)
# fit <- cascadeKM(gdc_data, 1, 20, iter = 100)
# plot(fit, sortg = TRUE, grpmts.plot = TRUE)

# calinski.best <- as.numeric(which.max(fit$results[2,]))
# cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

# Output: Calinski criterion optimal number of clusters: 2

# 4. Gap statistic

# set.seed(13)
# gap <- clusGap(gdc_data, kmeans, 20, B = 100, verbose = interactive())
# plot(gap, main = "Gap statistic")
# abline(v=which.max(gap$Tab[,3]), lty = 2)

# Output: recommends 17

###########################################

# Compute kmeans

# Compute k-means with k = 3
# set.seed(123)

gdc_data <- gdc_data[ , which(apply(gdc_data, 2, var) != 0)]


km_gdc <- kmeans(gdc_data, 3, nstart = 25)
print(km_gdc)
fviz_cluster(km_gdc, data = gdc_data)



#######################################################################
kal_data <- read.table(file = 'cb_kal_results.txt', sep = '\t', header = TRUE)

# Only keep samples with p-value <= 0.05
kal_data <- kal_data[kal_data$P.value <= 0.05,]

# Get rid of p-values, correlation, RMSE
kal_data <- kal_data[,1:(ncol(gdc_data)-3)]

# fviz_nbclust(kal_data, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2)
# Output: 3

# 2. Determine k using average silhouette width

# library(cluster)
# sil <- rep(0, 20)
# repeat k-means for 1:20 and extract silhouette:
# for(i in 2:20){
#   k1to20 <- kmeans(gdc_data, centers = i, nstart = 25, iter.max = 20)
#   ss <- silhouette(k1to20$cluster, dist(gdc_data))
#   sil[i] <- mean(ss[, 3])
# }
# Plot the  average silhouette width
# plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
# abline(v = which.max(sil), lty = 2)
# OUTPUT: 2

# 3. Calinsky Criterion

# install.packages("vegan")
# library(vegan)
# fit <- cascadeKM(gdc_data, 1, 20, iter = 100)
# plot(fit, sortg = TRUE, grpmts.plot = TRUE)

# calinski.best <- as.numeric(which.max(fit$results[2,]))
# cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
# OUTPUT: 2

kal.km <- kmeans(gdc_data, 3, nstart = 25)
print(kal.km)


fviz_cluster(kal.km, data = kal_data)
