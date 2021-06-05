require(data.table)
require(dplyr)
require(plyr)
library(ggplot2)
require(cluster) # sillhuette score
library(igraph)
library(psych) # cor2dist



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

# remove info columns
cb.gdc.val <- cb.gdc[,1:(ncol(cb.gdc)-3)]
cb.kal.val <- cb.kal[,1:(ncol(cb.kal)-3)]

# filter by p-val
cb.gdc.val.sig <- cb.gdc.val[cb.gdc$P.value <= 0.05, ]
cb.kal.val.sig <- cb.kal.val[cb.kal$P.value <= 0.05, ]

# filter cell type by variability
#   want cell types that a
#tmp <- apply(cb.gdc.val.sig.joint, 2, function(x) {sum(x)} )
#tmp[tmp < 1.0]

# filter for prevelance (if too low %, then remove)
# filter for variance (if too little variance, then remove)
#   - want to find/cluster on what is diff... so just remove those that aren't different????



# intersecting samples
overlap.samples <- intersect(rownames(cb.kal.val.sig), rownames(cb.gdc.val.sig))
print(paste(length(overlap.samples),  '/', nrow(cb.kal.val.sig)))

cb.gdc.val.sig.joint <- cb.gdc.val.sig[rownames(cb.gdc.val.sig) %in% overlap.samples, ]
cb.kal.val.sig.joint <- cb.kal.val.sig[rownames(cb.kal.val.sig) %in% overlap.samples, ]

# transpose and re-name columns
cb.gdc.val.sig.joint.t <- transpose(cb.gdc.val.sig.joint)
rownames(cb.gdc.val.sig.joint.t) <- colnames(cb.gdc.val.sig.joint)
colnames(cb.gdc.val.sig.joint.t) <- rownames(cb.gdc.val.sig.joint)
cb.kal.val.sig.joint.t <- transpose(cb.kal.val.sig.joint)
rownames(cb.kal.val.sig.joint.t) <- colnames(cb.kal.val.sig.joint)
colnames(cb.kal.val.sig.joint.t) <- rownames(cb.kal.val.sig.joint)

# sort samples by sample name
cb.gdc.val.sig.joint.t <- cb.gdc.val.sig.joint.t[,order(colnames(cb.gdc.val.sig.joint.t))]
cb.kal.val.sig.joint.t <- cb.kal.val.sig.joint.t[,order(colnames(cb.kal.val.sig.joint.t))]
colnames(cb.kal.val.sig.joint.t) == colnames(cb.gdc.val.sig.joint.t)


# cluster data
# https://psych-networks.com/r-tutorial-identify-communities-items-networks/
data <-      cb.gdc.val.sig.joint.t                                             ## dataset
data.un.t <- cb.gdc.val.sig.joint                                               ## dataset (untransposed)

data <-      cb.kal.val.sig.joint.t                                             ## dataset
data.un.t <- cb.kal.val.sig.joint                                               ## dataset (untransposed)


data.un.t <- data.un.t[order(rownames(data.un.t)), ]

data <- data[ , which(apply(data, 2, var) != 0)]
correlationmatrix = cor(data)
distancematrix <- cor2dist(correlationmatrix)
DM1 <- as.matrix(distancematrix)

# test correlations
#View(apply(correlationmatrix, 1, function(x) {sum(x)/ length(x)}))



## Zero out connections where there is low (absolute) correlation
## Keeps connection for cor ~ -1
## You may wish to choose a different threshhold
DM1[abs(correlationmatrix) < 0.33] = 0  # threshold 0.33
#DM1[correlationmatrix < 0.33] = 0      # no abs (filters out negative correlation)

G1 <- graph.adjacency(DM1, mode = "undirected", weighted = TRUE, diag = TRUE)
vcount(G1)
ecount(G1)

clusterlouvain <- cluster_louvain(G1)
plot(G1, vertex.color=rainbow(12, alpha=0.6)[clusterlouvain$membership])

num.clusters <- length(unique(clusterlouvain$membership))
num.clusters

# sillhuette score
clusterlouvain$names == rownames(data.un.t)
ss <- silhouette(clusterlouvain$membership, dist(data.un.t))
mean(ss[, 3])

# cb.gdc.val.sig.joint
#   ss: -0.05452075
#    k: 4





# testing order
#data.frame(clusterlouvain$membership, rownames(cb.kal.val.sig.joint))

#clusterlouvain$names









# explore LM22 matrix
lm22 <- read.table(
  file = "./bin/cibersort/LM_not_22.txt",
  sep = ",",
  header = TRUE
)
















