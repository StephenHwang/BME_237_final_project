# Kallisto sleuth on clusters

require(data.table)
require(dplyr)
require(plyr)
library(ggplot2)
require(cluster) # sillhuette score
library(igraph)
library(psych) # cor2dist

require(sleuth)

# gene name to ensemble name
# load data

home.dir <- '/home/stephen/Documents/classes/bme/237/final_project'
setwd(home.dir)

# load expression and meta data
exp.gdc = readRDS(paste0(home.dir, "/data/exp_gdc_recoded_aggr.rds"))
exp.kal = readRDS(paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))
meta.data = readRDS(paste0(home.dir, "/data/meta_data.rds"))

# annotate meta.data by cluster id
gdc_km_4 = readRDS(paste0(home.dir, "/analysis/gdc_km_4.rds"))
kal_km_4 = readRDS(paste0(home.dir, "/analysis/kal_km_4.rds"))

# sort expr by sample name
exp.gdc <- exp.gdc[,order(colnames(exp.gdc))]
exp.kal <- exp.kal[,order(colnames(exp.kal))]

gdc_km_4 <- gdc_km_4$cluster[order(names(gdc_km_4$cluster))]
kal_km_4 <- kal_km_4$cluster[order(names(kal_km_4$cluster))]

# filter metadata and expression data
meta.data <-  meta.data[meta.data$sample %in% intersect(meta.data$sample, names(gdc_km_4)), ]
meta.data <- meta.data[order(meta.data$sample), ]
meta.data$gdc_clusters <- gdc_km_4
meta.data$kal_clusters <- kal_km_4

# filter sort expression data
exp.kal <- exp.kal[, colnames(exp.kal) %in% intersect(colnames(exp.kal), names(kal_km_4))]
exp.kal <- exp.kal[, order(colnames(exp.kal))]




# load in data
exp.gdc = readRDS(paste0(home.dir, "/data/exp_gdc_recoded_aggr.rds"))
exp.kal = readRDS(paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))

# sort data
exp.kal <- exp.kal[, order(colnames(exp.kal))]
exp.gdc <- exp.gdc[, order(colnames(exp.gdc))]








# load expression and metadata by cluster
kal.cluster.1 <- read.table(
  file = "./analysis/kal_cluster_1.tsv",
  sep = "\t",
  header = TRUE
)

kal.cluster.2 <- read.table(
  file = "./analysis/kal_cluster_2.tsv",
  sep = "\t",
  header = TRUE
)

kal.cluster.3 <- read.table(
  file = "./analysis/kal_cluster_3.tsv",
  sep = "\t",
  header = TRUE
)

kal.cluster.4 <- read.table(
  file = "./analysis/kal_cluster_4.tsv",
  sep = "\t",
  header = TRUE
)

rownames(kal.cluster.1) <- kal.cluster.1$X
rownames(kal.cluster.2) <- kal.cluster.2$X
rownames(kal.cluster.3) <- kal.cluster.3$X
rownames(kal.cluster.4) <- kal.cluster.4$X


# sleuth analysis

so <- sleuth_prep(s2c)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

plot_pca(so, color_by = 'condition')






















###############################################################################











































