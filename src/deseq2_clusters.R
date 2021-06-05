# Differential Gene Analysis across clusters from GDC and Kallisto data
setwd("/home/tmpereir/final")

# Packages
library("DESeq2")
library("tidyverse")

# Load GDC clusters
gdc.clusters <- readRDS(file="gdc_km_4.rds")

coldata <- as.data.frame(gdc.clusters$cluster, colnames = "Cluster")
colnames(coldata) <- c("cluster")
coldata <- arrange(coldata, row.names(coldata))
coldata$cluster <- factor(coldata$cluster)

# Get sample names for each cluster

# gdc.cluster.1.samples <- subset(coldata, coldata$cluster == 1)
# gdc.cluster.2.samples <- subset(coldata, coldata$cluster == 2)
# gdc.cluster.3.samples <- subset(coldata, coldata$cluster == 3)
# gdc.cluster.4.samples <- subset(coldata, coldata$cluster == 4)

# Load GDC count
# i-th row and the j-th column of the matrix tells how many reads can be assigned to gene i in sample j.
gdc.counts <- readRDS(file="exp_gdc_recoded_aggr.rds")

# only retain the 105 variables included in the clusters 
gdc.counts <- gdc.counts[,rownames(coldata)]  

### Construct cluster comparisons ###

# Cluster 1 v 2

# gdc.1v2.samples <- subset(coldata, (coldata$cluster == 1 | coldata$cluster == 2))
# gdc.1v2.samples <- arrange(gdc.1v2.samples, row.names(gdc.1v2.samples))
# gdc.1v2.counts <- gdc.counts[,rownames(gdc.1v2.samples)]
# gdc.1v2.dds <- DESeqDataSetFromMatrix(countData = gdc.1v2.counts, colData = gdc.1v2.samples, design = ~ cluster)
# gdc.1v2.dds <- DESeq(gdc.1v2.dds)
# gdc.1v2.results <- results(gdc.1v2.dds, alpha=0.01)  # cutoff 0.01
# gdc.1v2.results.ordered <- gdc.1v2.results[order(gdc.1v2.results$pvalue),]
# > results(gdc.1v2.dds)

# log2 fold change (MLE): cluster 2 vs 1 
# Wald test p-value: cluster 2 vs 1 
# DataFrame with 38481 rows and 6 columns
#             baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#            <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# A1BG        26.10325     -0.5151610  0.299407 -1.720605 0.0853226  0.441247
# A1BG-AS1   117.02513     -0.1507249  0.293670 -0.513246 0.6077793  0.864460
# A1CF         1.40726     -0.2194781  0.758064 -0.289524 0.7721800        NA
# A2M      25388.58834      0.6339904  0.323554  1.959456 0.0500594  0.353206
# A2M-AS1    157.96678      0.0404063  0.285088  0.141733 0.8872912  0.969581
# ...              ...            ...       ...       ...       ...       ...
# ZYG11B       2173.97     -0.1126270  0.154095 -0.730892 0.4648450  0.792240
# ZYX         28817.25      0.1798498  0.244357  0.736012 0.4617236  0.790712
# ZYXP1           0.00             NA        NA        NA        NA        NA
# ZZEF1        2512.65     -0.0332544  0.291037 -0.114262 0.9090304  0.976560
# ZZZ3         2612.96     -0.2911954  0.174468 -1.669045 0.0951084  0.462155

# Number of genes with p-values lower than 0.1
# > sum(gdc.1v2.results$padj < 0.1, na.rm=TRUE)
# [1] 872
################################################################################

# Cluster 1 v 3
# gdc.1v3.samples <- subset(coldata, (coldata$cluster == 1 | coldata$cluster == 3))
# gdc.1v3.samples <- arrange(gdc.1v3.samples, row.names(gdc.1v3.samples))
# gdc.1v3.counts <- gdc.counts[,rownames(gdc.1v3.samples)]
# gdc.1v3.dds <- DESeqDataSetFromMatrix(countData = gdc.1v3.counts, colData = gdc.1v3.samples, design = ~ cluster)
# gdc.1v3.dds <- DESeq(gdc.1v3.dds)

# > results(gdc.1v3.dds)
# log2 fold change (MLE): cluster 3 vs 1 
# Wald test p-value: cluster 3 vs 1 
# DataFrame with 38481 rows and 6 columns
#             baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#            <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# A1BG        29.10141      -0.364429  0.392873 -0.927602  0.353614  0.829559
# A1BG-AS1   118.39794      -0.193437  0.383636 -0.504221  0.614106  0.927775
# A1CF         1.92996       0.462438  0.816576  0.566313  0.571181        NA
# A2M      20510.73659       0.200554  0.287217  0.698267  0.485010  0.893383
# A2M-AS1    213.42932       0.519517  0.386860  1.342907  0.179302  0.696322
# ...              ...            ...       ...       ...       ...       ...
# ZYG11B       2500.11      0.0885034  0.147928  0.598287  0.549648  0.912319
# ZYX         25543.49     -0.1014468  0.244390 -0.415102  0.678068  0.946222
# ZYXP1           0.00             NA        NA        NA        NA        NA
# ZZEF1        2551.86     -0.0653491  0.219697 -0.297451  0.766122  0.966647
# ZZZ3         2721.69     -0.2789569  0.196674 -1.418370  0.156083  0.669724

# Trying LFC shrinkage to reduce noise
# gdc.1v2.results.LFC <- lfcShrink(gdc.1v2.dds, coef="cluster_2_vs_1", type="apeglm")

# > gdc.1v2.results.LFC
# log2 fold change (MAP): cluster 2 vs 1 
# Wald test p-value: cluster 2 vs 1 
# DataFrame with 38481 rows and 5 columns
# baseMean log2FoldChange      lfcSE    pvalue      padj
# <numeric>      <numeric>  <numeric> <numeric> <numeric>
#   A1BG        26.10325   -2.26284e-06 0.00144268 0.0853226  0.441247
# A1BG-AS1   117.02513   -1.98549e-06 0.00144268 0.6077793  0.864460
# A1CF         1.40726   -3.78292e-07 0.00144269 0.7721800        NA
# A2M      25388.58834    5.58707e-06 0.00144269 0.0500594  0.353206
# A2M-AS1    157.96678    5.58140e-07 0.00144268 0.8872912  0.969581
# ...              ...            ...        ...       ...       ...
# ZYG11B       2173.97    4.14112e-05 0.00144293 0.4648450  0.792240
# ZYX         28817.25   -7.28618e-05 0.00144359 0.4617236  0.790712
# ZYXP1           0.00             NA         NA        NA        NA
# ZZEF1        2512.65   -5.55216e-05 0.00144321 0.9090304  0.976560
# ZZZ3         2612.96   -2.51557e-06 0.00144264 0.0951084  0.462155

# > sum(gdc.1v2.results.LFC$padj < 0.05, na.rm=TRUE)
# [1] 528
# > sum(gdc.1v2.results.LFC$padj < 0.01, na.rm=TRUE)
# [1] 225

################################################################################

# Cluster 1 v 4
# gdc.1v4.samples <- subset(coldata, (coldata$cluster == 1 | coldata$cluster == 4))
# gdc.1v4.samples <- arrange(gdc.1v4.samples, row.names(gdc.1v4.samples))
# gdc.1v4.counts <- gdc.counts[,rownames(gdc.1v4.samples)]
# gdc.1v4.dds <- DESeqDataSetFromMatrix(countData = gdc.1v4.counts, colData = gdc.1v4.samples, design = ~ cluster)
# gdc.1v4.dds <- DESeq(gdc.1v4.dds)

# > results(gdc.1v4.dds)
# log2 fold change (MLE): cluster 4 vs 1 
# Wald test p-value: cluster 4 vs 1 
# DataFrame with 38481 rows and 6 columns
#             baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#            <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# A1BG        26.26966      -0.643862  0.324902 -1.981709 0.0475118  0.543098
# A1BG-AS1   112.71402      -0.367239  0.358982 -1.023003 0.3063062  0.844595
# A1CF         1.81303       0.306977  0.778169  0.394486 0.6932220        NA
# A2M      21394.58640       0.191736  0.209791  0.913938 0.3607493  0.866559
# A2M-AS1    197.80561       0.300189  0.337782  0.888707 0.3741608  0.870264
# ...              ...            ...       ...       ...       ...       ...
# ZYG11B       2453.95     -0.0307184  0.171711 -0.178896  0.858020  0.984115
# ZYX         27038.22     -0.0792944  0.240274 -0.330016  0.741388  0.969741
# ZYXP1           0.00             NA        NA        NA        NA        NA
# ZZEF1        2875.03      0.0706517  0.222235  0.317915  0.750550  0.970621
# ZZZ3         2929.35     -0.2148229  0.147530 -1.456134  0.145356  0.725088

################################################################################

# Cluster 1 v 2,3,4 (all)
gdc.1vAllSamples <- coldata
gdc.1vAllSamples$cluster <- as.character(gdc.1vAllSamples$cluster)  # convert to char class
gdc.1vAllSamples$cluster[gdc.1vAllSamples$cluster != "1"] <- "all"
gdc.1vAllSamples$cluster <- as.factor(gdc.1vAllSamples$cluster)  # convert back to factor
gdc.1vAllSamples <- arrange(gdc.1vAllSamples, row.names(gdc.1vAllSamples))
gdc.1vAllCounts <- gdc.counts[,rownames(gdc.1vAllSamples)]
# 
# # Run DESeq on dataset
gdc.1vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.1vAllCounts, colData = gdc.1vAllSamples, design = ~ cluster)
gdc.1vAllSamples.dds <- DESeq(gdc.1vAllSamples.dds)
# 
# # Save the result, order them
gdc.1vAll.results <- results(gdc.1vAllSamples.dds)
gdc.1vAll.results.ordered <- gdc.1vAll.results[order(gdc.1vAll.results$pvalue),]

# Summary
# > summary(gdc.1vAll.results.ordered)
# 
# out of 36710 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 330, 0.9%
# LFC < 0 (down)     : 71, 0.19%
# outliers [1]       : 0, 0%
# low counts [2]     : 12780, 35%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# > gdc.1vAll.results.ordered
# log2 fold change (MLE): cluster all vs 1 
# Wald test p-value: cluster all vs 1 
# DataFrame with 38481 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   RNU1-21P    53.5070       23.07619  2.208495  10.44883 1.48340e-25 3.54994e-21
# PRDM16-DT   74.7964       -4.75515  0.690843  -6.88311 5.85590e-12 7.00688e-08
# DPEP3     1931.8278        6.45099  1.140610   5.65574 1.55175e-08 1.23783e-04
# CRYM        40.6738       -3.03837  0.542219  -5.60358 2.09974e-08 1.25622e-04
# PRDM16     129.7840       -3.25386  0.591058  -5.50515 3.68860e-08 1.76544e-04
# ...             ...            ...       ...       ...         ...         ...
# ZNF885P           0             NA        NA        NA          NA          NA
# ZNF886P           0             NA        NA        NA          NA          NA
# ZNRF3-AS1         0             NA        NA        NA          NA          NA
# ZNRF3-IT1         0             NA        NA        NA          NA          NA
# ZYXP1             0             NA        NA        NA          NA          NA

# # Number of significantly expressed genes (padj < 0.01)
# sum(gdc.1vAll.results$padj < 0.01, na.rm=TRUE)
# # [1] 79 (0.205%)
# sum(gdc.1vAll.results$padj < 0.05, na.rm=TRUE)
# # [1] 248 (0.644%)


################################################################################

# Cluster 2 v 3

# gdc.2v3.samples <- subset(coldata, (coldata$cluster == 2 | coldata$cluster == 3))
# gdc.2v3.samples <- arrange(gdc.2v3.samples, row.names(gdc.2v3.samples))
# gdc.2v3.counts <- gdc.counts[,rownames(gdc.2v3.samples)]
# gdc.2v3.dds <- DESeqDataSetFromMatrix(countData = gdc.2v3.counts, colData = gdc.2v3.samples, design = ~ cluster)
# gdc.2v3.dds <- DESeq(gdc.2v3.dds)

# > results(gdc.2v3.dds)
# log2 fold change (MLE): cluster 3 vs 2 
# Wald test p-value: cluster 3 vs 2 
# DataFrame with 38481 rows and 6 columns
#             baseMean log2FoldChange     lfcSE       stat    pvalue      padj
#            <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
# A1BG        25.71215      0.1482209  0.234864   0.631092 0.5279805  0.801847
# A1BG-AS1   115.40855     -0.0482492  0.233086  -0.207002 0.8360085  0.942474
# A1CF         1.74475      0.6497550  0.545363   1.191417 0.2334901  0.574452
# A2M      24696.74529     -0.4357007  0.218562  -1.993485 0.0462084  0.264525
# A2M-AS1    194.69864      0.4768372  0.219130   2.176048 0.0295517  0.207991
# ...              ...            ...       ...        ...       ...       ...
# ZYG11B       2347.62     0.19510404  0.104622  1.8648494 0.0622025  0.305705
# ZYX         27660.75    -0.28855620  0.157446 -1.8327257 0.0668434  0.316725
# ZYXP1           0.00             NA        NA         NA        NA        NA
# ZZEF1        2535.15    -0.03834353  0.181833 -0.2108727 0.8329866  0.941636
# ZZZ3         2558.67     0.00757182  0.122007  0.0620606 0.9505146  0.983855

################################################################################

# Cluster 2 v 4

# gdc.2v4.samples <- subset(coldata, (coldata$cluster == 2 | coldata$cluster == 4))
# gdc.2v4.samples <- arrange(gdc.2v4.samples, row.names(gdc.2v4.samples))
# gdc.2v4.counts <- gdc.counts[,rownames(gdc.2v4.samples)]
# gdc.2v4.dds <- DESeqDataSetFromMatrix(countData = gdc.2v4.counts, colData = gdc.2v4.samples, design = ~ cluster)
# gdc.2v4.dds <- DESeq(gdc.2v4.dds)
# 
# > results(gdc.2v4.dds)
# log2 fold change (MLE): cluster 4 vs 2 
# Wald test p-value: cluster 4 vs 2 
# DataFrame with 38481 rows and 6 columns
#             baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#            <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# A1BG        23.99985      -0.127693  0.198568 -0.643072 0.5201776  0.824421
# A1BG-AS1   111.71202      -0.214846  0.217495 -0.987820 0.3232409  0.701106
# A1CF         1.68103       0.505808  0.515553  0.981098 0.3265444  0.703998
# A2M      25062.55733      -0.441956  0.183332 -2.410692 0.0159223  0.208795
# A2M-AS1    185.33904       0.258700  0.191391  1.351680 0.1764778  0.564972
# ...              ...            ...       ...       ...       ...       ...
# ZYG11B       2325.22      0.0829704 0.1113444  0.745169 0.4561696  0.787094
# ZYX         28527.52     -0.2589548 0.1520666 -1.702904 0.0885861  0.436930
# ZYXP1           0.00             NA        NA        NA        NA        NA
# ZZEF1        2744.92      0.1039405 0.1744308  0.595884 0.5512526  0.841453
# ZZZ3         2703.54      0.0777118 0.0984926  0.789012 0.4301049  0.770505

################################################################################

# Cluster 2 v 1,3,4

gdc.2vAllSamples <- coldata
gdc.2vAllSamples$cluster <- as.character(gdc.2vAllSamples$cluster)  # convert to char class
gdc.2vAllSamples$cluster[gdc.2vAllSamples$cluster != "2"] <- "all"
gdc.2vAllSamples$cluster <- as.factor(gdc.2vAllSamples$cluster)  # convert back to factor
gdc.2vAllSamples <- arrange(gdc.2vAllSamples, row.names(gdc.2vAllSamples))
gdc.2vAllCounts <- gdc.counts[,rownames(gdc.2vAllSamples)]
# 
# Run DESeq on dataset
gdc.2vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.2vAllCounts, colData = gdc.2vAllSamples, design = ~ cluster)
gdc.2vAllSamples.dds <- DESeq(gdc.2vAllSamples.dds)
# 
# # Save the result, order them
gdc.2vAll.results <- results(gdc.2vAllSamples.dds)
gdc.2vAll.results.ordered <- gdc.2vAll.results[order(gdc.2vAll.results$pvalue),]

# Summary:

# > summary(gdc.2vAll.results.ordered)
# 
# out of 36710 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 612, 1.7%
# LFC < 0 (down)     : 866, 2.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 10650, 29%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# > gdc.2vAll.results.ordered
# log2 fold change (MLE): cluster all vs 2 
# Wald test p-value: cluster all vs 2 
# DataFrame with 38481 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   SLC13A2    129.8274       4.542635  0.622165   7.30134 2.84924e-13 3.85407e-09
# MYBPHL      81.2521       3.304045  0.452838   7.29631 2.95773e-13 3.85407e-09
# MUC2        13.4277      -3.528642  0.490807  -7.18947 6.50433e-13 5.65031e-09
# PYY        152.9412       3.936062  0.577730   6.81298 9.55945e-12 6.22822e-08
# CD33       347.4338      -0.893277  0.134677  -6.63273 3.29542e-11 1.71764e-07
# ...             ...            ...       ...       ...         ...         ...
# ZNF885P           0             NA        NA        NA          NA          NA
# ZNF886P           0             NA        NA        NA          NA          NA
# ZNRF3-AS1         0             NA        NA        NA          NA          NA
# ZNRF3-IT1         0             NA        NA        NA          NA          NA
# ZYXP1             0             NA        NA        NA          NA          NA

# sum(gdc.2vAll.results$padj < 0.01, na.rm=TRUE)
# # [1] 350 (0.9095%)
# sum(gdc.2vAll.results$padj < 0.05, na.rm=TRUE)
# # [1] 907 (2.36%)

################################################################################

# Cluster 3 v 4

# gdc.3v4.samples <- subset(coldata, (coldata$cluster == 3 | coldata$cluster == 4))
# gdc.3v4.samples <- arrange(gdc.3v4.samples, row.names(gdc.3v4.samples))
# gdc.3v4.counts <- gdc.counts[,rownames(gdc.3v4.samples)]
# gdc.3v4.dds <- DESeqDataSetFromMatrix(countData = gdc.3v4.counts, colData = gdc.3v4.samples, design = ~ cluster)
# gdc.3v4.dds <- DESeq(gdc.3v4.dds)

# > results(gdc.3v4.dds)
# log2 fold change (MLE): cluster 4 vs 3 
# Wald test p-value: cluster 4 vs 3 
# DataFrame with 38481 rows and 6 columns
#             baseMean log2FoldChange     lfcSE       stat    pvalue      padj
#            <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
# A1BG        25.80273    -0.27349406  0.237368 -1.1521925  0.249242  0.632097
# A1BG-AS1   112.35602    -0.16684787  0.252322 -0.6612490  0.508453  0.819292
# A1CF         2.01736    -0.14810478  0.523577 -0.2828709  0.777276  0.935460
# A2M      21981.56198    -0.00554852  0.167417 -0.0331418  0.973561  0.993582
# A2M-AS1    220.88350    -0.21760727  0.233499 -0.9319392  0.351368  0.719037
# ...              ...            ...       ...        ...       ...       ...
# ZYG11B       2533.11     -0.1130465  0.110317  -1.024744  0.305484  0.680411
# ZYX         26453.85      0.0304433  0.153144   0.198789  0.842428  0.955709
# ZYXP1           0.00             NA        NA         NA        NA        NA
# ZZEF1        2777.64      0.1426979  0.146783   0.972167  0.330968  0.702124
# ZZZ3         2772.16      0.0688792  0.108518   0.634729  0.525605  0.828846

################################################################################

# Cluster 3 v 1,2,4

gdc.3vAllSamples <- coldata
gdc.3vAllSamples$cluster <- as.character(gdc.3vAllSamples$cluster)  # convert to char class
gdc.3vAllSamples$cluster[gdc.3vAllSamples$cluster != "3"] <- "all"
gdc.3vAllSamples$cluster <- as.factor(gdc.3vAllSamples$cluster)  # convert back to factor
gdc.3vAllSamples <- arrange(gdc.3vAllSamples, row.names(gdc.3vAllSamples))
gdc.3vAllCounts <- gdc.counts[,rownames(gdc.3vAllSamples)]

# Run DESeq2 on dataset
gdc.3vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.3vAllCounts, colData = gdc.3vAllSamples, design = ~ cluster)
gdc.3vAllSamples.dds <- DESeq(gdc.3vAllSamples.dds)

# Save the result, order them
gdc.3vAll.results <- results(gdc.3vAllSamples.dds)
gdc.3vAll.results.ordered <- gdc.3vAll.results[order(gdc.3vAll.results$pvalue),]

# Summary:
# > summary(gdc.3vAll.results.ordered)
# 
# out of 36709 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1198, 3.3%
# LFC < 0 (down)     : 934, 2.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 9941, 27%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Total number of significant genes

sum(gdc.3vAll.results$padj < 0.01, na.rm=TRUE)
# [1] 721 (1.87% significant expressed genes)
sum(gdc.3vAll.results$padj < 0.05, na.rm=TRUE)
# [1] 1453 (3.77%)

# > gdc.3vAll.results.ordered
# log2 fold change (MLE): cluster all vs 3 
# Wald test p-value: cluster all vs 3 
# DataFrame with 38481 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   PYY         133.4244       -4.51289  0.495235  -9.11262 8.04141e-20 2.15269e-15
# EGF          84.2477       -2.64860  0.353776  -7.48667 7.06424e-14 9.45548e-10
# TAP1      10924.0518       -1.20340  0.172116  -6.99182 2.71351e-12 2.42135e-08
# SEZ6L        38.3774        2.80929  0.412352   6.81286 9.56766e-12 6.40315e-08
# COLEC12    1776.0094        1.44501  0.220191   6.56250 5.29135e-11 2.83299e-07
# ...              ...            ...       ...       ...         ...         ...
# ZNF885P            0             NA        NA        NA          NA          NA
# ZNF886P            0             NA        NA        NA          NA          NA
# ZNRF3-AS1          0             NA        NA        NA          NA          NA
# ZNRF3-IT1          0             NA        NA        NA          NA          NA
# ZYXP1              0             NA        NA        NA          NA          NA

################################################################################

# Cluster 4 v 1,2,3

gdc.4vAllSamples <- coldata
gdc.4vAllSamples$cluster <- as.character(gdc.4vAllSamples$cluster)  # convert to char class
gdc.4vAllSamples$cluster[gdc.4vAllSamples$cluster != "4"] <- "all"
gdc.4vAllSamples$cluster <- as.factor(gdc.4vAllSamples$cluster)  # convert back to factor
gdc.4vAllSamples <- arrange(gdc.4vAllSamples, row.names(gdc.4vAllSamples))
gdc.4vAllCounts <- gdc.counts[,rownames(gdc.4vAllSamples)]

# Run DESeq on dataset
gdc.4vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.4vAllCounts, colData = gdc.4vAllSamples, design = ~ cluster)
gdc.4vAllSamples.dds <- DESeq(gdc.4vAllSamples.dds)

# Save the result, order them
gdc.4vAll.results <- results(gdc.4vAllSamples.dds)
gdc.4vAll.results.ordered <- gdc.4vAll.results[order(gdc.4vAll.results$pvalue),]

# Summary:

# > summary(gdc.4vAll.results.ordered)
# 
# out of 36709 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 336, 0.92%
# LFC < 0 (down)     : 156, 0.42%
# outliers [1]       : 0, 0%
# low counts [2]     : 10651, 29%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Total number of significant genes

sum(gdc.4vAll.results$padj < 0.01, na.rm=TRUE)
# [1] 48 (0.125%)
sum(gdc.4vAll.results$padj < 0.05, na.rm=TRUE)
# [1] 215 (0.559%)

# > results(gdc.4vAllSamples.dds)
# log2 fold change (MLE): cluster all vs 4 
# Wald test p-value: cluster all vs 4 
# DataFrame with 38481 rows and 6 columns
#             baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#            <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# A1BG        25.95116      0.2677714  0.185617  1.442601  0.149133  0.571823
# A1BG-AS1   114.31391      0.2165177  0.191011  1.133532  0.256991  0.690339
# A1CF         1.77696     -0.1612413  0.420149 -0.383772  0.701148  0.926105
# A2M      23368.17184      0.1960864  0.157768  1.242878  0.213913  0.647855
# A2M-AS1    196.11169     -0.0376423  0.183301 -0.205358  0.837292  0.964668
# ...              ...            ...       ...       ...       ...       ...
# ZYG11B       2392.13      0.0194952 0.0877599  0.222142  0.824203  0.961103
# ZYX         27415.20      0.1189274 0.1289442  0.922317  0.356363  0.771042
# ZYXP1           0.00             NA        NA        NA        NA        NA
# ZZEF1        2675.75     -0.1145592 0.1353564 -0.846352  0.397356  0.797625
# ZZZ3         2710.68     -0.0326437 0.0908547 -0.359295  0.719374  0.931215

################################################################################

