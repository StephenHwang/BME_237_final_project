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

gdc.1v2.samples <- subset(coldata, (coldata$cluster == 1 | coldata$cluster == 2))
gdc.1v2.samples <- arrange(gdc.1v2.samples, row.names(gdc.1v2.samples))
gdc.1v2.counts <- gdc.counts[,rownames(gdc.1v2.samples)]
gdc.1v2.dds <- DESeqDataSetFromMatrix(countData = gdc.1v2.counts, colData = gdc.1v2.samples, design = ~ cluster)
gdc.1v2.dds <- DESeq(gdc.1v2.dds)
gdc.1v2.results <- results(gdc.1v2.dds, alpha=0.01)  # cutoff 0.01
gdc.1v2.results.ordered <- gdc.1v2.results[order(gdc.1v2.results$pvalue),]
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
gdc.1v2.results.LFC <- lfcShrink(gdc.1v2.dds, coef="cluster_2_vs_1", type="apeglm")

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
# gdc.1vAllSamples <- coldata
# gdc.1vAllSamples$cluster <- as.character(gdc.1vAllSamples$cluster)  # convert to char class
# gdc.1vAllSamples$cluster[gdc.1vAllSamples$cluster != "1"] <- "all"
# gdc.1vAllSamples$cluster <- as.factor(gdc.1vAllSamples$cluster)  # convert back to factor
# gdc.1vAllSamples <- arrange(gdc.1vAllSamples, row.names(gdc.1vAllSamples))
# gdc.1vAllCounts <- gdc.counts[,rownames(gdc.1vAllSamples)]
# gdc.1vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.1vAllCounts, colData = gdc.1vAllSamples, design = ~ cluster)
# gdc.1vAllSamples.dds <- DESeq(gdc.1vAllSamples.dds)

# > results(gdc.1vAllSamples.dds)
# log2 fold change (MLE): cluster all vs 1 
# Wald test p-value: cluster all vs 1 
# DataFrame with 38481 rows and 6 columns
# baseMean log2FoldChange     lfcSE        stat    pvalue      padj
# <numeric>      <numeric> <numeric>   <numeric> <numeric> <numeric>
#   A1BG        25.95116      -0.508099  0.307542   -1.652129 0.0985083  0.578930
# A1BG-AS1   114.31391      -0.237967  0.322285   -0.738375 0.4602865  0.869163
# A1CF         1.77696       0.225483  0.733876    0.307250 0.7586531        NA
# A2M      23368.17184       0.353718  0.265528    1.332132 0.1828168  0.694574
# A2M-AS1    196.11169       0.298742  0.308219    0.969251 0.3324202  0.810822
# ...              ...            ...       ...         ...       ...       ...
# ZYG11B       2392.13    -0.01580087  0.147894 -0.10683915 0.9149166  0.985357
# ZYX         27415.20     0.00423612  0.217979  0.01943362 0.9844952  0.998515
# ZYXP1           0.00             NA        NA          NA        NA        NA
# ZZEF1        2675.75    -0.00164311  0.228712 -0.00718421 0.9942679  0.999686
# ZZZ3         2710.68    -0.25780171  0.150884 -1.70860903 0.0875234  0.556166

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

# gdc.2vAllSamples <- coldata
# gdc.2vAllSamples$cluster <- as.character(gdc.2vAllSamples$cluster)  # convert to char class
# gdc.2vAllSamples$cluster[gdc.2vAllSamples$cluster != "2"] <- "all"
# gdc.2vAllSamples$cluster <- as.factor(gdc.2vAllSamples$cluster)  # convert back to factor
# gdc.2vAllSamples <- arrange(gdc.2vAllSamples, row.names(gdc.2vAllSamples))
# gdc.2vAllCounts <- gdc.counts[,rownames(gdc.2vAllSamples)]
# gdc.2vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.2vAllCounts, colData = gdc.2vAllSamples, design = ~ cluster)
# gdc.2vAllSamples.dds <- DESeq(gdc.2vAllSamples.dds)

# > results(gdc.2vAllSamples.dds)
# log2 fold change (MLE): cluster all vs 2 
# Wald test p-value: cluster all vs 2 
# DataFrame with 38481 rows and 6 columns
#             baseMean log2FoldChange     lfcSE      stat     pvalue      padj
#            <numeric>      <numeric> <numeric> <numeric>  <numeric> <numeric>
# A1BG        25.95116      0.0759335  0.193751  0.391912 0.69512331 0.8751890
# A1BG-AS1   114.31391     -0.0998295  0.198210 -0.503655 0.61450350 0.8329645
# A1CF         1.77696      0.5350430  0.449159  1.191209 0.23357137 0.5574271
# A2M      23368.17184     -0.4640584  0.157969 -2.937656 0.00330704 0.0735366
# A2M-AS1    196.11169      0.3181019  0.187491  1.696624 0.08976784 0.3696965
# ...              ...            ...       ...       ...        ...       ...
# ZYG11B       2392.13      0.1305730 0.0898558  1.453139  0.1461853  0.456270
# ZYX         27415.20     -0.2636288 0.1313102 -2.007680  0.0446773  0.272483
# ZYXP1           0.00             NA        NA        NA         NA        NA
# ZZEF1        2675.75      0.0365515 0.1403098  0.260506  0.7944734  0.921482
# ZZZ3         2710.68      0.0755771 0.0936856  0.806710  0.4198337  0.712412

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

# gdc.3vAllSamples <- coldata
# gdc.3vAllSamples$cluster <- as.character(gdc.3vAllSamples$cluster)  # convert to char class
# gdc.3vAllSamples$cluster[gdc.3vAllSamples$cluster != "3"] <- "all"
# gdc.3vAllSamples$cluster <- as.factor(gdc.3vAllSamples$cluster)  # convert back to factor
# gdc.3vAllSamples <- arrange(gdc.3vAllSamples, row.names(gdc.3vAllSamples))
# gdc.3vAllCounts <- gdc.counts[,rownames(gdc.3vAllSamples)]
# gdc.3vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.3vAllCounts, colData = gdc.3vAllSamples, design = ~ cluster)
# gdc.3vAllSamples.dds <- DESeq(gdc.3vAllSamples.dds)
# 
# > results(gdc.3vAllSamples.dds)
# log2 fold change (MLE): cluster all vs 3 
# Wald test p-value: cluster all vs 3 
# DataFrame with 38481 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat    pvalue      padj
# <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
#   A1BG        25.95116     -0.1320719  0.194368 -0.679492  0.496826  0.802868
# A1BG-AS1   114.31391     -0.0288409  0.200309 -0.143982  0.885515  0.965983
# A1CF         1.77696     -0.3621670  0.434985 -0.832596  0.405072  0.748467
# A2M      23368.17184      0.1752992  0.164930  1.062872  0.287840  0.657431
# A2M-AS1    196.11169     -0.3551225  0.188281 -1.886133  0.059277  0.318452
# ...              ...            ...       ...       ...       ...       ...
# ZYG11B       2392.13     -0.1426942 0.0905260 -1.576278  0.114962  0.440843
# ZYX         27415.20      0.1526490 0.1343022  1.136608  0.255702  0.625413
# ZYXP1           0.00             NA        NA        NA        NA        NA
# ZZEF1        2675.75      0.0911641 0.1414611  0.644447  0.519286  0.815904
# ZZZ3         2710.68      0.0667084 0.0946691  0.704648  0.481029  0.793515

################################################################################

# Cluster 4 v 1,2,3

# gdc.4vAllSamples <- coldata
# gdc.4vAllSamples$cluster <- as.character(gdc.4vAllSamples$cluster)  # convert to char class
# gdc.4vAllSamples$cluster[gdc.4vAllSamples$cluster != "4"] <- "all"
# gdc.4vAllSamples$cluster <- as.factor(gdc.4vAllSamples$cluster)  # convert back to factor
# gdc.4vAllSamples <- arrange(gdc.4vAllSamples, row.names(gdc.4vAllSamples))
# gdc.4vAllCounts <- gdc.counts[,rownames(gdc.4vAllSamples)]
# gdc.4vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.4vAllCounts, colData = gdc.4vAllSamples, design = ~ cluster)
# gdc.4vAllSamples.dds <- DESeq(gdc.4vAllSamples.dds)
# 
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

