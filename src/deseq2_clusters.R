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

gdc.2vAllSamples <- coldata
gdc.2vAllSamples$cluster <- as.character(gdc.2vAllSamples$cluster)  # convert to char class
gdc.2vAllSamples$cluster[gdc.2vAllSamples$cluster != "2"] <- "all"
gdc.2vAllSamples$cluster <- as.factor(gdc.2vAllSamples$cluster)  # convert back to factor
gdc.2vAllSamples <- arrange(gdc.2vAllSamples, row.names(gdc.2vAllSamples))
gdc.2vAllCounts <- gdc.counts[,rownames(gdc.2vAllSamples)]
gdc.2vAllSamples.dds <- DESeqDataSetFromMatrix(countData = gdc.2vAllCounts, colData = gdc.2vAllSamples, design = ~ cluster)
gdc.2vAllSamples.dds <- DESeq(gdc.2vAllSamples.dds)

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

