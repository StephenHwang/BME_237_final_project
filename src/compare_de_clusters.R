# Comparing DE Genes between STAR+HTSeq and Kallisto Clusters

# Load STAR-HTSeq DE clusters
gdc.de.1 <- readRDS('gdc_de_1vAll.rds')
gdc.de.2 <- readRDS('gdc_de_2vAll.rds')
gdc.de.3 <- readRDS('gdc_de_3vAll.rds')
gdc.de.4 <- readRDS('gdc_de_4vAll.rds')

# Load Kallisto DE Clusters
kal.de.1 <- readRDS('kal_de_1vAll.rds')
kal.de.2 <- readRDS('kal_de_2vAll.rds')
kal.de.3 <- readRDS('kal_de_3vAll.rds')
kal.de.4 <- readRDS('kal_de_4vAll.rds')

# Remove genes whose pvalue > 0.05
gdc.de.1 <- gdc.de.1[which(gdc.de.1$pvalue < 0.05),]
gdc.de.2 <- gdc.de.2[which(gdc.de.2$pvalue < 0.05),]
gdc.de.3 <- gdc.de.3[which(gdc.de.3$pvalue < 0.05),]
gdc.de.4 <- gdc.de.4[which(gdc.de.4$pvalue < 0.05),]

kal.de.1 <- kal.de.1[which(kal.de.1$pvalue < 0.05),]
kal.de.2 <- kal.de.2[which(kal.de.2$pvalue < 0.05),]
kal.de.3 <- kal.de.3[which(kal.de.3$pvalue < 0.05),]
kal.de.4 <- kal.de.4[which(kal.de.4$pvalue < 0.05),]

# Load results and order by log2FoldChange
gdc.de.1.ordered <- gdc.de.1[order(gdc.de.1$log2FoldChange),]
gdc.de.2.ordered <- gdc.de.2[order(gdc.de.2$log2FoldChange),]
gdc.de.3.ordered <- gdc.de.3[order(gdc.de.3$log2FoldChange),]
gdc.de.4.ordered <- gdc.de.4[order(gdc.de.4$log2FoldChange),]

kal.de.1.ordered <- kal.de.1[order(kal.de.1$log2FoldChange),]
kal.de.2.ordered <- kal.de.2[order(kal.de.2$log2FoldChange),]
kal.de.3.ordered <- kal.de.3[order(kal.de.3$log2FoldChange),]
kal.de.4.ordered <- kal.de.4[order(kal.de.4$log2FoldChange),]

# Get top 100 genes
gdc.de.1.top50 <- row.names(gdc.de.1.ordered)[1:100]
gdc.de.2.top50 <- row.names(gdc.de.2.ordered)[1:100]
gdc.de.3.top50 <- row.names(gdc.de.3.ordered)[1:100]
gdc.de.4.top50 <- row.names(gdc.de.4.ordered)[1:100]

kal.de.1.top50 <- row.names(kal.de.1.ordered)[1:100]
kal.de.2.top50 <- row.names(kal.de.2.ordered)[1:100]
kal.de.3.top50 <- row.names(kal.de.3.ordered)[1:100]
kal.de.4.top50 <- row.names(kal.de.4.ordered)[1:100]


# Compare Genes between STAR-HTSeq and Kallisto Clusters

# Compare genes between star-htseq cluster 1, kallisto cluster 1 
# 
gdc.1.vs.kal.1 <- intersect(gdc.de.1.top50, kal.de.1.top50)
# Little to no similarity, maybe similar cell fractions between NK cells (which 
# make up a very small percentage of either cluster)
# [1] "SNORA19" "IGHD" 

# Compare genes between star-htseq cluster 1, kallisto cluster 2 
# both enriched for cd4 t cells, would expect similarities in dge
gdc.1.vs.kal.2 <- intersect(gdc.de.1.top50, kal.de.2.top50)
# character(0)

# Compare genes between star-htseq cluster 1, kallisto cluster 3
# fairly similar in CD4 T cells cell fractions
gdc.1.vs.kal.3 <- intersect(gdc.de.1.top50, kal.de.3.top50)
# character(0)

# Compare genes between star-htseq cluster 1, kallisto cluster 4 
gdc.1.vs.kal.4 <- intersect(gdc.de.1.top50, kal.de.4.top50)
# [1] "AJAP1"

# Compare genes between star-htseq cluster 2, kallisto cluster 1 
gdc.2.vs.kal.1 <- intersect(gdc.de.2.top50, kal.de.1.top50)
# [1] "SLC7A10"  "IGLV1-36"

# Compare genes between star-htseq cluster 2, kallisto cluster 2 
gdc.2.vs.kal.2 <- intersect(gdc.de.2.top50, kal.de.2.top50)
# character(0)

# Compare genes between star-htseq cluster 2, kallisto cluster 3 
gdc.2.vs.kal.3 <- intersect(gdc.de.2.top50, kal.de.3.top50)
# [1] "VSTM2B"

# Compare genes between star-htseq cluster 2, kallisto cluster 4 
# both enriched for M2 Macrophages, some enrichment for CD4 T cells
gdc.2.vs.kal.4 <- intersect(gdc.de.2.top50, kal.de.4.top50)
# [1] "DNAI2"   "ANKRD66" "BCAR4" 

# Compare genes between star-htseq cluster 3, kallisto cluster 1 
gdc.3.vs.kal.1 <- intersect(gdc.de.3.top50, kal.de.1.top50)
# [1] "SSX3"  "KCNU1"

# Compare genes between star-htseq cluster 3, kallisto cluster 2 
# similar cell fractions of M2 macrophages
gdc.3.vs.kal.2 <- intersect(gdc.de.3.top50, kal.de.2.top50)
# [1] "LINC02403" "PIP"       "SLC13A2"   "AQP4"      "XIRP2"     "KCNU1"

# Compare genes between star-htseq cluster 3, kallisto cluster 3 
gdc.3.vs.kal.3 <- intersect(gdc.de.3.top50, kal.de.3.top50)
# character(0)

# Compare genes between star-htseq cluster 3, kallisto cluster 4 
gdc.3.vs.kal.4 <- intersect(gdc.de.3.top50, kal.de.4.top50)
# character(0)

# Compare genes between star-htseq cluster 4, kallisto cluster 1 
gdc.4.vs.kal.1 <- intersect(gdc.de.4.top50, kal.de.1.top50)
# character(0)

# Compare genes between star-htseq cluster 4, kallisto cluster 2 
gdc.4.vs.kal.2 <- intersect(gdc.de.4.top50, kal.de.2.top50)
# [1] "C5orf38"

# Compare genes between star-htseq cluster 4, kallisto cluster 3 
gdc.4.vs.kal.3 <- intersect(gdc.de.4.top50, kal.de.3.top50)
# character(0)

# Compare genes between star-htseq cluster 4, kallisto cluster 4 
gdc.4.vs.kal.4 <- intersect(gdc.de.4.top50, kal.de.4.top50)
# [1] "CELA3A" "SCRT2"