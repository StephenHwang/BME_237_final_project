# Differential Gene Analysis across clusters from kal and Kallisto data
library(DESeq2)
library(tidyverse)
require(ggplot2)
require(ggrepel)

home.dir <- '/home/stephen/Documents/classes/bme/237/final_project'
setwd(home.dir)

################################################################################
theme_tufte <- function (base_size = 19, base_family = "Arial", ticks = FALSE) #"serif"
{
  ret <- theme_bw(base_family = base_family, base_size = base_size) +
    theme(legend.background = element_blank(), legend.key = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(),
          strip.background = element_blank(), plot.background = element_blank(),
          axis.line  = element_line(colour = "black", size = rel(1)), #axis.line = element_blank(),
          panel.grid = element_blank(),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
          panel.spacing = unit(1, "lines"),
          plot.margin=unit(c(10,10,5,0),"mm"),
          axis.title = element_text(size = 19, face = "plain", family = "Arial"),
          strip.text.x = element_text(colour = "black", face = "plain", family = "Arial", size = 19))
  ret
}

################################################################################
## Load data

# Load kal clusters
kal.clusters <- readRDS(paste0(home.dir, "/analysis/kal_km_4.rds"))

coldata <- as.data.frame(kal.clusters$cluster, colnames = "Cluster")
colnames(coldata) <- c("cluster")
coldata <- arrange(coldata, row.names(coldata))
coldata$cluster <- factor(coldata$cluster)

# Load kal count
# i-th row and the j-th column of the matrix tells how many reads can be assigned to gene i in sample j.
kal.counts <- readRDS(paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))

# only retain the 105 variables included in the clusters
kal.counts <- kal.counts[,rownames(coldata)]

################################################################################
# Cluster 1 v 2,3,4 (all)
kal.1vAllSamples <- coldata
kal.1vAllSamples$cluster <- as.character(kal.1vAllSamples$cluster)  # convert to char class
kal.1vAllSamples$cluster[kal.1vAllSamples$cluster != "1"] <- "all"
kal.1vAllSamples$cluster <- as.factor(kal.1vAllSamples$cluster)  # convert back to factor
kal.1vAllSamples <- arrange(kal.1vAllSamples, row.names(kal.1vAllSamples))
kal.1vAllCounts <- kal.counts[,rownames(kal.1vAllSamples)]

# # Run DESeq on dataset
kal.1vAllSamples.dds <- DESeqDataSetFromMatrix(countData = kal.1vAllCounts, colData = kal.1vAllSamples, design = ~ cluster)
kal.1vAllSamples.dds <- DESeq(kal.1vAllSamples.dds)

# # Save the result, order them
kal.1vAll.results <- results(kal.1vAllSamples.dds)
kal.1vAll.results.ordered <- kal.1vAll.results[order(kal.1vAll.results$pvalue),]

################################################################################
# Cluster 2 v 1,3,4

kal.2vAllSamples <- coldata
kal.2vAllSamples$cluster <- as.character(kal.2vAllSamples$cluster)  # convert to char class
kal.2vAllSamples$cluster[kal.2vAllSamples$cluster != "2"] <- "all"
kal.2vAllSamples$cluster <- as.factor(kal.2vAllSamples$cluster)  # convert back to factor
kal.2vAllSamples <- arrange(kal.2vAllSamples, row.names(kal.2vAllSamples))
kal.2vAllCounts <- kal.counts[,rownames(kal.2vAllSamples)]
#
# Run DESeq on dataset
kal.2vAllSamples.dds <- DESeqDataSetFromMatrix(countData = kal.2vAllCounts, colData = kal.2vAllSamples, design = ~ cluster)
kal.2vAllSamples.dds <- DESeq(kal.2vAllSamples.dds)
#
# # Save the result, order them
kal.2vAll.results <- results(kal.2vAllSamples.dds)
kal.2vAll.results.ordered <- kal.2vAll.results[order(kal.2vAll.results$pvalue),]

################################################################################
# Cluster 3 v 1,2,4

kal.3vAllSamples <- coldata
kal.3vAllSamples$cluster <- as.character(kal.3vAllSamples$cluster)  # convert to char class
kal.3vAllSamples$cluster[kal.3vAllSamples$cluster != "3"] <- "all"
kal.3vAllSamples$cluster <- as.factor(kal.3vAllSamples$cluster)  # convert back to factor
kal.3vAllSamples <- arrange(kal.3vAllSamples, row.names(kal.3vAllSamples))
kal.3vAllCounts <- kal.counts[,rownames(kal.3vAllSamples)]

# Run DESeq2 on dataset
kal.3vAllSamples.dds <- DESeqDataSetFromMatrix(countData = kal.3vAllCounts, colData = kal.3vAllSamples, design = ~ cluster)
kal.3vAllSamples.dds <- DESeq(kal.3vAllSamples.dds)

# Save the result, order them
kal.3vAll.results <- results(kal.3vAllSamples.dds)
kal.3vAll.results.ordered <- kal.3vAll.results[order(kal.3vAll.results$pvalue),]

################################################################################
# Cluster 4 v 1,2,3

kal.4vAllSamples <- coldata
kal.4vAllSamples$cluster <- as.character(kal.4vAllSamples$cluster)  # convert to char class
kal.4vAllSamples$cluster[kal.4vAllSamples$cluster != "4"] <- "all"
kal.4vAllSamples$cluster <- as.factor(kal.4vAllSamples$cluster)  # convert back to factor
kal.4vAllSamples <- arrange(kal.4vAllSamples, row.names(kal.4vAllSamples))
kal.4vAllCounts <- kal.counts[,rownames(kal.4vAllSamples)]

# Run DESeq on dataset
kal.4vAllSamples.dds <- DESeqDataSetFromMatrix(countData = kal.4vAllCounts, colData = kal.4vAllSamples, design = ~ cluster)
kal.4vAllSamples.dds <- DESeq(kal.4vAllSamples.dds)

# Save the result, order them
kal.4vAll.results <- results(kal.4vAllSamples.dds)
kal.4vAll.results.ordered <- kal.4vAll.results[order(kal.4vAll.results$pvalue),]

################################################################################
# saving DESeq2 results as RDS files
#    baseMean = the average of the normalized counts taken over all samples
#    log2FoldChange = log2 fold change between the groups.
#       E.g. value 2 means that the expression has increased 4-fold
#    lfcSE = standard error of the log2FoldChange estimate
#    stat = Wald statistic
#    pvalue = Wald test p-value
#    padj = Benjamini-Hochberg adjusted p-value

#saveRDS(kal.1vAll.results.ordered, file = paste0(home.dir, '/analysis/kal_DE_results/kal_de_1vAll.rds'))
#saveRDS(kal.2vAll.results.ordered, file = paste0(home.dir, '/analysis/kal_DE_results/kal_de_2vAll.rds'))
#saveRDS(kal.3vAll.results.ordered, file = paste0(home.dir, '/analysis/kal_DE_results/kal_de_3vAll.rds'))
#saveRDS(kal.4vAll.results.ordered, file = paste0(home.dir, '/analysis/kal_DE_results/kal_de_4vAll.rds'))

# load DESeq2 data
kal.1vAll.results <- readRDS(paste0(home.dir, '/analysis/kal_DE_results/kal_de_1vAll.rds'))
kal.2vAll.results <- readRDS(paste0(home.dir, '/analysis/kal_DE_results/kal_de_2vAll.rds'))
kal.3vAll.results <- readRDS(paste0(home.dir, '/analysis/kal_DE_results/kal_de_3vAll.rds'))
kal.4vAll.results <- readRDS(paste0(home.dir, '/analysis/kal_DE_results/kal_de_4vAll.rds'))

kal.1vAll.results.ordered <- kal.1vAll.results[order(kal.1vAll.results$pvalue),]
kal.2vAll.results.ordered <- kal.2vAll.results[order(kal.2vAll.results$pvalue),]
kal.3vAll.results.ordered <- kal.3vAll.results[order(kal.3vAll.results$pvalue),]
kal.4vAll.results.ordered <- kal.4vAll.results[order(kal.4vAll.results$pvalue),]

################################################################################
de.colors <- c("blue", "red", "#9d9d9d")
names(de.colors) <- c("Down", "Up", "None")

## volcano plot
p.val.threshold <- 0.05
reg.threshold <- 0.5

data <- kal.1vAll.results.ordered
#data <- kal.2vAll.results.ordered
#data <- kal.3vAll.results.ordered
#data <- kal.4vAll.results.ordered

data <- data[!is.na(data$padj),]

#data <- cbind(data[,c('log2FoldChange', 'padj')], 'None')
data <- cbind(data[,c('log2FoldChange', 'pvalue')], 'None')


colnames(data) <- c('log2.fc', 'padj', 'diffExp')
data$neg.log10.pval <- -log10(data$padj)

# if pvalue < 0.05 and log2(fc) >0.6 or < -0.6 → up or down regulated
data$diffExp[data$log2.fc > reg.threshold & data$padj <= p.val.threshold] <- "Up"
data$diffExp[data$log2.fc < -reg.threshold & data$padj <= p.val.threshold] <- "Down"

# labels for all DE genes
data$de.label <- NA
data$de.label[data$diffExp != "None"] <- rownames(data)[data$diffExp != "None"]

# distance from origin
data$dist <- sqrt(data$log2.fc^2 +data$neg.log10.pval^2)

# top upregulated genes
upreg.genes <- data[data$diffExp == 'Up',]
upreg.genes <- upreg.genes[order(-upreg.genes$dist),][1:5,]
upreg.genes <- upreg.genes[order(upreg.genes$de.label),] # alphabetical ordering

# top downregulated genes
downreg.genes <- data[data$diffExp == 'Down',]
downreg.genes <- downreg.genes[order(-downreg.genes$dist),][1:5,] # getting top 10
downreg.genes <- downreg.genes[order(downreg.genes$de.label),] # alphabetical ordering

# label top up/down regulated genes
data$label <- NA
data$label[rownames(data) %in% rownames(downreg.genes)] <- downreg.genes$de.label
data$label[rownames(data) %in% rownames(upreg.genes)] <- upreg.genes$de.label

# convert data to data.frame
data <- as.data.frame(data)

## Volcano plot
ggplot(data=data,
       aes(x=log2.fc,
           y=neg.log10.pval,
           col=diffExp,
           label=label)) +
  xlab("log2(fold change)") +
  ylab("-log10(p-value)") +
  theme_tufte() +
  geom_point(alpha=0.7) +
  scale_colour_manual(values = de.colors) +
  geom_text_repel(point.padding = 0.1, box.padding = 0.75, show.legend = F, segment.colour = "black", color='black') +
  geom_vline(xintercept=c(reg.threshold, -reg.threshold), col="#3b3b3b") +
  geom_hline(yintercept=-log10(p.val.threshold), col="#3b3b3b") +
  guides(color = guide_legend(title = "Regulation", position = "right", direction  = "vertical")) +
  labs(title='Kallisto: Cluster 1 vs others') #+
  # xlim(-2.25, 2.25) +
  #ylim(-0.5, 0.75)


# plot using
# width: 800
# height: 675









