require(data.table)
require(dplyr)
require(plyr)
library(ggplot2)
require(ggpubr)

library(tximport)
library(biomaRt)
library(stringr)
library(forcats)


home.dir <- '/home/stephen/Documents/classes/bme/237/final_project'
setwd(home.dir)
# download survival data
# load in data
# filter kallisto for OV data
#    which samples is OV data
# intersect samples
# with data
#    PCA
#    expression boxplot
#

ds <- fread(
  file = './data/kallisto/tcga_Kallisto_est_counts',
  sep = "\t",
  check.names = FALSE,
  data.table = FALSE,
  header = TRUE
)
rownames(ds) <- ds$sample
ds <- ds[-1]
ds_samples <- colnames(ds)

meta.data <- fread(
  file = './data/kallisto/Survival_SupplementalTable_S1_20171025_xena_sp',
  sep = "\t",
  check.names = FALSE,
  data.table = FALSE,
  header = TRUE
)
meta.data <- meta.data[meta.data$`cancer type abbreviation` == 'OV',]
rownames(meta.data) <- meta.data$sample

# sample ids with meta data
ms_samples <- meta.data$sample
sample_ids <- intersect(ms_samples, ds_samples)


################################################################################
# loading GDC data

gdc.valid.files <- fread(
  file = './data/gdc/tcgaFileToSampleID.tsv',
  sep = "\t",
  check.names = FALSE,
  data.table = FALSE,
  header = FALSE
)
gdc.files.recode <- gdc.valid.files$V2
names(gdc.files.recode) <- gdc.valid.files$V1

# all files
gdc.files <- list.files(path='/home/stephen/Documents/classes/bme/237/final_project/data/gdc/individual/')
#gdc.files <- intersect(gdc.files, gdc.valid.files$V1)

# reading in data
data.dir <- '/home/stephen/Documents/classes/bme/237/final_project/data/gdc/individual/'
tmp.files <- lapply(gdc.files, function(x) {
  tmp.file.path <- Sys.glob(file.path(paste0(data.dir, x), "*.htseq.counts.gz"))
  tmp.df <- as.data.frame(read.table(gzfile(tmp.file.path)))
  rownames(tmp.df) <- tmp.df$V1
  tmp.df <- tmp.df[-1]
  colnames(tmp.df) <- c(x)
  return(tmp.df)
  })

# check if all rownames the same
length(unique(unlist(lapply(tmp.files, rownames)))) == nrow(tmp.files[[1]])

d <- unlist(tmp.files)
d <- matrix(d, nrow = nrow(tmp.files[[1]]), byrow = FALSE)
colnames(d) <- recode(unlist(lapply(tmp.files, names)), !!!gdc.files.recode)
rownames(d) <- rownames(tmp.files[[1]])

# filter for intersect
valid.samples <- intersect(sample_ids, colnames(d))  # in both datasets with metadata

d <- d[,colnames(d) %in% valid.samples]
ds <- ds[,colnames(ds) %in% valid.samples]
#ds <- ds[, colnames(ds) %in% colnames(exp.gdc)]
md <- meta.data[meta.data$sample %in% valid.samples, ]

# save as RDS
#saveRDS(d, file = paste0(home.dir, "/data/exp_gdc.rds"))
#saveRDS(ds, file = paste0(home.dir, "/data/exp_kal.rds"))
#saveRDS(md, file = paste0(home.dir, "/data/meta_data.rds"))

# read RDS
exp.gdc <- readRDS(paste0(home.dir, "/data/exp_gdc.rds"))
exp.kal <- readRDS(paste0(home.dir, "/data/exp_kal.rds"))
meta.data <- readRDS(paste0(home.dir, "/data/meta_data.rds"))

# TODO:
# check if datasets are normalized in anyway
#    GDC: not normalized, count data
#    Kal: est_counts (https://xenabrowser.net/datapages/?dataset=tcga_Kallisto_est_counts&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
# convert transcript ID to gene ID

# get hgnc_symbol gene names
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
enst = rownames(exp.gdc)
enst.no_version = sapply(strsplit(as.character(enst),"\\."),"[[",1)
g.name <- as.data.frame(getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values=enst.no_version, mart=ensembl))
head(g.name)
#
# #saveRDS(g.name, file = paste0(home.dir, "/data/kal_recode.rds"))
# #saveRDS(g.name, file = paste0(home.dir, "/data/gdc_recode.rds"))
#
#
# #g.name.list <- g.name$hgnc_symbol
# #names(g.name.list) <- g.name$ensembl_transcript_id
# #g.name.list <- g.name.list[!duplicated(names(g.name.list))]
#
#  # exp.gdc
# # change g.name getBM depending on dataset
# #   exp.gdc uses ensembl_gene_id
# #   exp.kal uses ensembl_transcript_id
#
# # convert to ensembl_transcript_id to hgnc_symbol names
# cts <- exp.kal
# cts <- as.data.frame(cbind(cts, enst.no_version))
# #cts$enst.no_version <- lvls_revalue(factor(cts$enst.no_version, levels = unique(g.name$ensembl_transcript_id)), g.name$hgnc_symbol)
# cts$enst.no_version <- recode(cts$enst.no_version, !!!g.name.list)
#
# cts <- na.omit(cts)
# head(cts)
#
# # aggregate duplicate gene and bind to integer count
# cts.tmp <- aggregate(apply(cts[-7], 2, as.numeric), cts["enst.no_version"], sum)
# cts <- cts.tmp[-1, -1]
# cts <- apply(cts, 2, as.integer)
# rownames(cts) <- cts.tmp$enst.no_version[-1]
# head(cts)

################################################################################
                    ### reload recoded gene name data  ###

# load GDC data
exp.gdc.csv <- fread(
  file = './data/exp_gdc.csv',
  sep = ",",
  check.names = FALSE,
  data.table = FALSE,
  header = TRUE
)
exp.gdc <- aggregate(apply(exp.gdc.csv[-1], 2, as.numeric), exp.gdc.csv["V1"], sum)
rownames(exp.gdc) <- exp.gdc$V1
exp.gdc <- exp.gdc[,-1]
exp.gdc <- exp.gdc[-1,]
#exp.gdc.tmp <- apply(exp.gdc, 2, as.integer)

# load Kallisto data
exp.kal.csv <- fread(
  file = './data/exp_kal.csv',
  sep = ",",
  check.names = FALSE,
  data.table = FALSE,
  header = TRUE
)

exp.kal <- aggregate(apply(exp.kal.csv[-1], 2, as.numeric), exp.kal.csv["V1"], sum)
rownames(exp.kal) <- exp.kal$V1
exp.kal <- exp.kal[,-1]
exp.kal <- exp.kal[-1,]
exp.kal.tmp <- as.data.frame(apply(exp.kal, 2, as.integer))
rownames(exp.kal.tmp) <- rownames(exp.kal)
exp.kal <- exp.kal.tmp

#saveRDS(exp.gdc, file = paste0(home.dir, "/data/exp_gdc_recoded_aggr.rds"))
#saveRDS(exp.kal, file = paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))





################################################################################
# scatter plot

# Read RDS
exp.gdc = readRDS(paste0(home.dir, "/data/exp_gdc_recoded_aggr.rds"))
exp.kal = readRDS(paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))

# plotting theme
theme_tufte <- function(base_size = 19, base_family = "Arial", ticks = FALSE) # "serif"
{
  ret <- theme_bw(base_family = base_family, base_size = base_size) +
    theme(
      legend.background = element_blank(), legend.key = element_blank(),
      panel.background = element_blank(), panel.border = element_blank(),
      strip.background = element_blank(), plot.background = element_blank(),
      axis.line = element_line(colour = "black", size = rel(1)),
      panel.grid = element_blank(),
      axis.text.y = element_text(colour = "black"),
      panel.spacing = unit(1, "lines"),
      plot.margin = unit(c(10, 10, 5, 0), "mm"),
      axis.title = element_text(size = 19, face = "plain", family = "Arial"),
      strip.text.x = element_text(colour = "black", face = "plain", family = "Arial", size = 19),
      plot.caption = element_text(hjust = 1, size = 9) # 0 or 1 for left or right align
    )
  ret
}

all.genes <- c(rownames(exp.gdc), rownames(exp.kal))
all.overlap.genes <- all.genes[duplicated(all.genes)]

exp.gdc.overlap.genes <- exp.gdc[rownames(exp.gdc) %in% all.overlap.genes, ]
exp.kal.overlap.genes <- exp.kal[rownames(exp.kal) %in% all.overlap.genes, ]

# average across rows
exp.gdc.overlap.genes.avg <- rowMeans(exp.gdc.overlap.genes)
exp.gdc.overlap.genes.avg <- rowMeans(exp.kal.overlap.genes)

# bind data
data <- cbind(rownames(exp.kal.overlap.genes), exp.kal.overlap.genes, exp.gdc.overlap.genes.avg)
colnames(data) <- c('gene_id', 'GDC_count', 'KAL_count')

data <- cbind(exp.kal.overlap.genes, exp.gdc.overlap.genes.avg)
colnames(data) <- c('KAL_count', 'GDC_count')
data <- as.data.frame(data)

# plot
ggplot(data, aes(x=log10(GDC_count), y=log10(KAL_count))) +
  theme_tufte() +
  geom_point() +
  labs(
    title = 'GDC_count vs KAL_count gene counts',
    x = 'GDC_count raw counts (log10)',
    y = 'KAL_count raw counts (log10)'
  ) +
  # stat_cor(
  #   method = 'spearman'
  # ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma)

cor.test(
  data$GDC_count,
  data$KAL_count,
  method = 'pearson'
)

ggplot(data, aes(x=GDC_count, y=KAL_count)) +
  theme_tufte() +
  geom_point() +
  labs(
    title = 'GDC_count vs KAL_count gene counts',
    x = 'GDC_count raw counts (log10)',
    y = 'KAL_count raw counts (log10)'
  ) +
  # stat_cor(
  #   method = 'spearman'
  # ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma)




# dataset expression boxplot
exp.kal
exp.gdc

##
data <- data.frame(v1=rnorm(100),v2=rnorm(100),v3=rnorm(100), v4=rnorm(100))
meltData <- melt(data)
ggplot(meltData, aes(factor(variable), value)) +
  geom_boxplot()
##

#
meltData <- melt(exp.gdc[,1:10])
ggplot(meltData, aes(factor(variable), value)) +
  geom_boxplot()
#






p.data <- transpose(data)
colnames(p.data) <- rownames(data) # set sample names
rownames(p.data) <- colnames(data) # set gene names

p.data <- melt(p.data)
p.data <- (cbind(p.data, m.data[rep(seq_len(nrow(m.data)), each = ncol(data)), ]))

m.data <- meta.data
#p.data <- transpose(exp.kal)
#colnames(p.data) <- rownames(exp.kal) # set sample names
#rownames(p.data) <- colnames(exp.kal) # set gene names

p.data <- melt(p.data)
p.data <- (cbind(p.data, m.data[rep(seq_len(ncol(m.data)), each = nrow(p.data)), ]))



ggplot(
  data = data,
  aes_string(
    x = data.x,
    y = data.y
  )
) +
  geom_boxplot(
    outlier.size = 0.5,
    outlier.alpha = 0.8,
    outlier.shape = outlier.shape,
    color = "grey50",
    alpha = 0.6,
    lwd = 0.4
  )







###
exp.gdc = readRDS(paste0(home.dir, "/data/exp_gdc_recoded_aggr.rds"))
exp.kal = readRDS(paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))


write.table(exp.gdc, file="data/exp_txt/exp_gdc_recoded_aggr.txt",sep = "\t", quote = F
            ,row.names = T, col.names = T)

write.table(exp.kal, file="data/exp_txt/exp_kal_recoded_aggr.txt",sep = "\t", quote = F
            ,row.names = T, col.names = T)


################################################################################
###                   Plotting
################################################################################

require(data.table)
require(dplyr)
require(plyr)
library(ggplot2)

home.dir <- '/home/stephen/Documents/classes/bme/237/final_project'
setwd(home.dir)

# plotting theme
theme_tufte <- function(base_size = 19, base_family = "Arial", ticks = FALSE) # "serif"
{
  ret <- theme_bw(base_family = base_family, base_size = base_size) +
    theme(
      legend.background = element_blank(), legend.key = element_blank(),
      panel.background = element_blank(), panel.border = element_blank(),
      strip.background = element_blank(), plot.background = element_blank(),
      axis.line = element_line(colour = "black", size = rel(1)),
      panel.grid = element_blank(),
      axis.text.y = element_text(colour = "black"),
      panel.spacing = unit(1, "lines"),
      plot.margin = unit(c(10, 10, 5, 0), "mm"),
      axis.title = element_text(size = 19, face = "plain", family = "Arial"),
      strip.text.x = element_text(colour = "black", face = "plain", family = "Arial", size = 19),
      plot.caption = element_text(hjust = 1, size = 9) # 0 or 1 for left or right align
    )
  ret
}


# load in data
exp.gdc = readRDS(paste0(home.dir, "/data/exp_gdc_recoded_aggr.rds"))
exp.kal = readRDS(paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))

# sort data
exp.kal <- exp.kal[, order(colnames(exp.kal))]
exp.gdc <- exp.gdc[, order(colnames(exp.gdc))]

# overlap genes
all.genes <- c(rownames(exp.gdc), rownames(exp.kal))
all.overlap.genes <- all.genes[duplicated(all.genes)]

exp.gdc.overlap.genes <- exp.gdc[rownames(exp.gdc) %in% all.overlap.genes, ]
exp.kal.overlap.genes <- exp.kal[rownames(exp.kal) %in% all.overlap.genes, ]


## scatter plot
# average across rows
exp.gdc.avg <- rowMeans(exp.gdc.overlap.genes)
exp.kal.avg <- rowMeans(exp.kal.overlap.genes)

data <- cbind(exp.kal.avg, exp.gdc.avg)
colnames(data) <- c('KAL_count', 'GDC_count')
data <- as.data.frame(data)

# plot
ggplot(data, aes(x=log10(GDC_count), y=log10(KAL_count))) +
  theme_tufte() +
  geom_point() +
  labs(
    title = 'STAR vs Kallisto: raw expression',
    x = 'STAR raw counts (log10)',
    y = 'Kallisto raw counts (log10)'
  ) +
  stat_cor(cor.coef.name = 'rho')

# spearman correlation coefficient of: 0.9207038
cor.test(
  data$GDC_count,
  data$KAL_count,
  # method = 'pearson'
  method = 'spearman'
)


################################################################################
# Dataset boxplot of (237 sig samples)
# data <- exp.gdc.overlap.genes
# data <- exp.kal.overlap.genes

# load in data
exp.gdc = readRDS(paste0(home.dir, "/data/exp_gdc_recoded_aggr.rds"))
exp.kal = readRDS(paste0(home.dir, "/data/exp_kal_recoded_aggr.rds"))
exp.kal <- exp.kal[, order(colnames(exp.kal))] # sort data
exp.gdc <- exp.gdc[, order(colnames(exp.gdc))]
gdc_km_4 = readRDS(paste0(home.dir, "/analysis/gdc_km_4.rds"))
kal_km_4 = readRDS(paste0(home.dir, "/analysis/kal_km_4.rds"))

# filter expression data by final samples
exp.kal <- exp.kal[, colnames(exp.kal) %in% intersect(colnames(exp.kal), names(kal_km_4$cluster))]
exp.kal <- exp.kal[, order(colnames(exp.kal))]
exp.gdc <- exp.gdc[, colnames(exp.gdc) %in% intersect(colnames(exp.gdc), names(gdc_km_4$cluster))]
exp.gdc <- exp.gdc[, order(colnames(exp.kal))]

# load in data
data <- exp.gdc
data <- exp.kal
data <- melt(data)

ggplot(
  data = data,
  aes(factor(variable),
      log10(value+1))) +
  geom_boxplot(
    outlier.size = 0.5,
    outlier.alpha = 0.8,
    outlier.shape = NA,
    color = "grey50",
    alpha = 0.6,
    lwd = 0.4
  ) +
  theme_tufte() +
  labs(
    title = 'Kallisto dataset expression',
    x = 'Dataset',
    y = 'log10 expression'
  ) +
  # scale_y_continuous(limits = quantile(log10(data$value+1), c(0.1, 0.9))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))


