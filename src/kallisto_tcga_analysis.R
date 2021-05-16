require(data.table)
require(dplyr)
require(plyr)

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

library(ggplot2)
library(scales)
#library(ggpubr)

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


# expression boplot
























