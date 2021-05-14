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
gdc.files <- intersect(gdc.files, gdc.valid.files$V1)

# reading in data
data.dir <- '/home/stephen/Documents/classes/bme/237/final_project/data/gdc/individual/'
tmp.files <- lapply(gdc.files$V1, function(x) {
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
#ds <- ds[,colnames(ds) %in% valid.samples]
ds <- ds[, colnames(ds) %in% colnames(exp.gdc)]
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
enst = rownames(exp.kal)
enst.no_version = sapply(strsplit(as.character(enst),"\\."),"[[",1)
g.name <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol'), filters = 'ensembl_transcript_id', values=enst.no_version, mart=ensembl)
head(g.name)

cts <- exp.kal # exp.gdc 
# change g.name getBM depending on dataset
# exp.gdc uses ensembl_gene_id
# exp.kal uses ensembl_transcript_id

# convert to ensembl_transcript_id to hgnc_symbol names
cts <- as.data.frame(cbind(cts, enst.no_version))
cts$enst.no_version <- lvls_revalue(factor(cts$enst.no_version, levels = g.name$ensembl_transcript_id), g.name$hgnc_symbol)

#cts$enst.no_version <- lvls_revalue(factor(cts$enst.no_version, levels = g.name$ensembl_transcript_id), g.name$hgnc_symbol)
# recode into gene symbol
# recode(unlist(lapply(tmp.files, names)), !!!gdc.files.recode)
# ^ something like above line, with gdc.files.recode as g.name as a list
# if works, do same for Kallisto data

cts <- na.omit(cts)
head(cts)

# aggregate duplicate gene and bind to integer count
cts.tmp <- aggregate(apply(cts[-7], 2, as.numeric), cts["enst.no_version"], sum)
cts <- cts.tmp[-1, -1]
cts <- apply(cts, 2, as.integer)
rownames(cts) <- cts.tmp$enst.no_version[-1]
head(cts)


################################################################################

enst = rownames(exp.gdc)
enst.no_version = sapply(strsplit(as.character(enst),"\\."),"[[",1)
g.name <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values=enst.no_version, mart=ensembl)
head(g.name)

cts <- exp.gdc

# convert to ensembl_gene_id to hgnc_symbol names
cts <- as.data.frame(cbind(cts, enst.no_version))
cts$enst.no_version <- lvls_revalue(factor(cts$enst.no_version, levels = g.name$ensembl_gene_id), g.name$hgnc_symbol)
cts <- na.omit(cts)
head(cts)

# aggregate duplicate gene and bind to integer count
cts.tmp <- aggregate(apply(cts[-7], 2, as.numeric), cts["enst.no_version"], sum)
cts <- cts.tmp[-1, -1]
cts <- apply(cts, 2, as.integer)
rownames(cts) <- cts.tmp$enst.no_version[-1]
head(cts)




