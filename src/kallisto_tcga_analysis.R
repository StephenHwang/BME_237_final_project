require(data.table)
require(dplyr)
require(plyr)

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
ds <- ds[,colnames(ds) %in% valid.samples]
md <- meta.data[meta.data$sample %in% valid.samples, ]

# save as RDS
saveRDS(d, file = paste0(home.dir, "/data/exp_gdc.rds"))
saveRDS(ds, file = paste0(home.dir, "/data/exp_kal.rds"))
saveRDS(md, file = paste0(home.dir, "/data/meta_data.rds"))

# read RDS
exp.gdc <- readRDS(paste0(home.dir, "/data/exp_gdc.rds"))
exp.kal <- readRDS(paste0(home.dir, "/data/exp_kal.rds"))
meta.data <- readRDS(paste0(home.dir, "/data/meta_data.rds"))





