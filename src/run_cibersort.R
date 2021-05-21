home.dir <- '/hb/home/sjhwang/bme237/BME_237_final_project'
setwd(home.dir)
source('bin/cibersort/CIBERSORT.R')

# run Cibersort w/
# (1) Signature file (for example LM22.txt)
# (2) Mixture file (the expression matrix for example ExampleMixtures-GEPs.txt)
# (3) 100 permutations (to get a p-value)
# (4) quantile normalization = FALSE (recommended for RNA-Seq data)
# rest default, which basically means to not run 'absolute' mode

signature_matrix <- 'bin/cibersort/LM22.txt'
mixture_file <- 'data/exp_kal_recoded_aggr.txt'
mixture_file <- 'data/exp_gdc_recoded_aggr.txt'

results <- CIBERSORT(signature_matrix, mixture_file, 1, FALSE)

# write result table to file
write.table(results, file="analysis/cb_gdc_resuts.txt", sep = "\t", quote = F,
            row.names = T, col.names = T)
