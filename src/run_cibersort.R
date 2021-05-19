
home.dir <- '/home/stephen/Documents/classes/bme/237/final_project'
setwd(home.dir)
source('bin/cibersort/CIBERSORT.R')

# run Cibersort w/ 
# (1) Signature file (for example LM22.txt)
# (2) Mixture file (the expression matrix for example ExampleMixtures-GEPs.txt)
# (3) 100 permutations (to get a p-value)
# (4) quantile normalization = FALSE (recommended for RNA-Seq data)
# rest default, which basically means to not run 'absolute' mode

signature_matrix <- '/home/stephen/Documents/classes/bme/237/final_project/bin/cibersort/LM22.txt' 
mixture_file <- '/home/stephen/Documents/classes/bme/237/final_project/bin/cibersort/LM22.txt' 


results <- CIBERSORT('./CIBERSORT_package/LM22.txt','./CIBERSORT_package/ExampleMixtures-GEPs.txt', 100, FALSE)

# write result table to file
write.table(results,file="./Cibersort_resuts.txt",sep = "\t", quote = F
            ,row.names = T, col.names = T)


