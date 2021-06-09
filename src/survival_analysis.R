# Survival Analysis
# Going off of this tutorial: https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html using survminer

# Packages
# install.packages('survival')
# install.packages('survminer')

# Load the library.
library(survival)
library(survminer)
library(ggplot2)

# Load metadata
metadata <- readRDS('meta_data.rds')

# Load the cluster sample names
gdc.cluster.1 <- readRDS('gdc_cluster_1.rds')
gdc.cluster.1.samples <- colnames(gdc.cluster.1)

gdc.cluster.2 <- readRDS('gdc_cluster_2.rds')
gdc.cluster.2.samples <- colnames(gdc.cluster.2)

gdc.cluster.3 <- readRDS('gdc_cluster_3.rds')
gdc.cluster.3.samples <- colnames(gdc.cluster.3)

gdc.cluster.4 <- readRDS('gdc_cluster_4.rds')
gdc.cluster.4.samples <- colnames(gdc.cluster.4)

kal.cluster.1 <- readRDS('kal_cluster_1.rds')
kal.cluster.1.samples <- colnames(kal.cluster.1)

kal.cluster.2 <- readRDS('kal_cluster_2.rds')
kal.cluster.2.samples <- colnames(kal.cluster.2)

kal.cluster.3 <- readRDS('kal_cluster_3.rds')
kal.cluster.3.samples <- colnames(kal.cluster.3)

kal.cluster.4 <- readRDS('kal_cluster_4.rds')
kal.cluster.4.samples <- colnames(kal.cluster.4)

# Add new column for gdc cluster
metadata$STAR_HTSeq_cluster <- vector(mode='integer', length = length(rownames(metadata)))
metadata$kal_cluster <- vector(mode='integer', length = length(rownames(metadata)))

# Set the cluster values
metadata$STAR_HTSeq_cluster[rownames(metadata) %in% gdc.cluster.1.samples] <- 1 
metadata$STAR_HTSeq_cluster[rownames(metadata) %in% gdc.cluster.2.samples] <- 2 
metadata$STAR_HTSeq_cluster[rownames(metadata) %in% gdc.cluster.3.samples] <- 3 
metadata$STAR_HTSeq_cluster[rownames(metadata) %in% gdc.cluster.4.samples] <- 4 

metadata$kal_cluster[rownames(metadata) %in% kal.cluster.1.samples] <- 1 
metadata$kal_cluster[rownames(metadata) %in% kal.cluster.2.samples] <- 2 
metadata$kal_cluster[rownames(metadata) %in% kal.cluster.3.samples] <- 3 
metadata$kal_cluster[rownames(metadata) %in% kal.cluster.4.samples] <- 4 

# Remove rows that do not belong to any cluster
metadata <- metadata[-c(which(metadata$STAR_HTSeq_cluster == 0 & metadata$kal_cluster == 0)),]

#################################
# Using survminer 
#################################

# Comparing survival times between groups
gdc.fit <- survdiff(Surv(metadata$death_days_to, metadata$vital_status == 'Dead') ~ STAR_HTSeq_cluster, data = metadata)
kal.fit <- survdiff(Surv(metadata$death_days_to, metadata$vital_status == 'Dead') ~ kal_cluster, data = metadata)

# Call:
#   survdiff(formula = Surv(metadata$death_days_to, metadata$vital_status == 
#                             "Dead") ~ gdc_cluster, data = metadata)
# 
# n=54, 51 observations deleted due to missingness.
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# gdc_cluster=1  3        3     1.57      1.31     1.402
# gdc_cluster=2 20       20    23.57      0.54     0.992
# gdc_cluster=3 14       14    19.46      1.53     2.576
# gdc_cluster=4 17       17     9.41      6.13     7.878
# 
# Chisq= 10.4  on 3 degrees of freedom, p= 0.02

gdc.fit.all <- survfit(Surv(metadata$death_days_to, metadata$vital_status == 'Dead') ~ STAR_HTSeq_cluster, data = metadata)
kal.fit.all <- survfit(Surv(metadata$death_days_to, metadata$vital_status == 'Dead') ~ kal_cluster, data = metadata)

#fit = gdc.fit.all
fit = kal.fit.all

ggsurvplot(
  fit = fit, 
  data = metadata, 
  xlab = "Days", 
  ylab = "Overall survival probability",
  legend = c(0.8, 0.75))
