# Cancer Phylogenies I
# Exercise 2: Exploring tumours from the Cancer Genome Atlas (TCGA)

# First, install and load the required packages:

#BiocManager::install("TCGAbiolinks")
#BiocManager::install("maftools")

library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(ggplot2)

# load TCGA mutation data for lung cancer:
maf <- GDCquery_Maf("LUAD", pipelines = "mutect2") %>% read.maf
snvs <- data.frame(subsetMaf(maf = maf, tsb =
                               c('TCGA-05-4382-01A-01D-1931-08',
                                 'TCGA-17-Z031-01A-01W-0746-08',
                                 'TCGA-55-8506-01A-11D-2393-08',
                                 'TCGA-55-8513-01A-11D-2393-08'), 
                             mafObj = FALSE))



R.version
updateR