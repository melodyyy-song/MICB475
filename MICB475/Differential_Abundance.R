#Differential Abundance

####Core Microbiome####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
sampdatFP <- "colombia_metadata.txt"
sampdat <- read.delim(sampdatFP)
otuFP <- "Aim 1 qiime2 data/table.qza"
otu <- read.delim(otuFP, skip=1)
taxaFP <- "fecal_taxonomy.txt"
taxa <- read.delim(taxaFP)
