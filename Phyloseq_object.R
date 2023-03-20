#Differential Abundance
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ape) 

#### Load data ####
metafp <- "colombia_metadata.txt"
meta <- read.delim(sampdatFP)
otufp <- "18_60-feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)
taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")
phylotreefp <- "rooted-tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$X.SampleID
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
colombia <- phyloseq(OTU, SAMP, TAX, phylotree)

#Rarefy samples to a depth of 1000. Use the rarefied table for all further plots.
#### Filter First ####
colombia_filt <- subset_taxa(colombia,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
colombia_filt_nolow <- filter_taxa(colombia_filt, function(x) sum(x)>5, prune = TRUE)
colombia_filt_nolow_samps <- prune_samples(sample_sums(colombia_filt_nolow)>100, colombia_filt_nolow)
colombia_final <- subset_samples(colombia_filt_nolow_samps, !is.na(ph) )

#### Rarefy samples to depth 1000 ####
atacama_rare <- rarefy_even_depth(atacama_final, rngseed = 1, sample.size = 1000)

save(atacama_final, file="atacama_final.RData")
save(atacama_rare, file="atacama_rare.RData")