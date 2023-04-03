#changing files into phyloseq object

# load packages

library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

# load data

metafp <- "metadata_edit.txt"   #file with calculated_HOMA_IR
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
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
rownames(samp_df)<- meta$`SampleID`
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####

colombia <- phyloseq(OTU, SAMP, TAX, phylotree)

# split sample data into the two age ranges, 18_40 and 41_62

sample_variables(colombia)
unique(SAMP$age_range)

colombia_young <- subset_samples(colombia, age_range == '18_40')
colombia_old <- subset_samples(colombia, age_range == '41_62')

######### ANALYZE ##########
# Remove non-bacterial sequences
colombia_young_final <- subset_taxa(colombia_young,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
colombia_old_final <- subset_taxa(colombia_old,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Rarefy samples
rarecurve(t(as.data.frame(otu_table(colombia_final))), cex=0.1)

#sample size set to 20,000
colombia_young_rare <- rarefy_even_depth(colombia_young_final, rngseed = 1, sample.size = 20000)
colombia_old_rare <- rarefy_even_depth(colombia_old_final, rngseed = 1, sample.size = 20000)


##### Save rarefied files #####
write.table()
save(colombia_young_final, file="colombia_young_final.RData")
save(colombia_young_rare, file="colombia_young_rare.RData")
save(colombia_old_final, file="colombia_old_final.RData")
save(colombia_old_rare, file="colombia_old_rare.RData")
