#changing files into phyloseq object

# load packages

library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

# load data

metafp <- "colombia_metadata.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$`#SampleID`
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
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
# Remove non-bacterial sequences, if any
colombia_young_filt <- subset_taxa(colombia_young,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
colombia_young_filt_nolow <- filter_taxa(colombia_young_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
colombia_final <- prune_samples(sample_sums(colombia_young_filt_nolow)>100, colombia_young_filt_nolow)


# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
rarecurve(t(as.data.frame(otu_table(colombia_final))), cex=0.1)
colombia_rare <- rarefy_even_depth(colombia_final, rngseed = 1, sample.size = 20000)


##### Saving #####
write.table()
save(colombia_final, file="colombia_final.RData")
save(colombia_rare, file="colombia_rare.RData")

# Alpha Diversity

plot_richness(atacama_rare) 

plot_richness(atacama_rare, measures = c("Shannon","Chao1")) 

gg_richness <- plot_richness(atacama_rare, x = "vegetation", measures = c("Shannon","Chao1")) +
  xlab("Vegetation present") +
  geom_boxplot()
gg_richness

ggsave(filename = "Module13/plot_richness.png"
       , gg_richness
       , height=4, width=6)

estimate_richness(atacama_rare)
# ggpubr