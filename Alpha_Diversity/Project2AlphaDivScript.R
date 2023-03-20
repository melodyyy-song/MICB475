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
rownames(samp_df)<- meta$`SampleID`
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
# Remove non-bacterial sequences
colombia_young_final <- subset_taxa(colombia_young,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
colombia_old_final <- subset_taxa(colombia_young,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
rarecurve(t(as.data.frame(otu_table(colombia_final))), cex=0.1)
colombia_young_rare <- rarefy_even_depth(colombia_young_final, rngseed = 1, sample.size = 20000)
colombia_old_rare <- rarefy_even_depth(colombia_old_final, rngseed = 1, sample.size = 20000)


##### Saving #####
#write.table()
#save(colombia_final, file="colombia_final.RData")
#save(colombia_rare, file="colombia_rare.RData")

# Alpha Diversity young population

gg_richness_young <- plot_richness(colombia_young_rare, x = "insulin_resistance", measures = c("Observed","Shannon","Chao1")) +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_young

ggsave(filename = "Module13/plot_richness_young.png"
       , gg_richness
       , height=4, width=6)

estimate_richness(colombia_young_rare)

gg_richness_old <- plot_richness(colombia_old_rare, x = "insulin_resistance", measures = c("Observed","Shannon","Chao1")) +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_old


# ggpubr

# t test/wilcoxon test

library(tidyverse)
library(phyloseq)

######## Load data ##########
# load("Module13/atacama_rare.RData")

######## Comparison of two means with t-test (Parametric) ##########
# Let's do very simple plot with t-test
plot_richness(colombia_young_rare, x = "insulin_resistance", measures="Shannon")
# Need to extract information
alphadiv <- estimate_richness(colombia_young_rare)
samp_dat <- sample_data(colombia_young_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)
# These are equivalent:
# t.test()
t.test(samp_dat_wdiv$Shannon ~ samp_dat_wdiv$insulin_resistance)


# Note: you can set variances to be equal for a "classic" t-test
# t.test(Shannon ~ vegetation, data=samp_dat_wdiv, var.equal=TRUE)
# If your data is paired, you can set paired=TRUE
# t.test(Shannon ~ vegetation, data=samp_dat_wdiv, paired=TRUE)

#### Microbial count data is generally NON-NORMAL ####
# In fact, it is even more complex because microbial data is usually in RELATIVE ABUNDANCE
allCounts <- as.vector(otu_table(colombia_young_rare))
allCounts <- allCounts[allCounts>0]
hist(allCounts)
hist(log(allCounts))

wilcox.test(Shannon ~ insulin_resistance, data=samp_dat_wdiv)

