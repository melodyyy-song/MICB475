#### Loading Packages ####
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(ggpubr)

#### Loading Data ####
metafp <- "metadata_edit.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU Table ####
# save everything except first column (OTU ID) into a matrix

otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`

# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$`#SampleID`
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

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

#### Create phyloseq object ####
colombia <- phyloseq(OTU, SAMP, TAX, phylotree)

# Remove non-bacterial sequences, if any
colombia_filt <- subset_taxa(colombia,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
rarecurve(t(as.data.frame(otu_table(colombia_filt))), cex=0.1)
colombia_rare <- rarefy_even_depth(colombia_filt, rngseed = 1, sample.size = 20000)

##### Beta Diversity #####
# Beta Diversity: Bray Curtis #
bc_dm <- distance(colombia_rare, method="bray")
pcoa_bc <- ordinate(colombia_rare, method="PCoA", distance=bc_dm)

A <- plot_ordination(colombia_rare, pcoa_bc, color = "age_insulin_resistance") +
  stat_ellipse() + labs(color = "Age Group & Insulin Sensitivity") + 
  scale_color_hue(labels = c("Young & Insulin Resistant", "Young & Insulin Sensitive", "Old & Insulin Resistant", "Old & Insulin Sensitive"))

# Beta Diversity: unweighted unifrac #

unifrac_dm <- distance(colombia_rare, method="unifrac")

pcoa_unifrac <- ordinate(colombia_rare, method="PCoA", distance=unifrac_dm)

B <- plot_ordination(colombia_rare, pcoa_unifrac, color = "age_insulin_resistance") +
  stat_ellipse() + labs(color = "Age Group & Insulin Sensitivity") + 
  scale_color_hue(labels = c("Young & Insulin Resistant", "Young & Insulin Sensitive", "Old & Insulin Resistant", "Old & Insulin Sensitive"))

# Beta Diversity: weighted unifrac #

w_unifrac_dm <- distance(colombia_rare, method="wunifrac")

pcoa_wunifrac <- ordinate(colombia_rare, method="PCoA", distance=w_unifrac_dm)

C <- plot_ordination(colombia_rare, pcoa_wunifrac, color = "age_insulin_resistance") +
  stat_ellipse() + labs(color = "Age Group & Insulin Sensitivity") + 
  scale_color_hue(labels = c("Young & Insulin Resistant", "Young & Insulin Sensitive", "Old & Insulin Resistant", "Old & Insulin Sensitive"))

# Beta Diversity: Jaccard #

jaccard_dm <- distance(colombia_rare, method="jaccard", binary = TRUE)

pcoa_jaccard <- ordinate(colombia_rare, method="PCoA", distance=jaccard_dm)

D <- plot_ordination(colombia_rare, pcoa_jaccard, color = "age_insulin_resistance") +
  stat_ellipse() + labs(color = "Age Group & Insulin Sensitivity") + 
  scale_color_hue(labels = c("Young & Insulin Resistant", "Young & Insulin Sensitive", "Old & Insulin Resistant", "Old & Insulin Sensitive"))

# Arrange and Save Figure #
ggarrange(A, B, C, D,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") %>%
  ggexport(filename = "beta_test.png", width=1440, height=817, res = 100)




##### PERMANOVA #####
samp_dat_wdiv <- data.frame(sample_data(colombia_rare), estimate_richness(colombia_rare))

##bray curtis
#independent and interactive effects
adonis2(bc_dm ~ age_range*insulin_resistance, data=samp_dat_wdiv)
#independent only
adonis2(bc_dm ~ age_range+insulin_resistance, data=samp_dat_wdiv)

## unweighted unifrac
#independent and interactive effects
adonis2(unifrac_dm ~ age_range*insulin_resistance, data=samp_dat_wdiv)
#independent only
adonis2(unifrac_dm ~ age_range+insulin_resistance, data=samp_dat_wdiv)

## weighted unifrac
#independent and interactive effects
adonis2(w_unifrac_dm ~ age_range*insulin_resistance, data=samp_dat_wdiv)
#independent only
adonis2(w_unifrac_dm ~ age_range+insulin_resistance, data=samp_dat_wdiv)

## jaccard
#independent and interactive effects
adonis2(jaccard_dm ~ age_range*insulin_resistance, data=samp_dat_wdiv)
#independent only
adonis2(jaccard_dm ~ age_range+insulin_resistance, data=samp_dat_wdiv)
