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

# Alpha Diversity Boxplots

# load packages

library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

##### Load RData #####
load("colombia_young_final.RData")
load("colombia_young_rare.RData")
load("colombia_old_final.RData")
load("colombia_old_rare.RData")

# plot young population boxplots
# gg_richness_young <- plot_richness(colombia_young_rare, x = "insulin_resistance", color="insulin_resistance", measures = c("Observed","Shannon","Chao1")) +
#   xlab("Insulin Resistance") +
#   geom_boxplot()
# gg_richness_young

# ggsave(filename = "plot_richness_young.png"
#        , gg_richness_young
#        , height=4, width=6)

gg_richness_young_shannon <- plot_richness(colombia_young_rare, x = "insulin_resistance", color="insulin_resistance", measures = "shannon") +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_young_shannon

ggsave(filename = "plot_richness_young_shannon.png"
       , gg_richness_young_shannon
       , height=4, width=6)


# plot old population boxplots
# gg_richness_old <- plot_richness(colombia_old_rare, x = "insulin_resistance", color="insulin_resistance", measures = c("Observed","Shannon","Chao1")) +
#   xlab("Insulin Resistance") +
#   geom_boxplot()
# gg_richness_old

# ggsave(filename = "plot_richness_old.png"
#        , gg_richness_old
#        , height=4, width=6)

gg_richness_old_shannon <- plot_richness(colombia_old_rare, x = "insulin_resistance", color="insulin_resistance", measures = "Shannon") +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_old_shannon

ggsave(filename = "plot_richness_old_shannon.png"
       , gg_richness_old_shannon
       , height=4, width=6)

# statistical analysis

# t test/wilcoxon test

library(tidyverse)
library(phyloseq)
library(ggpubr)

##### Load RData #####
load("colombia_young_final.RData")
load("colombia_young_rare.RData")
load("colombia_old_final.RData")
load("colombia_old_rare.RData")

# young population t-test (Parametric)

# Need to extract information
alphadiv_young <- estimate_richness(colombia_young_rare)
samp_dat_young <- sample_data(colombia_young_rare)
samp_dat_wdiv_young <- data.frame(samp_dat_young, alphadiv_young)

# t test
t.test(samp_dat_wdiv_young$Shannon ~ samp_dat_wdiv_young$insulin_resistance)

# data is not a normal distribution
allCounts <- as.vector(otu_table(colombia_young_rare))
allCounts <- allCounts[allCounts>0]
hist(allCounts)

# wilcoxon test (non-parametric)
wilcox.test(Shannon ~ insulin_resistance, data=samp_dat_wdiv_young)

# old population statistics
alphadiv_old <- estimate_richness(colombia_old_rare)
samp_dat_old <- sample_data(colombia_old_rare)
samp_dat_wdiv_old <- data.frame(samp_dat_old, alphadiv_old)

# t.test()
t.test(samp_dat_wdiv_old$Shannon ~ samp_dat_wdiv_old$insulin_resistance)

# Wilcoxon
wilcox.test(Shannon ~ insulin_resistance, data=samp_dat_wdiv_old)


# adding p-values to boxplots using wilcoxon test

# young population p-value + boxplots
gg_richness_young <- plot_richness(colombia_young_rare, x = "insulin_resistance", color="insulin_resistance", measures = "Shannon") +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_young

gg_richness_young_pvalue <- gg_richness_young + stat_compare_means(label.y = 5)
gg_richness_young_pvalue

ggsave(filename = "plot_richness_young_pvalue.png"
       , gg_richness_young_pvalue
       , height=6, width=4)

# gg_richness_young_signifvector <- gg_richness_young + stat_compare_means(label = "p.signif", label.x = 1.5) + stat_compare_means(label.y = 5)
# gg_richness_young_signifvector

# ggsave(filename = "plot_richness_young_signifvector.png"
#        , gg_richness_young_signifvector
#        , height=6, width=4)

# old population
gg_richness_old <- plot_richness(colombia_old_rare, x = "insulin_resistance", color="insulin_resistance", measures = "Shannon") +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_old

gg_richness_old_pvalue <- gg_richness_old + stat_compare_means(label.y = 5)
gg_richness_old_pvalue

ggsave(filename = "plot_richness_old_pvalue.png"
       , gg_richness_old_pvalue
       , height=6, width=4)

# gg_richness_old_signifvector <- gg_richness_old + stat_compare_means(label = "p.signif", label.x = 1.5) + stat_compare_means(label.y = 5)
# gg_richness_old_signifvector

# ggsave(filename = "plot_richness_old_signifvector.png"
#        , gg_richness_old_signifvector
#        , height=6, width=4)



