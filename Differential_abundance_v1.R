#version 1#
#### Load packages ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggVennDiagram)
library(indicspecies)

####Loading Meta Data with HOMA-IR calculated ####
metafp <- "metadata_edit.txt"
meta <- read_delim(metafp, delim="\t")

otuFP <- "feature-table.txt"
otu <- read.delim(otuFP, skip=1)

taxaFP <- "taxonomy.tsv"
taxa <- read.delim(taxaFP)

phylotreefp <- "tree.nwk"
phylotree <-  read_tree(phylotreefp)
#### Format files to be a phyloseq object ####
### Meta Sample data ###
meta_phylo <- as.data.frame(meta[,-1])
rownames(meta_phylo) <- meta$SampleID
SAMP <- sample_data(meta_phylo)

### OTU table ###
otu_phylo <- otu[,-1]
rownames(otu_phylo) <- otu$X.OTU.ID
OTU <- otu_table(otu_phylo, taxa_are_rows = TRUE)

### Taxonomy ###
taxa_phylo <- taxa %>% 
  separate(Taxon, sep="; ", into=c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  select(-Confidence, -Feature.ID) %>% as.matrix()
rownames(taxa_phylo) <- taxa$Feature.ID
TAXONOMY <- tax_table(taxa_phylo)

### Make phyloseq object ###
colombia_phyObj <- phyloseq(SAMP, OTU, TAXONOMY, phylotree)




######### ANALYZE ##########
# Remove non-bacterial sequences, if any
colombia_filt <- subset_taxa(colombia_phyObj,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
colombia_filt_nolow <- filter_taxa(colombia_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
colombia_final <- prune_samples(sample_sums(colombia_filt_nolow)>100, colombia_filt_nolow)


# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
rarecurve(t(as.data.frame(otu_table(colombia_final))), cex=0.1)
colombia_rare <- rarefy_even_depth(colombia_final, rngseed = 1, sample.size = 500)

# Convert to relative abundance
colombia_RA <- transform_sample_counts(colombia_phyObj, fun=function(x) x/sum(x))

colombia_phyObj <- colombia_RA 
#### Differential Abundance ####

#### Run indicator species analysis by Age_range 18-40 and 41-62 ####
isa_age <- multipatt(t(otu_table(colombia_phyObj)), cluster = sample_data(colombia_phyObj)$age_range)
# Look at results
summary(isa_age)

# Extract taxonomy table
taxtable <- tax_table(colombia_phyObj) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
res_age <- isa_age$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

# View results
View(res_age)

##### ISA insulin resistance ####
# Run indicator species analysis by insulin sensitivity
isa_insulin <- multipatt(t(otu_table(colombia_phyObj)), cluster = sample_data(colombia_phyObj)$insulin_resistance)
# Look at results
summary(isa_insulin)

# Extract taxonomy table
taxtable <- tax_table(colombia_phyObj) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
res_insulin <- isa_insulin$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

# View results
View(res_insulin)

#### DESeqs2 on age-range ####
library(DESeq2)

## NOTE: If you get a zeros error, then you need to add '1' count to all reads
colombia_plus1 <- transform_sample_counts(colombia_phyObj, function(x) x+1)
age_deseq <- phyloseq_to_deseq2(colombia_plus1 , ~age_range)
DESEQ_age <- DESeq(age_deseq)
DeSeqres_age <- results(DESEQ_age, tidy=TRUE)
View(DeSeqres_age)

## Volcano plot: effect size VS significance
ggplot(DeSeqres_age) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
## Make variable to color by whether it is significant + large change
DeSeqres_age %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

# To get table of results
sigASVs_age <- DeSeqres_age %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_age)
# Get only asv names
sigASVs_vec <- sigASVs_age %>%
  pull(ASV)

# Prune phyloseq file
age_DESeq <- prune_taxa(sigASVs_vec,colombia_phyObj)
sigASVs_age <- tax_table(age_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_age) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs_age) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


####DESeq_insulin ####
insulin_deseq <- phyloseq_to_deseq2(colombia_plus1 , ~insulin_resistance)
DESEQ_insulin <- DESeq(insulin_deseq)
DeSeqres_insulin <- results(DESEQ_insulin, tidy=TRUE)
View(DeSeqres_insulin)

## Volcano plot: effect size VS significance
ggplot(DeSeqres_insulin) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
## Make variable to color by whether it is significant + large change
DeSeqres_insulin %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

# To get table of results
sigASVs_insulin <- DeSeqres_insulin %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_age)
# Get only asv names
sigASVs_vec2 <- sigASVs_insulin %>%
  pull(ASV)

# Prune phyloseq file
insulin_DESeq <- prune_taxa(sigASVs_vec2,colombia_phyObj)
sigASVs_insulin <- tax_table(insulin_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_insulin) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs_insulin) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
