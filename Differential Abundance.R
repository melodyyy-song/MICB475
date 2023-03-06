#Differential Abundance
#For Core Microbiome analysis
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
#For ISA
#library(tidyverse)
#library(phyloseq)
library(indicspecies)
#For DESeq2
#library(tidyverse)
#library(phyloseq)
library(DESeq2)

#### Load data ####
load("colombia_final.RData")
# For this exercise, let's filter out Depth!=2
#atacama_final <- subset_samples(atacama_final, depth==2)

#### "core" microbiome ####

# Convert to relative abundance
colombia_RA <- transform_sample_counts(colombia_final, fun=function(x) x/sum(x))

# Filter dataset by vegetation
colombia_young <- subset_samples(colombia_RA, age_range=="18_40")
colombia_old <- subset_samples(colombia_RA, age_range=="41_62")

# What ASVs are found in more than 50% of samples in each vegetation category?
young_ASVs <- core_members(colombia_young, detection=0, prevalence = 0.5)
old_ASVs <- core_members(colombia_old, detection=0, prevalence = 0.5)

# What are these ASVs?
prune_taxa(young_ASVs,colombia_final) %>%
  tax_table()

tax_table(prune_taxa(old_ASVs,colombia_final))

prune_taxa(old_ASVs,colombia_RA) %>%
  plot_bar(, fill="Genus")+
  facet_wrap(.~age_range, scales="free")

# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
young_list <- core_members(colombia_young, detection=0.001, prevalence = 0.10)
old_list <- core_members(colombia_old, detection=0.001, prevalence = 0.10)

list(Young = young_list, Old = old_list)

ggVennDiagram(x=list(Young = young_list, Old = old_list)
              , filename = "venndiagram_young_and_old.png", output=TRUE) +
  labs(fill="Count") 

#### Indicator Species Analysis ####
# glom to Genus
colombia_genus <- tax_glom(colombia_final, "Genus", NArm = FALSE)
colombia_genus_RA <- transform_sample_counts(colombia_genus, fun=function(x) x/sum(x))
#ISA
isa_colombia <- multipatt(t(otu_table(colombia_genus_RA)), cluster = sample_data(colombia_genus_RA)$age_range)
summary(isa_colombia)
taxtable <- tax_table(colombia_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_colombia$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()

#### DESeq ####
colombia_deseq <- phyloseq_to_deseq2(colombia_final, ~age_range)
DESEQ_colombia <- DESeq(colombia_deseq)

## NOTE: If you get a zeros error, then you need to add '1' count to all reads
colombia_plus1 <- transform_sample_counts(colombia_final, function(x) x+1)
colombia_deseq <- phyloseq_to_deseq2(colombia_plus1, ~age_range)
DESEQ_colombia <- DESeq(colombia_deseq)
res <- results(DESEQ_colombia, tidy=TRUE)
View(res)

# Look at results 

## Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
colombia_DESeq <- prune_taxa(sigASVs_vec,colombia_final)
sigASVs <- tax_table(colombia_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
