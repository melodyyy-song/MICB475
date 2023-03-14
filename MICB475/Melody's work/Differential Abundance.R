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

# Filter dataset by age range
colombia_young <- subset_samples(colombia_RA, age_range=="18_40")
colombia_old <- subset_samples(colombia_RA, age_range=="41_62")

# What ASVs are found in more than 50% of samples in each age range category?
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

####Core Microbiome for young Insuling sensitivity####
# Filter dataset by insulin sensitivity
colombia_young_IR <- subset_samples(colombia_young, insulin_resistance=="yes")
colombia_young_IR_Neg <- subset_samples(colombia_young, insulin_resistance=="no")

# What ASVs are found in more than 50% of samples in each age range category?
young_IR_ASVs <- core_members(colombia_young_IR, detection=0, prevalence = 0.5)
young_IR_Neg_ASVs <- core_members(colombia_young_IR_Neg, detection=0, prevalence = 0.5)

# What are these ASVs?
prune_taxa(young_IR_ASVs,colombia_final) %>%
  tax_table()

tax_table(prune_taxa(young_IR_Neg_ASVs,colombia_final))

prune_taxa(young_IR_Neg_ASVs,colombia_young) %>%
  plot_bar(, fill="Genus")+
  facet_wrap(.~age_range, scales="free")

# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
young_IR_list <- core_members(colombia_young_IR, detection=0.001, prevalence = 0.10)
young_IR_Neg_list <- core_members(colombia_young_IR_Neg, detection=0.001, prevalence = 0.10)

list(Young_IR = young_IR_list, Young_IR_Neg = young_IR_Neg_list)

ggVennDiagram(x=list(Young_IR = young_IR_list, Young_IR_Neg = young_IR_Neg_list)
              , filename = "venndiagram_young_and_old.png", output=TRUE) +
  labs(fill="Count") 

####Core Microbiome for old Insuling sensitivity####
# Filter dataset by insulin sensitivity
colombia_old_IR <- subset_samples(colombia_old, insulin_resistance=="yes")
colombia_old_IR_Neg <- subset_samples(colombia_old, insulin_resistance=="no")

# What ASVs are found in more than 50% of samples in each age range category?
old_IR_ASVs <- core_members(colombia_old_IR, detection=0, prevalence = 0.5)
old_IR_Neg_ASVs <- core_members(colombia_old_IR_Neg, detection=0, prevalence = 0.5)

# What are these ASVs?
prune_taxa(old_IR_ASVs,colombia_final) %>%
  tax_table()

tax_table(prune_taxa(old_IR_Neg_ASVs,colombia_final))

prune_taxa(old_IR_Neg_ASVs,colombia_old) %>%
  plot_bar(, fill="Genus")+
  facet_wrap(.~age_range, scales="free")

# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
old_IR_list <- core_members(colombia_old_IR, detection=0.001, prevalence = 0.10)
old_IR_Neg_list <- core_members(colombia_old_IR_Neg, detection=0.001, prevalence = 0.10)

list(old_IR = old_IR_list, old_IR_Neg = old_IR_Neg_list)

ggVennDiagram(x=list(old_IR = old_IR_list, old_IR_Neg = old_IR_Neg_list)
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

##Need Taxa table to tell what the ASVs are?
#### DESeq by age####
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
  mutate(significant = padj<0.01 & abs(log2FoldChange)>1.5) %>%
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

#### DESeq for young by IR####
colombia_young_deseq <- phyloseq_to_deseq2(colombia_young, ~insulin_sensitivity)
DESEQ_colombia_young <- DESeq(colombia_young_deseq)

## NOTE: If you get a zeros error, then you need to add '1' count to all reads
colombia_young_plus1 <- transform_sample_counts(colombia_final, function(x) x+1)
colombia_young_deseq <- phyloseq_to_deseq2(colombia_young_plus1, ~insulin_sensitivity)
DESEQ_colombia_young <- DESeq(colombia_young_deseq)
res <- results(DESEQ_colombia_young, tidy=TRUE)
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
colombia_young_DESeq <- prune_taxa(sigASVs_vec,colombia_final)
sigASVs <- tax_table(colombia_young_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#### DESeq for old by IR####
colombia_old_deseq <- phyloseq_to_deseq2(colombia_old, ~insulin_sensitivity)
DESEQ_colombia_old <- DESeq(colombia_old_deseq)

## NOTE: If you get a zeros error, then you need to add '1' count to all reads
colombia_old_plus1 <- transform_sample_counts(colombia_final, function(x) x+1)
colombia_old_deseq <- phyloseq_to_deseq2(colombia_old_plus1, ~insulin_sensitivity)
DESEQ_colombia_old <- DESeq(colombia_old_deseq)
res <- results(DESEQ_colombia_old, tidy=TRUE)
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
colombia_old_DESeq <- prune_taxa(sigASVs_vec,colombia_final)
sigASVs <- tax_table(colombia_old_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
