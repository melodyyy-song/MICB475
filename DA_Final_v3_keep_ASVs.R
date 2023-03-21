#### Load packages ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(DESeq2)

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
rownames(meta_phylo) <- meta$`SampleID`
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

#### Filter samples ####
colombia_filt <- subset_taxa(colombia_phyObj,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
colombia_filt_nolow <- filter_taxa(colombia_filt, function(x) sum(x)>5, prune = TRUE)
colombia_filt_nolow_samps <- prune_samples(sample_sums(colombia_filt_nolow)>100, colombia_filt_nolow)
colombia_final <- subset_samples(colombia_filt_nolow_samps, !is.na(insulin_resistance))

##### Filter dataset by age range ####
colombia_young <- subset_samples(colombia_final, age_range=="18_40")
colombia_old <- subset_samples(colombia_final, age_range=="41_62")

#### DESeqs2 Analyses ####
#for Age differences
colombia_plus1 <- transform_sample_counts(colombia_final, function(x) x+1)
age_deseq <- phyloseq_to_deseq2(colombia_plus1 , ~age_range)
DESEQ_age <- DESeq(age_deseq, fitType = "local")
DeSeqres_age <- results(DESEQ_age, tidy=TRUE)
View(DeSeqres_age)

#All ages for Insulin
insulin_deseq <- phyloseq_to_deseq2(colombia_plus1 , ~insulin_resistance)
DESEQ_insulin <- DESeq(insulin_deseq)
DeSeqres_insulin <- results(DESEQ_insulin, tidy=TRUE)
View(DeSeqres_insulin)

#18-40 for insulin
colombia_young_plus1 <- transform_sample_counts(colombia_young, function(x) x+1)
colombia_young_deseq <- phyloseq_to_deseq2(colombia_young_plus1, ~insulin_resistance)
DESEQ_colombia_young <- DESeq(colombia_young_deseq)
DeSeqres_young <- results(DESEQ_colombia_young, tidy=TRUE)
View(DeSeqres_young)

#41-62 for insulin
colombia_old_plus1 <- transform_sample_counts(colombia_old, function(x) x+1)
colombia_old_deseq <- phyloseq_to_deseq2(colombia_old_plus1, ~insulin_resistance)
DESEQ_colombia_old <- DESeq(colombia_old_deseq)
DeSeqres_old <- results(DESEQ_colombia_old, tidy=TRUE)
View(DeSeqres_old)

####Volcano plots####
#For age
ggplot(DeSeqres_age) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
DeSeqres_age %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#For Insulin
ggplot(DeSeqres_insulin) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
DeSeqres_insulin %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#18-40 for Insulin
ggplot(DeSeqres_young) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
DeSeqres_young %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#41-62 for insulin
ggplot(DeSeqres_old) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
DeSeqres_old %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#### Table of results of significant ASVs####
#for age
sigASVs_age <- DeSeqres_age %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_age)
#For ASV names
sigASVs_age_vec <- sigASVs_age %>%
  pull(ASV)
#Prune phyloseq file
age_DESeq <- prune_taxa(sigASVs_age_vec,colombia_final)
sigASVs_age_final <- tax_table(age_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_age) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(!is.na(Genus))
  

#For insulin
sigASVs_insulin <- DeSeqres_insulin %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_insulin)
sigASVs_insulin_vec <- sigASVs_insulin %>%
  pull(ASV)
insulin_DESeq <- prune_taxa(sigASVs_insulin_vec,colombia_final)
sigASVs_insulin_final <- tax_table(insulin_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_insulin) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))%>%
  filter(!is.na(Genus))

#18-40 for insulin
sigASVs_young <- DeSeqres_young %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_young)
sigASVs_young_vec <- sigASVs_young %>%
  pull(ASV)
colombia_young_DESeq <- prune_taxa(sigASVs_young_vec,colombia_final)
sigASVs_young_final <- tax_table(colombia_young_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_young) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))%>%
  filter(!is.na(Genus))

#41-60 for insulin
sigASVs_old <- DeSeqres_old %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_old)
sigASVs_old_vec <- sigASVs_old %>%
  pull(ASV)
colombia_old_DESeq <- prune_taxa(sigASVs_old_vec,colombia_final)
sigASVs_old_final <- tax_table(colombia_old_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_old) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))%>%
  filter(!is.na(Genus))

#### log2FoldChange ####
#For age
ggplot(sigASVs_age_final) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  xlab("Genus of significant ASVs") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#For insulin
ggplot(sigASVs_insulin_final) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  xlab("Genus of significant ASVs") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#18-40 for insulin
ggplot(sigASVs_young_final) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  xlab("Genus of significant ASVs") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#41-62 for insulin
ggplot(sigASVs_old_final) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  xlab("Genus of significant ASVs") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#### Combined Bar graph ####
ggplot() +
  geom_bar(data = sigASVs_old_final, aes(x=Genus, y=log2FoldChange), stat="identity", fill="lightblue")+
  geom_errorbar(data = sigASVs_old_final, aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  geom_bar(data = sigASVs_young_final, aes(x=Genus, y=log2FoldChange), stat="identity", fill="lightgreen")+
  geom_errorbar(data = sigASVs_young_final, aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Log2 Fold Change") + 
  xlab("Genus of significant ASVs") +
  coord_flip()

####Counts of ASVs####
#For age
sigASVs_age_for_table <- tax_table(age_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_age) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(log2FoldChange = ifelse(log2FoldChange <0, paste ('-'),paste ('+'))) %>%
  filter(!is.na(Genus))
sigASVs_age_table <- as.data.frame(table(sigASVs_age_for_table[c("Genus", "log2FoldChange")]))
sigASVs_age_table <- sigASVs_age_table[sigASVs_age_table$Freq %in% 1:10, ] %>%
  arrange(Genus)
view(sigASVs_age_table)

#For insulin
sigASVs_insulin_for_table <- tax_table(insulin_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_insulin) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(log2FoldChange = ifelse(log2FoldChange <0, paste ('-'),paste ('+'))) %>%
  filter(!is.na(Genus))
sigASVs_insulin_table <- as.data.frame(table(sigASVs_insulin_for_table[c("Genus", "log2FoldChange")]))
sigASVs_insulin_table <- sigASVs_insulin_table[sigASVs_insulin_table$Freq %in% 1:10, ] %>%
  arrange(Genus)
view(sigASVs_insulin_table)
#18-40 for insulin
sigASVs_young_for_table <- tax_table(colombia_young_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_young) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(log2FoldChange = ifelse(log2FoldChange <0, paste ('-'),paste ('+'))) %>%
  filter(!is.na(Genus))
sigASVs_young_table <- as.data.frame(table(sigASVs_young_for_table[c("Genus", "log2FoldChange")]))
sigASVs_young_table <- sigASVs_young_table[sigASVs_young_table$Freq %in% 1:10, ] %>%
  arrange(Genus)
view(sigASVs_young_table)

#41-62 for insulin
sigASVs_old_for_table <- tax_table(colombia_old_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_old) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(log2FoldChange = ifelse(log2FoldChange <0, paste ('-'),paste ('+'))) %>%
  filter(!is.na(Genus))
sigASVs_old_table <- as.data.frame(table(sigASVs_old_for_table[c("Genus", "log2FoldChange")]))
sigASVs_old_table <- sigASVs_old_table[sigASVs_old_table$Freq %in% 1:10, ] %>%
  arrange(Genus)
view(sigASVs_old_table)
