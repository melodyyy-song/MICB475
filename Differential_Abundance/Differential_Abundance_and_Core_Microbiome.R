#### Load packages ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(DESeq2)
library(ggpubr)
library(microbiome)
library(ggVennDiagram)

####Loading Meta Data with HOMA-IR calculated ####
metafp <- "metadata_edit.txt"
meta <- read_delim(metafp, delim="\t")

otuFP <- "feature-table.txt"
otu <- read.delim(otuFP, skip=1)

taxaFP <- "taxonomy.tsv"
taxa <- read.delim(taxaFP)

phylotreefp <- "rooted-tree.nwk"
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

#### Combine Genus ####
(genus_colombia_final <- tax_glom(colombia_final, taxrank="Genus") )
ntaxa(colombia_final); ntaxa(genus_colombia_final)

#### Filter dataset by age range ####
genus_colombia_young <- subset_samples(genus_colombia_final, age_range=="18_40")
genus_colombia_old <- subset_samples(genus_colombia_final, age_range=="41_62")

#### DESeqs2 Analyses ####

#18-40 for insulin
colombia_young_plus1 <- transform_sample_counts(genus_colombia_young, function(x) x+1)
colombia_young_deseq <- phyloseq_to_deseq2(colombia_young_plus1, ~insulin_resistance)
DESEQ_colombia_young <- DESeq(colombia_young_deseq)
DeSeqres_young <- results(DESEQ_colombia_young, tidy=TRUE)
View(DeSeqres_young)

#41-62 for insulin
colombia_old_plus1 <- transform_sample_counts(genus_colombia_old, function(x) x+1)
colombia_old_deseq <- phyloseq_to_deseq2(colombia_old_plus1, ~insulin_resistance)
DESEQ_colombia_old <- DESeq(colombia_old_deseq)
DeSeqres_old <- results(DESEQ_colombia_old, tidy=TRUE)
View(DeSeqres_old)

####Volcano plots####
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
#18-40 for insulin
sigASVs_young <- DeSeqres_young %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1) %>%
  dplyr::rename(ASV=row)
View(sigASVs_young)
sigASVs_young_vec <- sigASVs_young %>%
  pull(ASV)
colombia_young_DESeq <- prune_taxa(sigASVs_young_vec,colombia_final)
sigASVs_young <- tax_table(colombia_young_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_young) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(highlight_young = ifelse(abs(log2FoldChange)>1.5, "YT", "YF"))

#41-60 for insulin
sigASVs_old <- DeSeqres_old %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1 ) %>%
  dplyr::rename(ASV=row)
View(sigASVs_old)
sigASVs_old_vec <- sigASVs_old %>%
  pull(ASV)
colombia_old_DESeq <- prune_taxa(sigASVs_old_vec,colombia_final)
sigASVs_old <- tax_table(colombia_old_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_old) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(highlight_old = ifelse(abs(log2FoldChange)>1.5, "OT", "OF"))

####Making the figures have same genuses####
#Find shared Genus names between yung and old age groups
sigASVs_combined <- merge(sigASVs_young, sigASVs_old, by = "Genus")
view(sigASVs_combined) #Find the 2 shared Genus names
#Get dataframe of sigASVs_young without shared Genus names
to_combine_to_old <- sigASVs_young %>%
  filter((!Genus %in% c("g__Butyrivibrio", "f__Erysipelotrichaceae + g__uncl")))
#Dataframe listing all ASVs and genuses shown with sigASVs from both groups
all_combined <- rbind(sigASVs_old[,-15], to_combine_to_old[,-15])

#Updated list of sigASVs_young with all ASVs 
sigASVs_young2 <- DeSeqres_young %>% 
  dplyr::rename(ASV=row) %>%
  filter(ASV %in% all_combined$ASV)
sigASVs_young_vec2 <- sigASVs_young2 %>%
  pull(ASV) 
colombia_young_DESeq2 <- prune_taxa(sigASVs_young_vec2,colombia_final)
sigASVs_young2 <- tax_table(colombia_young_DESeq2) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_young2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(highlight_young = ifelse(abs(log2FoldChange)>1.5, "YT", "YF"))

#List of sigASVs_old with all ASVs
sigASVs_old2 <- DeSeqres_old %>% 
  dplyr::rename(ASV=row) %>%
  filter(ASV %in% all_combined$ASV)
sigASVs_old_vec2 <- sigASVs_old2 %>%
  pull(ASV) 
colombia_old_DESeq2 <- prune_taxa(sigASVs_old_vec2,colombia_final)
sigASVs_old2 <- tax_table(colombia_old_DESeq2) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_old2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(highlight_old = ifelse(abs(log2FoldChange)>1.5, "OT", "OF"))

#### Plotting log2FoldChange ####
#18-40 for insulin
young_bar_plot <- ggplot(sigASVs_young2) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill=highlight_young), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  scale_fill_manual(values = c("YT" = "lightgreen","YF" = "gray"),
                    name = "Significance (18-40 years)",
                    labels = c("no", "yes")) +
  theme(axis.text.x = element_text(hjust=0.5, vjust=0.5)) +
  ylab("log2FoldChange 
       (insulin-resistant/insulin-sensitive)") +
  coord_flip()
young_bar_plot

#41-62 for insulin
old_bar_plot <- ggplot(sigASVs_old2) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill = highlight_old), stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  scale_fill_manual(values = c("OT" = "lightblue","OF" = "gray"), 
                    name = "Significance (41-62 years)",
                    labels = c("no", "yes")) +
  ylab("log2FoldChange 
       (insulin-resistant/insulin-sensitive)") +
  theme(axis.text.x = element_text(hjust=0.5, vjust=0.5) ) +
  coord_flip()
old_bar_plot

####Relative Abundance Boxplots####
list_genus_colombia_old <- genus_colombia_old %>%
  #ps_mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  microbiome::transform("compositional") %>% 
  psmelt()
list_genus_colombia_young <- genus_colombia_young %>%
  microbiome::transform("compositional") %>% 
  psmelt()
#list of significant genuses (to match barplot)
ASV_of_interest_old <- sigASVs_old2$ASV
ASV_of_interest_young <- sigASVs_young2$ASV
#Boxplots
old_box_plot <- list_genus_colombia_old %>% filter(OTU %in% ASV_of_interest_old) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(Abundance = Abundance + min(.$Abundance[.$Abundance>0])/2) %>%
  ggplot(aes(Genus,Abundance,fill=insulin_resistance)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_log10() +
  scale_fill_manual(values = c("yes" = "cornflowerblue","no" = "tomato"),
                     name="Insulin Resistance (41-62 years)")+
  ylab("
       Abundance")
old_box_plot

young_box_plot <- list_genus_colombia_young %>% filter(OTU %in% ASV_of_interest_young) %>%
  mutate(Genus = ifelse(Genus =='g__uncultured', paste(Family,'g__uncl',sep=' + '),Genus)) %>%
  mutate(Abundance = Abundance + min(.$Abundance[.$Abundance>0])/2) %>%
  ggplot(aes(Genus,Abundance,fill=insulin_resistance)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_log10()+
  scale_fill_manual(values = c("yes" = "cornflowerblue","no" = "tomato"),
                    name="Insulin Resistance (18-40 years)")+
  ylab("
       Abundance")
young_box_plot

####Combined Figure####
#Version Final
ggarrange(young_bar_plot + rremove("y.ticks"),
          young_box_plot + rremove("y.text") + rremove("ylab") +rremove("y.ticks"),
          old_bar_plot, #+ rremove("y.text") + rremove("ylab") +rremove("y.ticks"),
          old_box_plot + rremove("y.text") + rremove("ylab") +rremove("y.ticks"),
          labels = "AUTO",
          ncol = 2,
          nrow = 2,
          widths = c(1.5,1),
          legend = "top")
#Use min pixels width:875, and height:1000 to make sure are labels are shown

####VENN DIAGRAM####
colombia_RA <- transform_sample_counts(colombia_final, fun=function(x) x/sum(x))
# Filter dataset by age range
colombia_young_RA <- subset_samples(colombia_RA, age_range=="18_40")
colombia_old_RA <- subset_samples(colombia_RA, age_range=="41_62")

####Core Microbiome for young Insulin sensitivity####
colombia_young_IR <- subset_samples(colombia_young_RA, insulin_resistance=="yes")
colombia_young_IR_Neg <- subset_samples(colombia_young_RA, insulin_resistance=="no")

young_IR_list <- core_members(colombia_young_IR, detection=0.001, prevalence = 0.3)
young_IR_Neg_list <- core_members(colombia_young_IR_Neg, detection=0.001, prevalence = 0.3)
list(Young_IR = young_IR_list, Young_IR_Neg = young_IR_Neg_list)

young_venn <- ggVennDiagram(x=list("Insulin Sensitive" = young_IR_Neg_list, "Insulin Resistant" = young_IR_list)
                            , filename = "venndiagram_young_and_old.png", output=TRUE) +
  labs(fill="Count") +
  scale_fill_distiller(palette = "Blues")+
  ggtitle("Young Age Group (18-40)")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

####Core Microbiome for old Insulin sensitivity####
colombia_old_IR <- subset_samples(colombia_old_RA, insulin_resistance=="yes")
colombia_old_IR_Neg <- subset_samples(colombia_old_RA, insulin_resistance=="no")

old_IR_list <- core_members(colombia_old_IR, detection=0.001, prevalence = 0.3)
old_IR_Neg_list <- core_members(colombia_old_IR_Neg, detection=0.001, prevalence = 0.3)

list(old_IR = old_IR_list, old_IR_Neg = old_IR_Neg_list)

old_venn <- ggVennDiagram(x=list("Insulin Sensitive" = old_IR_Neg_list, "Insulin Resistant" = old_IR_list)
                          , filename = "venndiagram_young_and_old.png", output=TRUE) +
  labs(fill="Count")+
  scale_fill_distiller(palette = "Blues")+
  ggtitle("Old Age Group (41-62)")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#Combined figure
ggarrange(young_venn,
          old_venn,
          common.legend = TRUE,
          legend = "right",
          labels = "AUTO")
