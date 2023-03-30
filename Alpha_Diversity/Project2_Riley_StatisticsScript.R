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

gg_richness_young_signifvector <- gg_richness_young + stat_compare_means(label = "p.signif", label.x = 1.5) + stat_compare_means(label.y = 5)
gg_richness_young_signifvector

ggsave(filename = "plot_richness_young_signifvector.png"
       , gg_richness_young_signifvector
       , height=6, width=4)

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

gg_richness_old_signifvector <- gg_richness_old + stat_compare_means(label = "p.signif", label.x = 1.5) + stat_compare_means(label.y = 5)
gg_richness_old_signifvector

ggsave(filename = "plot_richness_old_signifvector.png"
       , gg_richness_old_signifvector
       , height=6, width=4)




