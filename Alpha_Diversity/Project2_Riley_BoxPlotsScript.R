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
gg_richness_young <- plot_richness(colombia_young_rare, x = "insulin_resistance", color="insulin_resistance", measures = c("Observed","Shannon","Chao1")) +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_young

ggsave(filename = "plot_richness_young.png"
       , gg_richness_young
       , height=4, width=6)

gg_richness_young_shannon <- plot_richness(colombia_young_rare, x = "insulin_resistance", color="insulin_resistance", measures = "shannon") +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_young_shannon

ggsave(filename = "plot_richness_young_shannon.png"
       , gg_richness_young_shannon
       , height=4, width=6)


# plot old population boxplots
gg_richness_old <- plot_richness(colombia_old_rare, x = "insulin_resistance", color="insulin_resistance", measures = c("Observed","Shannon","Chao1")) +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_old

ggsave(filename = "plot_richness_old.png"
       , gg_richness_old
       , height=4, width=6)

gg_richness_old_shannon <- plot_richness(colombia_old_rare, x = "insulin_resistance", color="insulin_resistance", measures = "Shannon") +
  xlab("Insulin Resistance") +
  geom_boxplot()
gg_richness_old_shannon

ggsave(filename = "plot_richness_old_shannon.png"
       , gg_richness_old_shannon
       , height=4, width=6)
