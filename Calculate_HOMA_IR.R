#### Loading Packages ####
library(tidyverse)

#### Loading Data ####
metafp <- "colombia_metadata.txt"
meta <- read_delim(metafp, delim="\t")

# Calculate HOMA-IR #
HOMAIR <- (meta$insulin)*(meta$glucose)/405
meta$HOMAIR <- HOMAIR
meta

# Filter high/low HOMA-IR #

meta$insulin_resistance <- ifelse(meta$HOMAIR > 2.5, "yes", "no")

sum(meta$insulin_resistance == "yes") # 195 individuals have HOMA-IR > 2.5
sum(meta$insulin_resistance == "no") # 246 individuals have HOMA-IR < 2.5


# Save edited metadata file
write.table(meta, file="metadata_edit.txt", sep="\t", quote=FALSE, row.names=FALSE)

