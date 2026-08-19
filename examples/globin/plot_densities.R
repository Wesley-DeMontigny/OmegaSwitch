library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

tip_df <- read.csv("./example_results/tips1.log", sep="\t") %>% select(-c(Iteration))

sample_subset <- sample(colnames(tip_df), 15)
df_subset <- tip_df %>% select(sample_subset)
dnds_pivot <- df_subset %>% pivot_longer(cols = everything(), names_to = "TaxonSite", values_to = "dNdS")

ggplot(dnds_pivot, aes(dNdS)) + geom_density(aes(group = TaxonSite), alpha=0.9) +
  theme_classic() +
  theme(legend.position = "none") + ylab("Density") + xlab("dN/dS")

