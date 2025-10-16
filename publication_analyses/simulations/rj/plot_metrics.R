library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

df <- read.csv("./posterior_stats.csv")
pivot <- df %>% pivot_longer(cols = starts_with("probs"), names_to = "model", values_to = "probs") %>%
  group_by(true_labels, model) %>% 
  summarize(mean = mean(probs), sd = sd(probs))

ggplot(pivot, aes(x = true_labels, y = mean, fill = model)) + 
  geom_col(position = position_dodge(width = 0.9)) + theme_classic() + 
  labs(fill = "Inferred Model") + xlab("True Model") + ylab("Mean Posterior Probability")

tip_df <- read.csv("./tips.2.log5", sep="\t") %>% select(-c(Iteration))
tip_sim <- read.table(text = readLines("./simulation.2.log5", warn=FALSE), sep="\t", header=TRUE)
sample_subset <- sample(colnames(tip_df), 20)
df_subset <- tip_df %>% select(sample_subset)
dnds_pivot <- df_subset %>% pivot_longer(cols= starts_with("Taxon"), names_to = "TaxonSite", values_to = "dNdS")
sim_tip_subset <- tip_sim %>% select(sample_subset)
sim_pivot <- sim_tip_subset %>% pivot_longer(cols= starts_with("Taxon"), names_to = "TaxonSite", values_to = "dNdS")

ggplot(dnds_pivot, aes(dNdS)) + geom_density(aes(fill=as.factor(TaxonSite)), alpha=0.15) +
  geom_rug(data=sim_pivot, aes(x=dNdS, color=TaxonSite), alpha=0.25) + theme_classic() +
  theme(legend.position = "none")

col_pairs <- combn(names(df_subset), 2, simplify=FALSE)
pairwise_est <- as.data.frame(
  t(apply(df_subset, 1, function(row){
    sapply(col_pairs, function(cols) abs(row[cols[1]] - row[cols[2]]))
  })) 
)
pairwise_pivot <- pairwise_est %>% pivot_longer(cols= starts_with("Taxon"), names_to = "TaxonSite", values_to = "dNdS_diff")


ggplot(pairwise_pivot, aes(dNdS_diff)) + geom_density(aes(fill=as.factor(TaxonSite)), alpha=0.15) + theme_classic() +
  theme(legend.position = "none")
