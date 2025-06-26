library(ggplot2)
library(dplyr)

tips <- read.csv("./tips1.log", sep="\t")

hbb <- grep("HBB", colnames(tips), value=TRUE)
hba <- grep("HBA", colnames(tips), value=TRUE)

row_samples <- sample(nrow(tips), 10000, replace=TRUE)

hba_samples <- tips[row_samples, hba]
colnames(hba_samples)<-gsub("HBA", "", hba)
hbb_samples <- tips[row_samples, hbb]
colnames(hbb_samples)<-gsub("HBB", "", hbb)

differences <- hba_samples - hbb_samples[, colnames(hba_samples)]
suffixes <- sapply(strsplit(colnames(differences), "\\."), tail, 1)
grouped_cols <- split(colnames(differences), suffixes)
site_difference_means <- sapply(grouped_cols, function(cols) {
  rowMeans(differences[, cols, drop = FALSE])
})
sorted_diff <- apply(site_difference_means, 2, sort)

# Consider the minimum credibility interval required for zero difference to not be included
minimum_interval <- data.frame(rep(-1, ncol(site_difference_means)), seq(0, ncol(site_difference_means)-1))
colnames(minimum_interval) <- c("Min_Interval", "Site")
minimum_interval$Magnitude <- abs(colMeans(sorted_diff))
minimum_interval$Diff <- colMeans(sorted_diff)

for(i in 0:100){
  interval <- (100 - i) * 0.01
  q <- (1-interval)/2
  quantile_lower <- apply(sorted_diff, 2, quantile, probs=q)
  quantile_upper <- apply(sorted_diff, 2, quantile, probs=1-q)
  
  zero_excluding_sites <- !(quantile_lower <= 0 & quantile_upper >= 0)
  for(site in names(zero_excluding_sites)){
    print(site)
    if(zero_excluding_sites[site]){
      if(minimum_interval[minimum_interval$Site == site,"Min_Interval"] < 100 - i){
        minimum_interval[minimum_interval$Site == site,"Min_Interval"] <- 100 - i
      }
    }
  }
}

ggplot(data = minimum_interval, mapping = aes(Min_Interval, Magnitude)) + geom_point() +
  xlim(0, 100) + geom_vline(xintercept = 95, color="red", size=0.75, linetype=2) +
  theme_classic() + geom_hline(yintercept = 0.5, color="red", size=0.75, linetype=2) + 
  ylab("Absolute Mean Difference in dN/dS") + xlab("Minimum Credibility Interval Excluding Zero") +
  annotate("rect", xmin = -Inf, xmax = 95, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
  annotate("rect", xmin=95, xmax=Inf, ymin=-Inf, ymax=0.5, fill = "red", alpha = 0.1)

print(minimum_interval[minimum_interval$Min_Interval >= 95 & minimum_interval$Magnitude > 0.5,])
