library(ggplot2)
library(dplyr)

tips <- read.csv("./tips_prelim.log", sep="\t")

human_columns <- grep("Hsapiens", colnames(tips), value=TRUE)
human_hbb <- grep("HBB", human_columns, value=TRUE)
human_hba <- grep("HBA", human_columns, value=TRUE)

human_differences <- sample_n(tips[human_hba], 10000, replace = TRUE) - sample_n(tips[human_hbb], 10000, replace = TRUE)
sorted_diff <- apply(human_differences, 2, sort)

# Consider the minimum credibility interval required for zero difference to not be included
minimum_interval <- data.frame(rep(-1, length(human_hba)), seq(1, length(human_hba)), row.names = human_hba)
colnames(minimum_interval) <- c("Min_Interval", "Site")

for(i in 0:100){
  interval <- (100 - i) * 0.01
  q <- (1-interval)/2
  quantile_lower <- apply(sorted_diff, 2, quantile, probs=q)
  quantile_upper <- apply(sorted_diff, 2, quantile, probs=1-q)
  
  zero_excluding_sites <- !(quantile_lower <= 0 & quantile_upper >= 0)
  site_labels <- human_hba[array(zero_excluding_sites)]
  for(site in site_labels){
    if(minimum_interval[site,"Min_Interval"] < 100 - i){
      minimum_interval[site,"Min_Interval"] <- 100 - i
    }
  }
}

ggplot(data = minimum_interval, mapping = aes(Site, Min_Interval)) + geom_col() +
  ylim(0, 100) + geom_hline(yintercept = 90, color="red", size=0.75, linetype=2) +
  theme_classic()
