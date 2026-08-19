library(ape)
library(ggtree)
library(treeio)
library(dplyr)
library(stringr)
library(readr)
library(scales)
library(ggplot2)

newick_string <- "((((HBBDrerio[&index=20]:0.0901906,HBBCcarpio[&index=16]:0.15413)[&index=26]:0.616273,HBBSsalar[&index=24]:0.868368)[&index=27]:2.03548,((((HBBAindicus[&index=13]:0.021474,(HBBCminor[&index=17]:0.153898,HBBGgallus[&index=21]:0.224223)[&index=28]:0.0540914)[&index=29]:0.282224,(HBBCniloticus[&index=19]:2.01042,(HBBCmydas[&index=18]:0.304532,HBBPcastaneus[&index=23]:0.552562)[&index=30]:0.0166525)[&index=31]:0.0157695)[&index=32]:0.391196,(HBBBtaurus[&index=15]:0.540356,HBBHsapiens[&index=22]:0.32199)[&index=33]:1.45039)[&index=34]:0.461628,(HBBBbombina[&index=14]:0.99007,HBBXborealis[&index=25]:2.377)[&index=35]:1.54353)[&index=36]:0.0239444)[&index=37]:2.49336,(((HBADrerio[&index=7]:0.0390574,HBACcarpio[&index=3]:0.375545)[&index=38]:0.95593,HBASsalar[&index=11]:1.43851)[&index=39]:2.30223,((HBABbombina[&index=1]:0.89233,HBAXborealis[&index=12]:1.44808)[&index=40]:1.64901,((HBABtaurus[&index=2]:0.49234,HBAHsapiens[&index=9]:0.210743)[&index=41]:0.636702,((HBAAindicus[&index=0]:0.166592,(HBACminor[&index=4]:0.130561,HBAGgallus[&index=8]:0.1698)[&index=42]:0.0974355)[&index=43]:0.34458,((HBACmydas[&index=5]:0.514189,HBAPcastaneus[&index=10]:0.14252)[&index=44]:0.183971,HBACniloticus[&index=6]:0.847336)[&index=45]:0.179804)[&index=46]:0.268633)[&index=47]:0.733996)[&index=48]:0.600118)[&index=49]:1.17557)[&index=50]:0.0;"
trace_file <- "./example_results/ancestral_states1.log"
site_to_plot <- 35
burnin_frac <- 0

tree <- read.tree(text = newick_string)

# Make the tree ultrametric
tree <- chronos(tree, quiet = TRUE)

all_indices <- as.integer(
  str_match_all(newick_string, "\\[&index=([0-9]+)\\]")[[1]][,2]
)

node_map <- tibble(
  node = 1:(Ntip(tree) + tree$Nnode),
  trace_node_id = all_indices
)


trace <- read_tsv(trace_file, show_col_types = FALSE)

if ("Iteration" %in% names(trace)) {
  trace <- select(trace, -Iteration)
}

if (burnin_frac > 0) {
  trace <- trace[(floor(nrow(trace) * burnin_frac) + 1):nrow(trace), ]
}

pattern <- paste0("^([0-9]+)\\[", site_to_plot, "\\]$") 
cols <- names(trace)[str_detect(names(trace), pattern)] 

site_means <- tibble( 
  trace_node_id = as.integer(str_match(cols, pattern)[,2]), 
  value = colMeans(trace[cols], na.rm = TRUE) 
) 

node_values <- left_join(node_map, site_means, by = "trace_node_id") 

p <- ggtree(tree) %<+% 
  node_values + 
  geom_tree(aes(color = value), size = 1.2) + 
  geom_tiplab(offset = 0.1, size = 3) + 
  scale_color_gradient( 
    low = "yellow", 
    high = "red3", 
    limits = c(0, 2.0), 
    oob = squish, 
    name = paste0("Posterior mean\nsite ", site_to_plot) ) + 
  theme_tree2() 

p <- p + expand_limits(x = max(p$data$x) * 1.2) 

print(p)