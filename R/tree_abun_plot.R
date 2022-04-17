# library(ape)
# library(ggtree)
# library(tidyr)
# source("functions/cores.R")

# There must be Sample_ID and SampleType in map file.

##' @title tree_abun_plot
##' 
##' @description Plots phylogenetic tree with a heatmap denoting the average occurrence frequency of each taxa.
##' @param otu a community count data matrix,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param tree a phylogenetic tree. The tip labels must match the otu names in the otu,
##' @param mini_abun a value indicating whether the otu is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
##' @param sample_name the name of a column that containing the sample name,
##' @param sample_group the name of a column that containing the sample type/group,
##' @param ...
##' @details Plots phylogenetic tree with a heatmap denoting the average occurrence frequency of each taxa.
##' @return A phylogenetic tree with a heatmap denoting the average occurrence frequency of each taxa with color denoting and core or non-core otu.
##' @examples
##'  tree_abun_plot(otu, sample, tree, mini_abun=0, threshold=0.02)
##' @export

tree_abun_plot <- function(otu, sample, tree, mini_abun, threshold, sample_name, sample_group){
  
  # otu = read.csv(otu, sep = ",", header = T, row.names = 1)
  # sample = read.csv(sample, sep = ",", header = T, row.names = 1)
  # tree = read.tree(tree)
  
  BC_ranked <- common_core(otu, sample, mini_abun = mini_abun, threshold = threshold, sample_name, sample_group)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun
  otu_ranked <- BC_ranked$otu_ranked
  
  annotation_row = data.frame(
    ID = BC_ranked_abun$otu,
    Type = BC_ranked_abun$fill
  )
  cols <- c('core' = 'deeppink4', 'non-core' = 'darkseagreen3')
  p <- ggtree(tree) +  geom_tiplab(size = 2, align = TRUE) 
  p1 <- p %<+% annotation_row + 
    geom_tippoint(aes(color = Type), size = 1.5) + 
    scale_color_manual(values = cols) 
  
  otu_relaocc_plot = get_avg_occ(otu, sample, otu_ranked, BC_ranked_abun, sample_name, sample_group)
  
  gheatmap(p1, otu_relaocc_plot, offset = 0.1, width = 0.5, low = 'white', high = 'lightcoral',color = "black", colnames_angle = 270, hjust = 0, font.size = 4, legend_title = "Occupancy frequency")
  
}

