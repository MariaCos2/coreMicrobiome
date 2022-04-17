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
##' @param ttip the node size of the tree,
##' @param tlab the font size of the tree tip label,
##' @param hlab the font size of the heatmap label,
##' @param offset the offset between each subfigure,
##' @param width the width of each heatmap cell,
##' @param core_col the color of core tips,
##' @param noncore_col the color of non-core tips,
##' @param high_col the highest color of the heatmap,
##' @param low_col the lowest color of the heatmap,
##' @param ...
##' @details Plots phylogenetic tree with a heatmap denoting the average occurrence frequency of each taxa.
##' @return A phylogenetic tree with a heatmap denoting the average occurrence frequency of each taxa with color denoting and core or non-core otu.
##' @examples
##'  tree_abun_plot(otu, sample, tree, mini_abun=0, threshold=0.02, sample_name, sample_group, offset = 0.1, width = 0.5, core_col = 'deeppink4', noncore_col = 'darkseagreen3', low_col = 'white', high_col = 'lightcoral', hlab = 2, ttip = 2, tlab = 2)
##' @export

tree_abun_plot <- function(otu, sample, tree, mini_abun, threshold, sample_name, sample_group, offset, width, core_col, noncore_col, low_col, high_col, hlab, ttip, tlab){

  BC_ranked <- common_core(otu, sample, mini_abun = mini_abun, threshold = threshold, sample_name, sample_group)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun
  otu_ranked <- BC_ranked$otu_ranked

  annotation_row = data.frame(
    ID = BC_ranked_abun$otu,
    Type = BC_ranked_abun$fill
  )
  cols <- c('core' = core_col, 'non-core' = noncore_col)
  p <- ggtree(tree) +  geom_tiplab(size = tlab, align = TRUE)
  p1 <- p %<+% annotation_row +
    geom_tippoint(aes(color = Type), size = ttip) +
    scale_color_manual(values = cols)

  otu_relaocc_plot = get_avg_occ(otu, sample, otu_ranked, BC_ranked_abun, sample_name, sample_group)

  gheatmap(p1, otu_relaocc_plot, offset = offset, width = width, low = low_col, high = high_col, color = "black", colnames_angle = 90, hjust = 0, font.size = hlab, legend_title = "Occupancy frequency")

}

