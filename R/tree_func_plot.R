##' @title tree_func_plot
##'
##' @description Plots phylogenetic tree with a heatmap denoting the average occurrence frequency and average tpm of the functions of each taxa.
##' @param functional_profile a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param tree a phylogenetic tree. The tip labels must match the otu names in the otu,
##' @param mini_tpm a value indicating whether the function is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
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
##' @details Plots phylogenetic tree with a heatmap denoting the average occurrence frequency and average tpm of the functions of each taxa.
##' @return A phylogenetic tree with a heatmap denoting the average occurrence frequency and average tpm of the functions of each taxa with color denoting and core or non-core functional otu.
##' @examples
##'  tree_func_plot(functional_profile, sample, tree, mini_tpm=1000, threshold=0.5, offset = 0.1, width = 0.5, core_col = 'deeppink4', noncore_col = 'darkseagreen3', low_col = 'white', high_col = 'lightcoral', hlab = 2, ttip = 2, tlab = 2, hangle = 270, hoffset = -1)
##' @export

tree_func_plot <- function(functional_profile, sample, tree, mini_tpm, threshold, offset, width, core_col, noncore_col, low_col, high_col, hlab, ttip, tlab, hangle, hoffset){

  fun_cores <- function_core(functional_profile, sample, mini_tpm, threshold)
  fun_tpm_sum <- fun_cores[[1]]
  fun_core <- fun_cores[[2]]

  fun_occ <- fun_core %>% dplyr::select('fun', 'fun_occ', 'genome')
  fun_occ <- tidyr::spread(fun_occ, key = "fun", value = "fun_occ")
  fun_occ[is.na(fun_occ)] <- 0
  fun_occ <- fun_occ[, colSums(fun_occ) != 0]

  fun_rel <- fun_core %>% dplyr::select('fun', 'fun_rel', 'genome')
  fun_rel <- tidyr::spread(fun_rel, key = "fun", value = "fun_rel")
  fun_rel[is.na(fun_rel)] <- 0
  fun_rel <- fun_rel[, colnames(fun_occ)]

  annotation_row = data.frame(
    ID = fun_tpm_sum$genome,
    Type = fun_tpm_sum$fill
  )
  cols <- c('functional core' = core_col, 'non-functional core' = noncore_col)
  p <- ggtree(tree) +  geom_tiplab(size = 3, align = TRUE)
  p1 <- p %<+% annotation_row +
    geom_tippoint(aes(color = Type), size = 3) +
    scale_color_manual(values = cols)

  fun_occ_df <- data.frame(fun_occ)
  row.names(fun_occ_df) <- fun_occ_df$genome
  fun_occ_df <- fun_occ_df[,-grep("genome",colnames(fun_occ_df))]
  h <- fun_occ_df
  pt1 <- gheatmap(p1, h, offset = offset, width = width, low = low_col, high = high_col, color = "black", colnames_angle = 90, hjust = 0, font.size = hlab, legend_title = "Occupancy frequency")

  fun_rel_df <- data.frame(fun_rel)
  row.names(fun_rel_df) <- fun_rel_df$genome
  fun_rel_df <- fun_rel_df[,-grep("genome",colnames(fun_rel_df))]
  h <- fun_rel_df
  pt2 <- gheatmap(p1, h, offset = offset, width = width, low = low_col, high = high_col, color = "black", colnames_angle = hangle, colnames_offset_y = hoffset, hjust = 0, font.size = hlab, legend_title = "Avg TPM")

  list(pt1, pt2)
}
