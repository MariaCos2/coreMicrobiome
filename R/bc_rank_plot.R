##' @title bc_rank_plot
##'
##' @description Plots BC similarity vs ranked otu.
##' @param otu a community count data matrix,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_abun a value indicating whether the otu is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
##' @param sample_name the name of a column that containing the sample name,
##' @param sample_group the name of a column that containing the sample type/group,
##' @param ...
##' @details Plots BC similarity vs ranked otu.
##' @return A scatter plot showing the relationship between BC similarity and ranked otu with color denoting and core or non-core otu and a line indicating the threshold.
##' @examples
##'  bc_rank_plot(otu, sample, mini_abun=0, threshold=0.02, sample_name, sample_group)
##' @export

bc_rank_plot <- function(otu, sample, mini_abun, threshold, sample_name, sample_group){

  BC_ranked <- common_core(otu, sample, mini_abun = mini_abun, threshold = threshold, sample_name, sample_group)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun
  lastCall <- BC_ranked$lastCall

  p <- ggplot() +
    geom_point(data = BC_ranked_abun, aes(x = dim, y = proportionBC, size = otu_rel, fill = fill), pch = 21, color = 'black', alpha = .8) +
    #geom_point(data = BC_ranked_abun[BC_ranked_abun$fill == 'non-core',], aes(x = dim, y = proportionBC, size = otu_rel), pch = 21, fill = 'darkseagreen3', alpha = .8) +
    #geom_point(data = BC_ranked_abun[BC_ranked_abun$fill != 'non-core',], aes(x = dim, y = proportionBC, size = otu_rel), pch = 21, fill = 'deeppink4') +
    geom_vline(xintercept = lastCall, lty = 2, col = 'black', cex = 2) +
    labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
    scale_fill_manual(values = c("deeppink4", "darkseagreen3")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, color = "black", hjust = 0.5, vjust = 1, lineheight = 0.2),
          axis.title.x = element_text(size = 14, color = "black", hjust = 0.5),
          axis.title.y = element_text(size = 14,color = "black", hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 14),
          text = element_text(size =14),
          legend.position = "right",
          legend.title = element_text(colour = "black", size = 14),
          legend.text = element_text(colour ="black", size = 14),
          panel.background = element_rect(colour = "black", size = 1)
    ) + labs(size = "Relative abun", fill = "Common core")


  #ggplotly(p, height = height, width = weight, showlegend = T)
  p

}

