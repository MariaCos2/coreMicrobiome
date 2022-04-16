# library(dplyr)
# source("functions/cores.R")

# There must be Sample_ID and SampleType in map file.

##' @title abun_occ_plot
##' 
##' @description Plots occupancy vs abunbance.
##' @param otu a community count data matrix,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_abun a value indicating whether the otu is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
##' @param height the height of the figure,
##' @param weight the width of the figure,
##' @param ...
##' @details Plots occupancy vs abunbance.
##' @return A scatter plot showing the relationship between occupancy and abundance with color denoting and core or non-core otu.
##' @examples
##'  abun_occ_plot(otu, sample, mini_abun=0, threshold=0.02, height=400, weight=400)
##' @export

abun_occ_plot <- function(otu, sample, mini_abun, threshold, height, weight){
  
  otu = read.csv(otu, sep = ",", header = T, row.names = 1)
  sample = read.csv(sample, sep = ",", header = T, row.names = 1)
  
  BC_ranked <- common_core(otu, sample, mini_abun = mini_abun, threshold = threshold)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun
  
  p <- ggplot() +
    geom_point(data=BC_ranked_abun, aes(x = log10(otu_rel), y = otu_occ, size = otu_rel, fill = fill), pch = 21, color = 'black', alpha = .8)+
    #geom_point(data=BC_ranked_abun[BC_ranked_abun$fill == 'non-core',], aes(x = log10(otu_rel), y = otu_occ, size = otu_rel), pch = 21, fill = 'darkseagreen3', alpha = .8)+
    #geom_point(data=BC_ranked_abun[BC_ranked_abun$fill != 'non-core',], aes(x = log10(otu_rel), y = otu_occ, size = otu_rel), pch = 21, fill = 'deeppink4') +
    labs(x= "log10(mean relative abundance)", y = "Occupancy") +
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
  p
  #ggplotly(p, height = height, width = weight, showlegend = T) 
}

