# library(dplyr)
# library(ggpubr)
# library(plotly)
# library(gridExtra)
# source("functions/cores.R")

# There must be Sample_ID and SampleType in map file.

##' @title preve_abun_dis_plot
##' 
##' @description Plots average occupancy and average relative abundance between core and non-core otu. 
##' @param otu a community count data matrix,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_abun a value indicating whether the otu is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
##' @param sample_name the name of a column that containing the sample name,
##' @param sample_group the name of a column that containing the sample type/group,
##' @param ...
##' @details Plots average occupancy and average relative abundance between core and non-core otu. 
##' @return Two density plots showing theaverage occupancy and average relative abundance between core and non-core otu. 
##' @examples
##'  preve_abun_dis_plot(otu, sample, mini_abun=0, threshold=0.02, height=400, weight=400)
##' @export

preve_abun_dis_plot <- function(otu, sample, mini_abun, threshold, sample_name, sample_group){ # , height, weight
  
  # otu = read.csv(otu, sep = ",", header = T, row.names = 1)
  # sample = read.csv(sample, sep = ",", header = T, row.names = 1)
  
  otu_PA <- 1*((otu>mini_abun) == 1)                                               # presence: otu>0
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
  otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     # mean relative abundance
  occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance data frame
  occ_abun$rank <- as.factor(occ_abun$otu)
  
  BC_ranked <- common_core(otu, sample, mini_abun = mini_abun, threshold = threshold, sample_name, sample_group)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun
  
  occ_abun_noncore <- occ_abun %>% filter(!rank %in% BC_ranked_abun$otu[BC_ranked_abun$fill == 'core'])
  occ_abun_noncore$type <- 'non-core'
  occ_abun_core <- occ_abun %>% filter(rank %in% BC_ranked_abun$otu[BC_ranked_abun$fill == 'core'])
  occ_abun_core$type <- 'core'
  occ_abun_all <- rbind(occ_abun_noncore, occ_abun_core)
  
  p1 <- ggdensity(occ_abun_all, x = "otu_occ",
            add = "mean", rug = TRUE,
            color = "type", fill = "type",
            palette = c("deeppink4", "darkseagreen3")) +
    labs(title = 'Avg occupancy of otu') +
    xlab('avg occupancy') + 
    ylab('') +
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
    ) 
  
  p2 <- ggdensity(occ_abun_all, x = "otu_rel",
            add = "mean", rug = TRUE,
            color = "type", fill = "type",
            palette = c("deeppink4", "darkseagreen3")) +
    labs(title = 'Avg relative abundance of otu') +
    xlab('avg relative abundance') + 
    ylab('') +
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
    ) 
  
  
  #p <- subplot(p1, p2, nrows = 2)
  
  # ggarrange(p1, p2,
  #           ncol = 1, nrow = 2)
  #ggplotly(p, height = height, width = weight, showlegend = T) 
  list(p1, p2)
}

