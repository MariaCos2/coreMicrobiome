##' @title preve_abun_dis_plot
##'
##' @description Plots average occupancy and average relative abundance between core and non-core otu.
##' @param otu a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_abun a value indicating whether the otu is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
##' @param sample_name the name of a column that containing the sample name,
##' @param sample_group the name of a column that containing the sample type/group,
##' @param core_col the color of the core,
##' @param noncore_col the color of the non-core,
##' @param ...
##' @details Plots average occupancy and average relative abundance between core and non-core otu.
##' @return Two density plots showing the average occupancy and average relative abundance between core and non-core otu.
##' @examples
##'  preve_abun_dis_plot(otu, sample, mini_abun=0, threshold=0.02, sample_name, sample_group, core_col = "deeppink4", noncore_col = "darkseagreen3")
##' @export

preve_abun_dis_plot <- function(otu, sample, mini_abun, threshold, sample_name, sample_group, core_col, noncore_col){ # , height, weight

  otu <- data.frame(otu)
  otu_PA <- 1*((otu>mini_abun) == 1)                                               # presence: otu>0
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
  otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     # mean relative abundance
  occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance data frame
  occ_abun$rank <- as.factor(occ_abun$otu)

  BC_ranked <- common_core(otu, sample, mini_abun = mini_abun, threshold = threshold, sample_name, sample_group)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun

  occ_abun_noncore <- occ_abun %>% filter(!rank %in% BC_ranked_abun$otu[BC_ranked_abun$fill == 'common core'])
  occ_abun_noncore$type <- 'non-common core'
  occ_abun_core <- occ_abun %>% filter(rank %in% BC_ranked_abun$otu[BC_ranked_abun$fill == 'non-common core'])
  occ_abun_core$type <- 'common core'
  occ_abun_all <- rbind(occ_abun_noncore, occ_abun_core)

  p1 <- ggdensity(occ_abun_all, x = "otu_occ",
            add = "mean", rug = TRUE,
            color = "type", fill = "type",
            palette = c(core_col, noncore_col)) +
    labs(title = 'Avg occupancy frequency of taxa') +
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
            palette = c(core_col, noncore_col)) +
    labs(title = 'Avg relative abundance of taxa') +
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


  list(p1, p2)
  #p <- subplot(p1, p2, nrows = 2)

  #ggarrange(p1, p2, ncol = 2, nrow = 1)
  #ggplotly(p, height = height, width = weight, showlegend = T)
  #list(p1, p2)
}

